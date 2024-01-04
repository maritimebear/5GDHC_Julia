module GraphParsing # submodule, included in DHG.jl

export parse_gml, ComponentDict


import ParserCombinator
import Graphs
import DataStructures as ds
import StaticArrays as sa

import ..Utilities as utils
import ..NetworkComponents as nc



# Custom exception type, only to be thrown by functions inside this module
struct GraphParsingError <: Exception
    msg::String
end
## Method to pretty-print
Base.showerror(io::IO, e::GraphParsingError) = print(io, "GraphParsingError: ", e.msg)


# Custom dictionary type, to be returned from parse_gml()
Base.@kwdef struct ComponentDict{IdxType, # Vector index type
                                 ComponentType <: Union{nc.Node, nc.Edge}}
    # Holds a vector of network components and dictionaries to look up the indices of each component type
    # Usage example: _, _, edge_dict::ComponentDict{Int64, nc.Edge} = graph_parse()
    # edge_dict.indices[:prosumer] => Vector{Int64} with indices of prosumer edges,
    # edge_dict.components will correspond to these indices.

    indices::Dict{Symbol, Vector{IdxType}}
    components::Vector{ComponentType}
end


# "private" functions to check if required attributes and parameters are specified in the input network

# TODO: If argtype check doesn't work, replace with _ = component[symbol]
function _check_attribute(datastructure, key, ::Type{ExpectedType}) :: Nothing where {ExpectedType}
    # Check if an attribute is accessible by 'key' and is a subtype of 'ExpectedType'
    # Throws KeyError if access fails, TypeError if type assertion fails
    datastructure[key]::ExpectedType # Inline type assertion
    return nothing
end


function _checkattrs_allnodes(node) :: Nothing
    # Check attributes common to all nodes
    common_attrs = sa.@SVector [(:id, Int)]
    for (symbol, type) in common_attrs
        _check_attribute(node, symbol, type)
    end
    return nothing
end


function _checkattrs_alledges(edge) :: Nothing
    # Check attributes common to all edges
    common_attrs = sa.@SVector [(:source, Int), (:target, Int)]
    for (symbol, type) in common_attrs
        _check_attribute(edge, symbol, type)
    end
    return nothing
end


## Intended to be called by parse_gml()

function _checkparams_fixednode(node) :: Nothing
    # Parameters required for fixed nodes
    attrs = sa.@SVector [(:pressure, Real), (:temperature, Real)]
    _check_attribute(node, :fixed, AbstractDict)
    for (symbol, type) in attrs
        _check_attribute(node[:fixed], symbol, type)
    end
    return nothing
end


function _checkparams_prosumeredge(edge) :: Nothing
    # Parameters required for prosumer edges
    _check_attribute(edge, :prosumer, AbstractDict)
    _check_attribute(edge[:prosumer], :deltaT, Real)
    k = keys(edge[:prosumer])
    if !(xor((:deltaP in k), (:massflow in k)))
        throw(GraphParsingError("prosumer edge must have either deltaP or massflow specified, but not both together"))
    end
    if :deltaP in k
        _check_attribute(edge[:prosumer], :deltaP, Real)
    else
        _check_attribute(edge[:prosumer], :massflow, Real)
    end
    return nothing
end


function _checkparams_pipeedge(edge) :: Nothing
    # Parameters required for pipe edges
    _check_attribute(edge, :pipe, AbstractDict)
    for (symbol, type) in sa.@SVector [(:diameter, Real), (:length, Real), (:dx, Real)]
        _check_attribute(edge[:pipe], symbol, type)
        if !(edge[:pipe][symbol] > 0)
            throw(GraphParsingError("$(symbol) must be positive"))
        end
    end
    return nothing
end


function _construct_graph(nodes_vec, edges_vec, ::Type{IndexType}) where {IndexType <: Integer}
    # Returns Graphs.jl SimpleDiGraph, constructed from nodes_vec and edges_vec
    # !!! nodes_vec must be sorted in ascending order of node id !!!
    # This is required for the edge ordering to be correct, as Graphs.jl seems to implement
    # graphs based on just the number of nodes.
    # Based on GraphIO.jl, also see for unordered nodes implementation:
    #   _gml_read_one_graph(), https://github.com/JuliaGraphs/GraphIO.jl/blob/master/ext/GraphIOGMLExt.jl

    graph = Graphs.SimpleDiGraph{IndexType}(length(nodes_vec))
    for edge in edges_vec
        if !(Graphs.add_edge!(graph, edge[:source], edge[:target])) # add_edge! -> true/false to indicate success
            throw(GraphParsingError("Graphs.add_edge! failure: edge $(edge[:source]) => $(edge[:target])"))
        end
    end
    return graph
end


# "public" functions, to be exported

function parse_gml(filename::AbstractString)
    # Parse a network stored as a text file in GML format:
    #   https://en.wikipedia.org/wiki/Graph_Modelling_Language
    # GML supports specification of attributes for nodes and edges, this is used to specify
    # the types and parameters of network components.
    #   eg: junction node, prosumer edge with fixed massflow = 1.0
    #
    # Parameters:
    #   filename: text file to be parsed
    # Return:
    #   Tuple(Graphs.jl graph, ComponentDict{Integer, nc.Node}, ComponentDict{Integer, nc.Edge}

    parser_dict = ParserCombinator.Parsers.GML.parse_dict(read(filename, String)) # Read entire contents of filename as String |> parse_dict()

    # Each attribute in the GML hierarchy is parsed as a symbol, and is mapped to a Vector, String or number
    #   eg: parser_dict[:graph] is a Vector,
    #       parser_dict[:graph][1][:node] => Vector{'node' objects}, each object is a Dict{Symbol, value}

    nodes_vec = parser_dict[:graph][1][:node]
    edges_vec = parser_dict[:graph][1][:edge]

    # Check for unsupported/incompatible inputs
    if length(parser_dict[:graph]) != 1
        throw(GraphParsingError("GML file defines multiple graphs, only a single graph definition is supported"))
    end

    if parser_dict[:graph][1][:directed] != 1
        throw(GraphParsingError("Graph must be directed, ie. 'directed 1' attribute must be present in graph [ ... ] scope in GML file"))
    end

    # All nodes and edges must have a 'type' attribute, allowed values for type:
    #   Nodes: junction, fixed
    #   Edges: pipe, prosumer

    edgechecks = ds.DefaultDict((_) -> throw(GraphParsingError("unrecognised 'type' attribute")), # default value
                                Dict("pipe" =>  _checkparams_pipeedge,
                                     "prosumer" => _checkparams_prosumeredge,
                                    )
                               )
    
    # Check nodes
    fixednode_found::Bool = false # There can only be a single 'fixed' node in the network
    # Not using dict like edgechecks for nodes in order to mutate fixednode_found
    for (i, node) in enumerate(nodes_vec)
        try # Handle missing 'type' attribute or unexpected parameters for specified 'type'
            _checkattrs_allnodes(node) # Check that attributes common to all nodes are present
            if node[:type] == "junction"
                nothing # No parameters required for junction nodes
            elseif node[:type] == "fixed"
                if fixednode_found; throw(GraphParsingError("multiple fixed nodes specified")); end
                _checkparams_fixednode(node)
                fixednode_found = true
            else
                throw(GraphParsingError("unrecognised 'type' attribute"))
            end
        catch exc # Insert node index into stacktrace
            if typeof(exc) <: Union{KeyError, GraphParsingError} # KeyError if no 'type' attribute, _checkparams functions may throw GraphParsingErrors
                throw(GraphParsingError("GML node $i specification")) # Stacktrace should also show caught exception
            else
                rethrow(exc)
            end
        end
    end # loop

    # Check edges
    # Using edgechecks dict for cleaner control flow
    for (i, edge) in enumerate(edges_vec)
        try
            _checkattrs_alledges(edge)
            edgechecks[edge[:type]](edge)
        catch exc
            if typeof(exc) <: Union{KeyError, GraphParsingError}
                throw(GraphParsingError("GML edge $i specification"))
            else
                rethrow(exc)
            end
        end
    end

    # Sort nodes_vec in ascending order of node id, to maintain consistency with Graphs.jl ordering
    sort!(nodes_vec; lt = (lhs, rhs) -> (lhs[:id] < rhs[:id]))

    # Check if duplicate nodes are present (ie. multiple nodes with the same id)
    duplicate_idx = utils.adjacent_find((lhs, rhs) -> (lhs[:id] == rhs[:id]), nodes_vec)

    if duplicate_idx != length(nodes_vec) # adjacent_find() works similar to C++ std::adjacent_find()
        throw(GraphParsingError("multiple nodes with same id: $(duplicate_idx)"))
    end

    # At this point, the nodes are sorted in ascending order of node id, with no duplicate nodes.
    # This is required, as Graphs.jl seems to implement graphs based on just the number of nodes.
    # nodes_vec being sorted in ascending order guarantees that its ordering will match the
    # ordering of nodes in the graph structure from Graphs.jl.

    # The edge ordering is left up to Graphs.jl, and will be read later from the constructed
    # graph (the edge ordering scheme used by Graphs.jl does not appear to be documented, possibly
    # a sparse structure?)

    IndexType = Int # TODO: Change for optimisation?
    graph = _construct_graph(nodes_vec, edges_vec, IndexType)

    # Create Dicts and Vectors to be returned
    nodes_dict = ComponentDict(Dict(sym => Vector{IndexType}(undef, 0) for sym in [:junction, :fixed]),
                               Vector{nc.Node}(undef, 0)
                              )
    edges_dict = ComponentDict(Dict(sym => Vector{IndexType}(undef, 0) for sym in [:pipe, :deltaP, :massflow]),
                               Vector{nc.Edge}(undef, 0)
                              )

    for (i, node) in enumerate(nodes_vec)
        if node[:type] == "junction"
            push!(nodes_dict.indices[:junction], i)
            push!(nodes_dict.components,
                  nc.Junction())

        else # node[:type] == "fixed" guaranteed after previous checks
            push!(nodes_dict.indices[:fixed], i)
            push!(nodes_dict.components,
                  nc.FixedNode(node[:fixed][:pressure],
                               node[:fixed][:temperature])
                 )
        end
    end

    for (i, edge) in enumerate(edges_vec)
        if edge[:type] == "pipe"
            push!(edges_dict.indices[:pipe], i)
            push!(edges_dict.components,
                  nc.Pipe(length=edge[:pipe][:length],
                          diameter=edge[:pipe][:diameter],
                          dx=edge[:pipe][:dx])
                 )

        else # edge[:type] == "prosumer" guaranteed after previous checks
            if :deltaP in keys(edge[:prosumer])
                push!(edges_dict.indices[:deltaP], i)
                push!(edges_dict.components,
                      nc.PressureChangeProsumer(deltaP=edge[:prosumer][:deltaP],
                                                deltaT=edge[:prosumer][:deltaT])
                     )
            else # :massflow in keys(edge[:prosumer]) guaranteed after previous checks
                push!(edges_dict.indices[:massflow], i)
                push!(edges_dict.components,
                      nc.MassflowProsumer(massflow=edge[:prosumer][:massflow],
                                          deltaT=edge[:prosumer][:deltaT])
                     )
            end
        end
    end

    return (graph, nodes_dict, edges_dict)
end


end # (sub)module
