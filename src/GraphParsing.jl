module GraphParsing

import ParserCombinator
import Graphs
import DataStructures as ds
import StaticArrays as sa

# TODO: exports
export parse_gml


# Custom exception type, only to be thrown by functions inside this module
struct GraphParsingError <: Exception
    msg::String
end
## Method to pretty-print
Base.showerror(io::IO, e::GraphParsingError) = print(io, "GraphParsingError: ", e.msg)


# "private" functions to check if required attributes and parameters are specified in the input network

function _checkattrs_common(component, common_symbols) :: Nothing
    # Check attributes common to all nodes and edges
    # Try to access required attributes, throw KeyErrors if these attributes are not present
    for symbol in common_symbols
        _ = component[symbol]
    end
    return nothing
end


function _checkattrs_allnodes(node) :: Nothing
    # Check attributes common to all nodes
    common_syms = sa.@SVector [:id]
    _checkattrs_common(node, common_syms) # TODO: @static useful here?
    return nothing
end


function _checkattrs_alledges(edge) :: Nothing
    # Check attributes common to all edges
    common_syms = sa.@SVector [:src, :dst]
    _checkattrs_common(edge, common_syms)
    return nothing
end


## Intended to be called by parse_gml()

function _checkparams_fixednode(node) :: Nothing
    # Try to access required attributes, throw KeyErrors if these attributes are not present
    # KeyErrors will be caught by caller parse_gml()
    _ = node[:fixed][:pressure]
    _ = node[:fixed][:temperature]
    return nothing
end


function _checkparams_prosumeredge(edge) :: Nothing
    # Try to access required attributes, throw KeyErrors if these attributes are not present
    # KeyErrors will be caught by caller parse_gml()
    k = keys(edge[:prosumer])
    _ = edge[:prosumer][:delta_T]
    if !(xor((:delta_p in k), (:massflow in k)))
        throw(GraphParsingError("prosumer edge must have either delta_p or massflow specified, but not both together"))
    end
    return nothing
end


function _checkparams_pipeedge(edge) :: Nothing
    # Try to access required attributes, throw KeyErrors if these attributes are not present
    # KeyErrors will be caught by caller parse_gml()
    for sym in sa.@SVector [:diameter, :length, :dx]
        if !(edge[:pipe][sym] > 0)
            throw(GraphParsingError("$(sym) must be positive"))
        end
    end
    return nothing
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
    #   Tuple(Graphs.jl graph, Dict(node id => AttributeDict), Dict(edge key => AttributeDict))
    #       where AttributeDict >: Dict(Symbol => attribute value type)

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

    if parser_dict[:graph][:directed] != 1
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






end


end # module GraphParsing
