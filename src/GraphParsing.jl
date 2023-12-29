module GraphParsing

import ParserCombinator
import Graphs
import ArgCheck: @argcheck
import DataStructures as ds

# TODO: exports


function parse_gml(filename::AbstractString)
    # Parse a graph stored as a text file in GML format:
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

    nodes_vec = parser_dict[:graph][1][:node]
    edges_vec = parser_dict[:graph][1][:edge]

    # Check for unsupported/incompatible inputs
    @argcheck length(parser_dict[:graph]) == 1 "GML file defines multiple graphs,
    only a single graph definition is supported"

    @argcheck parser_dict[:graph][:directed] == 1 "Graph must be directed,
    ie. 'directed 1' attribute must be present in graph [ ... ] scope in GML file"

    # All nodes and edges must have a 'type' attribute, allowed values for type:
    #   Nodes: junction, fixed
    #   Edges: pipe, prosumer

    edgechecks = ds.DefaultDict((_) -> throw(ArgumentError("unrecognised 'type' attribute")), # default value
                                Dict("pipe" =>  _checkparams_pipeedge,
                                     "prosumer" => _checkparams_prosumeredge,
                                    )
                               )
    
    # Check nodes
    num_fixednodes::Int = 0 # There can only be a single 'fixed' node in the network
    # Not using dict like edgechecks for nodes to increment num_fixednodes
    for (i, node) in enumerate(nodes_vec)
        try # Handle missing 'type' attribute or unexpected parameters for specified 'type'
            if node[:type] == "junction"
                nothing # No parameters required for junction nodes
            elseif node[:type] == "fixed"
                if num_fixednodes > 0; throw(ArgumentError("multiple fixed nodes specified")); end
                _checkparams_fixednode(node)
                num_fixednodes += 1
            else
                throw(ArgumentError("unrecognised 'type' attribute"))
            end
        catch exc # Insert node index into stacktrace
            if typeof(exc) <: Union{KeyError, ArgumentError} # KeyError if no 'type' attribute, _checkparams functions may throw ArgumentErrors
                throw(ArgumentError("GML node $i specification")) # Stacktrace should also show caught exception
            else
                rethrow(exc)
            end
        end
    end # loop

    # Check edges
    # Using edgechecks dict for cleaner control flow
    for (i, edge) in enumerate(edges_vec)
        try
            edgechecks[edge[:type]](edge)
        catch exc
            if typeof(exc) <: Union{KeyError, ArgumentError}
                throw(ArgumentError("GML edge $i specification"))
            else
                rethrow(exc)
            end
        end
    end






end


end # module GraphParsing
