module DHG

include("./Utilities.jl")               # submodule Utilities
include("./FVM.jl")                     # provides submodule FVM
include("./NetworkComponents.jl")       # submodule NetworkComponents
include("./TransportProperties.jl")     # submodule TransportProperties
include("./DynamicalFunctions.jl")      # submodule DynamicalFunctions
include("./WrapperFunctions.jl")        # submodule WrapperFunctions
include("./GraphParsing.jl")            # submodule GraphParsing
include("./ParameterStructs.jl")        # submodule ParameterStructs

# TODO: export parse_gml, others?
# TODO: global size_t, but NetworkDynamics takes dims::Int

import Graphs: SimpleGraphs.SimpleDiGraph
import NetworkDynamics
import DifferentialEquations: SciMLBase

# import and re-export from submodules for easier access
import .GraphParsing: parse_gml
import .ParameterStructs: GlobalParameters
import .TransportProperties: TransportCoefficients

export GlobalParameters, TransportProperties, parse_gml

export DHGStruct


Base.@kwdef struct DHGStruct{IndexType <: Integer}
    node_functions::Vector{NetworkDynamics.DirectedODEVertex}
    edge_functions::Vector{NetworkDynamics.ODEEdge}
    parameters::ParameterStructs.Parameters
    graph::SimpleDiGraph{IndexType}
    edges::Vector{NetworkComponents.Edge}
end


function DHGStruct(graph_parser::Function,
                    global_parameters::ParameterStructs.GlobalParameters,
                    transport_coeffs::TransportProperties.TransportCoefficients,
                )
    # Constructor

    graph, node_dict, edge_dict = graph_parser()
    edgevec::Vector{NetworkComponents.Edge} = edge_dict.components
    parameters = ParameterStructs.Parameters(global_parameters = global_parameters,
                                             node_parameters = ParameterStructs.NodeParameters(node_dict),
                                             edge_parameters = ParameterStructs.EdgeParameters(edge_dict)
                                            )
    # Assemble dynamical functions
    node_fns = Vector{NetworkDynamics.DirectedODEVertex}(undef, length(node_dict.components))
    edge_fns = Vector{NetworkDynamics.ODEEdge}(undef, length(edgevec))

    for (i, node) in enumerate(node_dict.components)
        node_fns[i] = WrapperFunctions.node_fn(node)
    end

    for (i, edge) in enumerate(edgevec)
        edge_fns[i] = WrapperFunctions.edge_fn(i, edge, transport_coeffs)
    end

    return DHGStruct(node_functions=node_fns, edge_functions=edge_fns, parameters=parameters, graph=graph, edges=edgevec)
end

end # module DHG
