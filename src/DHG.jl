module DHG

include("./Utilities.jl")               # submodule Utilities
include("./FVM.jl")                     # provides submodule FVM
include("./NetworkComponents.jl")       # submodule NetworkComponents
include("./Transport.jl")               # submodule TransportProperties
include("./DynamicalFunctions.jl")      # submodule DynamicalFunctions
include("./WrapperFunctions.jl")        # submodule WrapperFunctions
include("./GraphParsing.jl")            # submodule GraphParsing
include("./ParameterStruct.jl")         # submodule ParameterStructs

# TODO: global size_t, but NetworkDynamics takes dims::Int

# import Graphs: SimpleGraphs.SimpleDiGraph
# import NetworkDynamics
# import DifferentialEquations: SciMLBase

# import and re-export from submodules for easier access
import .GraphParsing: parse_gml
import .ParameterStruct: Parameters
import .Transport: TransportProperties

export Parameters, TransportProperties, parse_gml
export Prosumer, Prosumer_PressureChange, Prosumer_Massflow


abstract type Prosumer end
struct Prosumer_PressureChange <: Prosumer end
struct Prosumer_Massflow <: Prosumer end


# export DHGStruct


# Base.@kwdef struct DHGStruct{IndexType <: Integer, FunctionType <: SciMLBase.AbstractODEFunction}
#     f::FunctionType
#     parameters::ParameterStructs.Parameters
#     graph::SimpleDiGraph{IndexType}
#     n_states::NamedTuple{(:nodes, :edges),
#                          Tuple{Vector{UInt16}, Vector{UInt16}} # Determines max number of states, change in constructor as well if modified
#                         }
# end


# function DHGStruct(graph_parser::Function,
#                     global_parameters::ParameterStructs.GlobalParameters,
#                     transport_coeffs::TransportProperties.TransportCoefficients,
#                 )
#     # Constructor

#     graph, node_dict, edge_dict = graph_parser()
#     parameters = ParameterStructs.Parameters(global_parameters = global_parameters,
#                                              node_parameters = ParameterStructs.NodeParameters(node_dict),
#                                              edge_parameters = ParameterStructs.EdgeParameters(edge_dict)
#                                             )
#     # Assemble dynamical functions
#     node_fns = Vector{NetworkDynamics.DirectedODEVertex}(undef, length(node_dict.components))
#     edge_fns = Vector{NetworkDynamics.ODEEdge}(undef, length(edge_dict.components))
#     n_states = (nodes=similar(node_fns, UInt16), edges=similar(edge_fns, UInt16))

#     for (i, node) in enumerate(node_dict.components)
#         node_fns[i] = WrapperFunctions.node_fn(node)
#         n_states.nodes[i] = node_fns[i].dim
#     end

#     for (i, edge) in enumerate(edge_dict.components)
#         edge_fns[i] = WrapperFunctions.edge_fn(i, edge, transport_coeffs)
#         n_states.edges[i] = edge_fns[i].dim
#     end

#     ode_fn = NetworkDynamics.network_dynamics(node_fns, edge_fns, graph)

#     return DHGStruct(f=ode_fn, parameters=parameters, graph=graph, n_states=n_states)
# end

end # module DHG
