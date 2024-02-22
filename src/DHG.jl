module DHG

include("Discretisation.jl")
include("Fluids.jl")
include("Transport.jl")
include("Miscellaneous.jl")
include("NetworkComponents.jl")
include("DynamicalFunctions.jl")
include("WrapperFunctions.jl")

# Imports from submodules to top-level namespace, for easier access from scripts
import .NetworkComponents: Prosumer, JunctionNode, ReferenceNode, Pipe, PressureChange, Massflow
import .Transport: TransportModels

export assemble

import Graphs
import NetworkDynamics


Sequence{T} = Union{Tuple{Vararg{T}}, Vector{T}}
    # type alias for convenience, intended use:
    # Vector{<:T} or Tuple{<:T, <:T, ..., <:T} where T is some abstract type


function assemble(nodes::Sequence{<:NetworkComponents.Node},
                  edges::Sequence{<:NetworkComponents.Edge},
                  transport_models::Transport.TransportModels,
                  discretisation::Discretisation.DiscretisationScheme,
                  ::Type{fluid}
                 ) where {fluid <: Fluids.Fluid}
    # -> Tuple{T, Graphs.SimpleDiGraph} where T <: DifferentialEquations.SciMLBase.ODEFunction
    # Wrapper around NetworkDynamics.network_dynamics() call,
    # passes through return from network_dynamics().
    # Creates Graphs.SimpleDiGraph from parameters nodes, edges

    # Check if edges are sorted according to Graphs.jl ordering
    ## Example of correctly sorted edges: [(1 => 3), (2 => 1)] where (e.src => e.dst) for e in edges
    ## Incorrectly sorted edges: [(1 => 3), (1 => 2)] or [(2 => 1), (1 => 3)]
    idx = Miscellaneous.adjacent_find((e1, e2) -> ( (e1.src > e2.src) ||
                                                    ((e1.src == e2.src) && (e1.dst > e2.dst))
                                                  ),
                                      edges
                                     )
    if idx != length(edges)
        throw("'edges' not sorted at indices $(idx), $(idx + 1): "*
              "Elements of 'edges' must be sorted in lexicographic order of (src, dst) for "*
              "uniformity with Graphs.jl data structure"
             )
    end

    # Create graph
    num_nodes = length(nodes)
    graph = Graphs.SimpleDiGraph(num_nodes)
    for e in edges
        # gr.add_edge!() -> bool to indicate success
        Graphs.add_edge!(graph, e.src, e.dst) ? nothing : throw("Failed to add edge "*
                                                                "$(e.src) => $(e.dst)"
                                                               )
    end

    # Call network_dynamics()
    nodes_gen = (WrapperFunctions.node(x) for x in nodes) # Generator
    edges_gen = (WrapperFunctions.edge(x, transport_models, discretisation, fluid)
                 for x in edges)

    # network_dynamics() only takes Vectors for nodes and edges
    return (NetworkDynamics.network_dynamics(collect(nodes_gen), collect(edges_gen), graph),
            graph
           ) # Tuple
end


end
