module DHG

import NetworkDynamics as nd
import LinearAlgebra as la

include("./Physics.jl") # provides module Physics

import .Physics

export pipe_edge, prosumer_edge, junction_node, fixed_node

## Wrapper functions around NetworkDynamics structs to calculate arguments to their constructors

function pipe_edge(length::Real, dx::Real) :: nd.ODEEdge
    # Edge with friction-induced pressure drop, temperature loss to surroundings
    # state 1 => mass flow rate
    # states 2:end => temperature in finite-volume cells
    # dims == num(states) = 1 + num(cells)

    n_cells = length / dx # finite-volume cells
    @assert (n_cells > 0 && n_cells != Inf)

    # Calculate arguments to NetworkDynamics.ODEEdge
    dims = n_cells + 1
    diagonal = la.Diagonal([1 for _ in 1:dims]) # Diagonal of mass matrix
    diagonal[1] = 0 # state 1 corresponds to an algebraic constraint
    symbols = [Symbol("T$i") for i in 0 : n_cells]
    symbols[1] = :m

    return nd.ODEEdge(f=Physics.pipe!, dim=dims, coupling=:directed,
                      mass_matrix=diagonal, sym=symbols)
end


function prosumer_edge() :: nd.ODEEdge
    # TODO: choose between (delta_p, delta_T) or (massflow_fixed, delta_T)
    # Prosumer edges always have dims == 3
    # state 1 => mass flow rate,
    # state 2 => temperature at edge start
    # state 3 => temperature at edge end

    dims = 3
    diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
    symbols = [:m, :T_start, :T_end]

    return nd.ODEEdge(f=Physics.prosumer!, dim=dims, coupling=:directed,
                      mass_matrix=diagonal, sym=symbols)
end


function junction_node() :: nd.DirectedODEVertex
    # dims == 2
    # state 1 => node pressure
    # state 2 => node temperature

    dims = 2
    diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
    symbols = [:p, :T]

    return nd.DirectedODEVertex(f=Physics.junction!, dim=dims,
                                mass_matrix=diagonal, sym=symbols)
end


function fixed_node() :: nd.DirectedODEVertex
    # dims == 2
    # state 1 => node pressure
    # state 2 => node temperature

    dims = 2
    diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
    symbols = [:p, :T]

    return nd.DirectedODEVertex(f=Physics.reference_node!, dim=dims,
                                mass_matrix=diagonal, sym=symbols)
end


end # module DHG
