module WrapperFunctions # submodule, included in DHG.jl

import NetworkDynamics as nd
import LinearAlgebra as la
import ..DynamicalFunctions
import ..Transport.TransportProperties
import ..NetworkComponents

export pipe_edge, prosumer_massflow, prosumer_deltaP, junction_node, fixed_node

## Wrapper functions around NetworkDynamics structs to calculate arguments to their constructors

# TODO: diagonals to Bool instead of Integer
# TODO: StaticVector for symbols, but NetworkDynamics only takes Vector{Symbol}


# # Multiple dispatch on NetworkComponents structs

# function node_fn(_::NetworkComponents.Junction)
#     return junction_node()
# end


# function node_fn(_::NetworkComponents.FixedNode)
#     return fixed_node()
# end


# function edge_fn(_::Integer, edge::NetworkComponents.Pipe, transport_coeffs::TransportProperties)
#     return pipe_edge(edge.diameter, edge.length, edge.dx, transport_coeffs)
# end


# function edge_fn(idx::Integer, _::NetworkComponents.MassflowProsumer, _::TransportProperties)
#     return prosumer_massflow(idx)
# end


# function edge_fn(idx::Integer, _::NetworkComponents.PressureChangeProsumer, _::TransportProperties)
#     return prosumer_deltaP(idx)
# end


# Methods for each NetworkComponents struct

function pipe_edge(diameter::Real, length::Real, dx::Real, coeff_fns::TransportProperties)
    # -> nd.ODEEdge
    # Edge with friction-induced pressure drop, temperature loss to surroundings
    # state 1 => mass flow rate
    # states 2:end => temperature in finite-volume cells
    # dims == num(states) = 1 + num(cells)

    n_cells::Int = div(length, dx, RoundNearest) # number of finite-volume cells, integer division
    if (n_cells <= 0 || n_cells == Inf)
        # Change n_cells to UInt to catch these errors, except n_cells == 0
        throw(ArgumentError("calculated n_cells: $n_cells")) # TODO: add idx to error msg
    end

    # Calculate arguments to NetworkDynamics.ODEEdge
    f = DynamicalFunctions.pipe(diameter, dx, coeff_fns)
    dims = n_cells + oneunit(n_cells) # type-stable, type-agnostic increment
    diagonal = la.Diagonal([1 for _ in 1:dims]) # Diagonal of mass matrix
    diagonal[1] = 0 # state 1 corresponds to an algebraic constraint
    symbols = [Symbol("T$i") for i in 0 : n_cells]
    symbols[1] = :m

    return nd.ODEEdge(f=f, dim=dims, coupling=:directed, mass_matrix=diagonal, sym=symbols)
end


# function prosumer_massflow(index::Integer) # -> nd.ODEEdge
#     # TODO: choose between (delta_p, delta_T) or (massflow_fixed, delta_T)
#     # Prosumer edges always have dims == 3
#     # state 1 => mass flow rate,
#     # state 2 => temperature at edge start
#     # state 3 => temperature at edge end

#     f = DynamicalFunctions.prosumer_massflow(index)
#     dims = 3
#     diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
#     symbols = [:m, :T_start, :T_end]

#     return nd.ODEEdge(f=f, dim=dims, coupling=:directed, mass_matrix=diagonal, sym=symbols)
# end


# function prosumer_deltaP(index::Integer) # -> nd.ODEEdge
#     # TODO: choose between (delta_p, delta_T) or (massflow_fixed, delta_T)
#     # Prosumer edges always have dims == 3
#     # state 1 => mass flow rate,
#     # state 2 => temperature at edge start
#     # state 3 => temperature at edge end

#     f = DynamicalFunctions.prosumer_deltaP(index)
#     dims = 3
#     diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
#     symbols = [:deltaP, :T_start, :T_end]

#     return nd.ODEEdge(f=f, dim=dims, coupling=:directed, mass_matrix=diagonal, sym=symbols)
# end


function prosumer(pressure_control::Function,
                    heatrate_control::Function,
                    hydraulic_characteristic::Function,
                    coeff_fns::TransportProperties
    )
    # -> nd.ODEEdge
    # Prosumer edge, pressure change and heat transfer are functions of time
    # state 1 => mass flow rate
    # state 2 => inlet temperature
    # state 3 => outlet temperature

    # Calculate arguments to NetworkDynamics.ODEEdge
    f = DynamicalFunctions.prosumer(pressure_control, heatrate_control, hydraulic_characteristic, coeff_fns)
    dims = 3
    diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
    symbols = [:m, :T_in, :T_out]

    return nd.ODEEdge(f=f, dim=dims, coupling=:directed, mass_matrix=diagonal, sym=symbols)
end


function junction_node() # -> nd.DirectedODEVertex
    # dims == 2
    # state 1 => node pressure
    # state 2 => node temperature

    f = DynamicalFunctions.junction!
    dims = 2
    diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
    symbols = [:p, :T]

    return nd.DirectedODEVertex(f=f, dim=dims, mass_matrix=diagonal, sym=symbols)
end


function fixed_node() # -> nd.DirectedODEVertex
    # dims == 2
    # state 1 => node pressure
    # state 2 => node temperature

    f = DynamicalFunctions.fixed_node!
    dims = 2
    diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
    symbols = [:p, :T]

    return nd.DirectedODEVertex(f=f, dim=dims, mass_matrix=diagonal, sym=symbols)
end


end # module
