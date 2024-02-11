module WrapperFunctions

import NetworkDynamics as nd
import LinearAlgebra as la
import ..DynamicalFunctions
import ..Transport.TransportProperties
import ..NetworkComponents as nc

export edge, node, pipe_edge, prosumer_edge, junction_node, reference_node

function edge(edge_struct::nc.Pipe, transport_coeffs::TransportProperties)
    return pipe_edge(edge_struct, transport_coeffs)
end

function edge(edge_struct<:nc.Prosumer, transport_coeffs::TransportProperties)
    return prosumer_edge(edge_struct, transport_coeffs)
end

function node(::nc.JunctionNode)
    return junction_node()
end

function node(node_struct::nc.ReferenceNode)
    return reference_node(node_struct)
end


function pipe_edge(pipe_struct::nc.Pipe, transport_coeffs::TransportProperties)
    # -> nd.ODEEdge
    # Edge with friction-induced pressure drop, temperature loss to surroundings
    # state 1 => mass flow rate
    # states 2:end => temperature in finite-volume cells
    # dims == num(states) = 1 + num(cells)

    n_cells::Int = div(pipe_struct.length, pipe_struct.dx, RoundNearest) # number of finite-volume cells, integer division
    if (n_cells <= 0 || n_cells == Inf)
        # Change n_cells to UInt to catch these errors, except n_cells == 0
        throw(ArgumentError("calculated n_cells: $n_cells"))
    end

    # Calculate arguments to NetworkDynamics.ODEEdge
    f = DynamicalFunctions.pipe(pipe_struct, transport_coeffs)
    dims = n_cells + oneunit(n_cells) # type-stable, type-agnostic increment
    diagonal = la.Diagonal([1 for _ in 1:dims]) # Diagonal of mass matrix
    diagonal[1] = 0 # state 1 corresponds to an algebraic constraint
    symbols = [Symbol("T$i") for i in 0 : n_cells]
    symbols[1] = :m

    return nd.ODEEdge(f=f, dim=dims, coupling=:directed, mass_matrix=diagonal, sym=symbols)
end


function prosumer_edge(prosumer_struct <: nc.Prosumer, transport_coeffs::TransportProperties)
    # -> nd.ODEEdge
    # Prosumer edges always have dims == 2
    # state 1 => mass flow rate,
    # state 2 => temperature at prosumer outlet

    if typeof(prosumer_struct) == nc.PressureChange 
        dyn_fn = DynamicalFunctions.prosumer_deltaP
    else # i.e. type == nc.Massflow
        dyn_fn = DynamicalFunctions.prosumer_massflow
    end

    f = dyn_fn(prosumer_struct, transport_coeffs)
    dims = 2
    diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
    symbols = [:m, :T_outlet]

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


function reference_node(node_struct::nc.ReferenceNode) # -> nd.DirectedODEVertex
    # dims == 2
    # state 1 => node pressure
    # state 2 => node temperature

    f = DynamicalFunctions.reference_node(node_struct)
    dims = 2
    diagonal = la.Diagonal([0 for _ in 1:dims]) # Diagonal of mass matrix
    symbols = [:p, :T]

    return nd.DirectedODEVertex(f=f, dim=dims, mass_matrix=diagonal, sym=symbols)
end


end # module
