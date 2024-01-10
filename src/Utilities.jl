module Utilities # submodule, included in DHG.jl
# Generic utility functions

export adjacent_find, initialiser

import DifferentialEquations as de
import Graphs


function adjacent_find(binary_predicate, array::AbstractArray)
    # Compare adjacent elements of 'array' using 'binary_predicate'.
    # Returns index of first 'true' from 'binary_predicate'.
    # Returns last index if no match found
    # Following C++ std::adjacent_find(), https://en.cppreference.com/w/cpp/algorithm/adjacent_find
    # No find-like function in Julia that takes a binary predicate?

    if length(array) == 0
        throw(DomainError("invalid argument: empty array"))
    end

    idxs = LinearIndices(array)
    for i in idxs[1:end-1]
        @inbounds if binary_predicate(array[idxs[i]], array[idxs[i+1]])
            return idxs[i]
        end
    end
    return idxs[end] # No match found
end


function initialiser(solver::de.SciMLBase.AbstractSteadyStateAlgorithm)
    # Returns closure: (f, x, parameters) -> Nothing,
    #   x modified in-place: x .= solve(SteadyStateProblem(f, x, parameters), solver)
    #   solve, SteadyStateProblem from DifferentialEquations
    function initialise!(f, x, parameters)
        # Can be used to get initial guess or to re-initialise the solution after a callback
        let solver = solver
            initial_prob = de.SteadyStateProblem(f, x, parameters)
            x .= de.solve(initial_prob, solver).u
            return nothing
        end
    end
    return initialise!
end


function find_edge_idx(graph::Graphs.SimpleGraphs.SimpleDiGraph{T}, src_node_idx::T, dst_node_idx::T) where {T <: Integer}
    # -> T (type used to store node indices in 'graph')
    # Find the index of the first (and assumed only) edge in 'graph' from 'src_node' to 'dst_node'
    # Returns zero(T) if no matching edge is found
    for (i, edge) in enumerate(Graphs.edges(graph))
        if (edge.src == src_node_idx) && (edge.dst == dst_node_idx)
            return T(i)
        end
    end
    return zero(T)
end


end # (sub)module
