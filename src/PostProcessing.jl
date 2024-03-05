module PostProcessing

# Convenience functions for plotting etc.

import Graphs

import GLMakie as mk
import GraphMakie

export syms_to_idxs, edge_T_idxs, edge_idxs, node_T_idxs, node_idxs
export cell_xs
export plot_graph


# Functions to find indices of states in ODESolution from vector of symbols
# Does the same job as NetworkDynamics.idx_containing(), without needing a network_dynamics object
function syms_to_idxs(syms::Vector{Symbol}, pattern::Union{AbstractString, Regex})
    # -> Vector{Int}
    return [idx for (idx, str) in enumerate(String.(syms)) if occursin(pattern, str)]
end

function syms_to_idxs(syms::Vector{Symbol}, patterns::Vector{P}) where {P <: Union{AbstractString, Regex}}
    # -> Vector{Vector{Int}}
    return [syms_to_idxs(syms, pattern) for pattern in patterns]
end


# Functions to get indices of states in edges and nodes

function edge_T_idxs(syms::Vector{Symbol}, edge_idx::Integer)
    # -> Vector{Int}
    # Return idxs of temperature states in the specifed edge
    return syms_to_idxs(syms, Regex("T_\\w+_$(edge_idx)"))
end


function edge_idxs(syms::Vector{Symbol}, edge_idx::Integer)
    # -> Vector{Int}
    # Returns idxs of [massflow, T_1 : T_end] states in specified edge
    return cat(syms_to_idxs(syms, "m_$(edge_idx)"), edge_T_idxs(syms, edge_idx), dims=1)
end


function node_T_idxs(syms::Vector{Symbol}, node_idx::Integer)
    # -> Vector{Int}
    # Return index of node-temperature state
    # Return type is vector for consistency with edge-state access
    return syms_to_idxs(syms, Regex("T_$(node_idx)\$"))
end


function node_idxs(syms::Vector{Symbol}, node_idx::Integer)
    # -> Vector{Int}
    # Return index of [pressure, temperature] states in specified node
    return cat(syms_to_idxs(syms, "p_$(node_idx)"), node_T_idxs(syms, node_idx), dims=1)
end


function cell_xs(length::Float64, dx::Float64, x_start::Float64=0.0)
    # -> Vector{Float64}
    # Calculates x-locations of FVM cell-centres in pipe edges
    return collect(range(x_start + dx/2, x_start + length, step=dx))
end


# Plotting functions

function plot_graph(g::Graphs.SimpleDiGraph)
    # -> handle to plot
    # Displays graph g with nodes and edges numbers (following implicit Graphs.jl ordering)
    return GraphMakie.graphplot(g; ilabels=repr.(1:Graphs.nv(g)), elabels=repr.(1:Graphs.ne(g)))
end


# function plot_pipe_temperature(axis::mk.Axis,
#                                 solution::Vector{Float64},
#                                 state_syms::Vector{Symbol},
#                                 edge_idx::Integer,
#                                 dx::Float64
#                                 ; plot_kwargs...
#                             )
#     # -> Makie.Lines
#     # Plots temperature distribution along a pipe edge
#     # plot_kwargs forwarded to Makie.lines!()

#     # return mk.lines!(axis, x, y; plot_kwargs...)
# end


end # module
