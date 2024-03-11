# Intended to be used only from other scripts in this directory ("discretisation_comparison")

module utils

import Plots as plt


Sequence{T} = Union{Tuple{Vararg{T}}, Vector{T}}
    # type alias for convenience, intended use:
    # Vector{<:T} or Tuple{<:T, <:T, ..., <:T} where T is some abstract type


function get_states(getter_fn, results_dict, elem_idx)
    # -> Dict{String, Vector{T}}: disc. scheme name => vector of element states at each dx
    # getter_fn(syms, elem_idx) converts syms to idxs, elem can be node or edge
    return Dict(scheme_name => [result.sol[getter_fn(result.syms, elem_idx)]
                                for result in results_vec
                               ]
                for (scheme_name, results_vec) in results_dict
               )
end


function relative_error(exact::T, approx::T) where {T <: Real}
    # -> T
    return abs((exact - approx) / exact)
end


function order_convergence(soln_vals::Sequence{Float64}, refinement_ratio::Real)
    # -> Float64
    # Requires:
    #   length(soln_vals) == 3
    #   soln_vals be ordered in decreasing order of grid sizes
    if length(soln_vals) != 3
        error("length(soln_vals) == $(length(soln_vals)), must == 3")
    end
    numerator = soln_vals[2] - soln_vals[1]
    denominator = soln_vals[3] - soln_vals[2]
    return log(numerator / denominator) / log(refinement_ratio)
end


function GCI_fine(error::Float64,
                  refinement_ratio::Real,
                  order_convergence::Float64,
                  factor_safety::Float64 = 1.25
                 )
    # -> Float64
    return factor_safety * abs(error) / ((refinement_ratio^order_convergence) - 1)
end


function plot_pipe_temp(scheme_name, dx_idx, edge_idx)
    dx = dxs[dx_idx]
    x = DHG.PostProcessing.cell_xs(edge_structs[edge_idx].length, dx)
    # x = 0.0:dx:edge_structs[edge_idx].length
    result = results[scheme_name][dx_idx]
    return plt.plot!(x, result.sol[DHG.PostProcessing.edge_T_idxs(result.syms, edge_idx)],
                     label="$scheme_name, dx=$dx", legendposition=:inline,
                     xlabel="x [m]", ylabel="Temperature [K]",
                     title="Temperature profile: edge #$edge_idx",
                    )
end


end # module
