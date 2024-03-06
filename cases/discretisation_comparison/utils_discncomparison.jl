# Intended to be used only from other scripts in this directory ("discretisation_comparison")

module utils

import Plots as plt


function get_states(getter_fn, results_dict, elem_idx)
    # -> Dict{String, Vector{Float64}}: disc. scheme name => vector of element states at each dx
    # getter_fn(syms, elem_idx) converts syms to idxs, elem can be node or edge
    return Dict(scheme_name => [result.sol[getter_fn(result.syms, elem_idx)]
                                for result in results_vec
                               ]
                for (scheme_name, results_vec) in results_dict
               )
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
