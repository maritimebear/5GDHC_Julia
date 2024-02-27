# Utility functions to study time-evolution of perturbed states

module Utils_Perturb

import DifferentialEquations as de
import LinearAlgebra as la
import Plots as plt

export error_norms, plot_errornorms

function error_norms(dynamic_solutions::Vector{Vector{Float64}}, reference_solution::Vector{Float64})
    # -> Vector{Float64}
    # Calculates 2-norms of error for each Vector{Float64} in dynamic_solutions
    return [la.norm(x - reference_solution) for x in dynamic_solutions]::Vector{Float64}
end


function plot_errornorms(time_vec::Vector{Float64}, errornorms::Vector{Float64}, title::AbstractString)
    # -> handle to plot
    # Convenience function adding labels etc. to plots
    p = plt.plot(time_vec, errornorms,
                 xlabel="time (s)",
                 ylabel="2-norm of error",
                 title=title,
                 yaxis=:log,
                )
    return p
end

end # module
