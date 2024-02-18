# Utility functions to study time-evolution of perturbed states

module Utils_Perturb

import DifferentialEquations as de
import LinearAlgebra as la
import Plots as plt

export error_norms, solve_steadystate, solve_dynamic, plot_errornorms

function error_norms(dynamic_solutions::Vector{Vector{Float64}}, reference_solution::Vector{Float64})
    # -> Vector{Float64}
    # Calculates 2-norms of error for each Vector{Float64} in dynamic_solutions
    return [la.norm(x - reference_solution) for x in dynamic_solutions]::Vector{Float64}
end


function check_retcode(solution, solvername)
    # Checks return code from de.solve() result, raise warning if not successful
    if solution.retcode !== de.ReturnCode.Success
        @warn "Unsuccessful retcode from $(solvername): $(solution.retcode)"
    end
end


function solve_steadystate(f, x0, p, solver)
    # -> de.solve(de.SteadyStateProblem, solver) type
    # Wrapper to make main script more readable
    prob = de.SteadyStateProblem(f, x0, p)
    sol = de.solve(prob, solver)
    check_retcode(sol, "steady-state solver")
    return sol
end


function solve_dynamic(f, x0, time_interval, p, solver, save_times=[])
    # -> de.solve(de.ODEProblem, solver, saveat=save_times) type
    # Wrapper to make main script more readable
    prob = de.ODEProblem(f, x0, time_interval, p)
    sol = de.solve(prob, solver, saveat=save_times)
    check_retcode(sol, "dynamic solver")
    return sol
end


function plot_errornorms(time_vec::Vector{Float64}, errornorms::Vector{Float64}, title::AbstractString)
    # -> handle to plot
    # Convenience function adding labels etc. to plots
    p = plt.plot(time_vec, errornorms,
                 xlabel="time (s)",
                 ylabel="2-norm of error",
                 title=title,
                )
    return p
end

end # module
