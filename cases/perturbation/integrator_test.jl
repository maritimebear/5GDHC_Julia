import DifferentialEquations as de
import Plots as plt

include("../../src/DHG.jl")
import .DHG.Miscellaneous: solve_dynamic

include("./utils_perturb.jl")
import .Utils_Perturb as utils


# Forced spring-mass system to test time integrator

params = (k = 1.0, # spring constant
          m = 1.0, # mass
         )

solver = de.Tsit5

timespan = (0.0, 200.0) # s
saveinterval = 0.1 # s

# Refine timesteps to simulate decreasing max. dt to maintain CFL num. = 1, as dx
# decreases with mesh refinement
dts = [1.0, 0.5, 0.25] .* 0.01

initial_state = [1.0, 0.0] # [position, velocity]

# Forcing function with ramp
ramp_width = 1.0
(x1, y1) = (100.0, 10.0)
(x2, y2) = (x1 + ramp_width, 20.0)
forcing = utils.piecewise_linear((-1.0, y1), (x1, y1), (x2, y2), (x2 + 1.0, y2))
# forcing = (t) -> (0.0)

function ode(du, u, p, t)
    # x'' = (-k/m)x + (F/m), forced spring-mass system
    x = u[1]
    v = u[2]
    du[1] = v
    du[2] = ((-p.k/p.m) * x) + forcing(t)
    return nothing
end

sols = [solve_dynamic(ode, initial_state, params, solver(), timespan,
                     # saveat=saveinterval,
                     saveat=collect(timespan[1]:dts[1]:timespan[end]),
                     dtmax=dt,
                     # tstops=[x1, x2],
                    ) for dt in dts
       ]

xs = [[u[1] for u in sol.u] for sol in sols]
vs = [[u[2] for u in sol.u] for sol in sols]


# plt.plot(sols[1].t, xs[1])

errors = [abs.(xs[i+1] - xs[i]) for i in 1:(length(xs) - 1)]

convergence = log.(errors[1] ./ errors[2]) ./ log(2)

plt.plot(sols[1].t, convergence)
