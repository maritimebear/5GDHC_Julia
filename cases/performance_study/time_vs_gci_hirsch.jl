# Script to benchmark performance: computation time vs GCI for the Hirsch test case

import DifferentialEquations as de
import BenchmarkTools as bt

include("./setup_hirsch.jl") # provides input parameters (except discretisation) and imports DHG

include("./utils_gci.jl")
import .utils


## Time integration parameters
time_interval = (0.0, 24 * 60.0 * 60.0) # [s]
save_interval = 1 * 60.0 # [s]
save_times = time_interval[1]:save_interval:time_interval[2]
max_CFL = nothing

solver_kwargs = Dict(:saveat => save_times)

## Discretisation
initial_dx = 50.0 # [m]
refinement_ratio = 2
n_refinement_levels = 5
dxs = [initial_dx / (refinement_ratio ^ r) for r in 0:n_refinement_levels]

convection_schemes = Dict("Upwind" => DHG.Discretisation.upwind,
                          "Linear Upwind" => DHG.Discretisation.linear_upwind,
                          "van Leer" => DHG.Discretisation.TVD_vanLeer,
                         )


## Benchmark parameters
benchmarkparameters = Dict(:samples => 10,
                           :seconds => 100.0,
                           :evals => 1,
                          )


# Compare discretisation schemes
for (name, scheme) in convection_schemes

    # println("Convection scheme: $name")
    results[name] = Vector{result_tuple_T}()

    for (idx_dx, dx) in enumerate(dxs)
        discretisation = DHG.Discretisation.FVM(dx=dx, convection=scheme)

        ## Set up problem
        nd_fn, g = DHG.assemble(node_structs, edge_structs, transport_models, discretisation, fluid_T)
            # nd_fn = NetworkDynamics.network_dynamics(collect(node_structs), collect(edge_structs), g)

        ## Initialise state vector: massflow states must be nonzero, required for node temperature calculation
        #   number of massflow and pressure states are constant; == n_nodes
        #   number of temperature states depends on dx
        initial_guess = DHG.Miscellaneous.initialise(nd_fn,
                                                     (_) -> init_massflows, # massflows to init_massflow
                                                     p_ref,                 # pressures to p_ref
                                                     T_ambient              # temperatures to T_ambient
                                                    )

        println("\nConvection scheme: $name, dx = $dx")
        if max_CFL !== nothing
            max_dt = max_CFL * dx / expected_velocity
            println("max dt = $max_dt")
            solver_kwargs[:dtmax] = max_dt
        end

        println("Starting benchmark\n")
        sol, benchmark = utils.benchmark_solve(nd_fn, initial_guess, params, solver(), time_interval,
                                               solver_kwargs, benchmarkparameters
                                              )

        dump(benchmark)

        push!(results[name], (syms=nd_fn.syms, sol=sol.u[end]))
        println("\n --- done ---\n")
    end
end


# Post-processing

node_Ts = [utils.get_states(DHG.PostProcessing.node_T_idxs, results, node_idx) |>
           dict -> Dict(k => [v[1] for v in vs] for (k, vs) in dict) # Unpack 1-element Vector{Float64}
           for (node_idx, _) in enumerate(node_structs)
          ] # Temperature at each node across discn. schemes and dxs


## Grid Convergence -- relative errors, order of convergence, GCI

node_errors = [Dict(scheme => [utils.relative_error(Ts_at_dx[i+1], Ts_at_dx[i]) # relative_error(finer, coarser)
                                for (i, _) in enumerate(Ts_at_dx[1:end-1])
                              ]
                    for (scheme, Ts_at_dx) in Ts_dict
                   )
                for Ts_dict in node_Ts
              ] # Relative error at each node across discn. schemes and (coarser dx, finer dx)

println("\n\nGrid convergence report")
println("-------------------------\n")
println("Grid base size: $initial_dx\n\n")

for (node_idx, _) in enumerate(node_Ts)
    println("Node $node_idx:\n")
    Ts_dict = node_Ts[node_idx]
    errors_dict = node_errors[node_idx]
    for (scheme, _) in Ts_dict
        println("Scheme: $scheme\n")
        for i in 1:3:(3 * div(length(dxs), 3, RoundDown)) # Order of convergence works on sets of 3 dxs
            println("Refinement levels: $(collect(i-1:i+1))")
            println("Grid sizes: $(dxs[i:i+2])")
            p = utils.order_convergence(Ts_dict[scheme][i:i+2], refinement_ratio)
            GCIs = [utils.GCI_fine(errors_dict[scheme][i], refinement_ratio, p),
                    utils.GCI_fine(errors_dict[scheme][i+1], refinement_ratio, p)
                   ]
            println("Order of convergence: $p")
            @show GCIs
            println("\nAsymptotic convergence check: \
                    (GCI ratio / r^p) = $(GCIs[1] / (GCIs[2] * refinement_ratio^p))\n")
        end
        println("\n")
    end
end
