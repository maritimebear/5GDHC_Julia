# Script to study the time-evolution of the solution between set-points for various
# discretisation schemes and grid sizings
#
# Find steady-state solution at initial set-point, then dynamic solution to new set-point
#
# References:
#
# Hirsch and Nicolai, "An efficient numerical solution method for detailed modelling of large
#   5th generation district heating and cooling networks", 2022, Section 4.1, Case 1
#
# Licklederer et al, "Thermohydraulic model of Smart Thermal Grids with bidirectional power flow
#   between prosumers", 2021
#
# Rocha et al, "Internal surface roughness of plastic pipes for irrigation", 2017
#
# Dang et al, "Fifth generation district heating and cooling: A comprehensive survey", 2024


import DifferentialEquations as de
import SciMLNLSolve
import Random
import LinearAlgebra as la

import HDF5
import Plots as plt

include("../../src/DHG.jl")
import .DHG

include("./utils_perturb.jl")
import .Utils_Perturb as utils


# --- Parameters ---- #

Random.seed!(93851203598)

export_results = true # Export results to HDF5
export_filename = "perturbation_hirsch.h5"
writemode = "cw" # https://juliaio.github.io/HDF5.jl/stable/#Creating-a-group-or-dataset

solver_steady = SciMLNLSolve.NLSolveJL
solver_dynamic = de.Rodas5

## Spatial discretisation
initial_dx = 20.0 # [m]
refinement_ratio = 2
n_refinement_levels = 5

## Temporal discretisation
time_interval = (0.0, 1 * 60 * 60.0) # [s]
saveinterval = 1.0 # [s]
max_CFL = 1.0 # limits maximum timstep, set to nothing to disable constraint


## Fixed/reference values
n_nodes = 4
init_massflows = rand(Float64, (n_nodes, )) # Initial values for massflow states
T_ambient = 273.15 + 5.0 # [K], Hirsch
p_ref = 0.0 # [Pa], reference pressure


## Pipe geometry, same for all pipes
# Values from Hirsch and Nicolai:
pipe_innerdiameter = 40.8e-3 # [m]
pipe_outerdiameter = 50e-3 # [m]
pipe_length = 100.0 # [m]
wall_conductivity = 0.4 # [W/m-K]


## Prosumers: massflow, thermal power
massflow = 0.3 # [kg/s], Hirsch and Nicolai


## Set-point control
consumer_power_sp1 = -2.7e3 # [W] Assuming temperature change across consumer = -4 K [Hirsch]
consumer_power_sp2 = -5e3 # [W] New set-point
sp_change_time = 15 * 60.0 # [s]
ramp_duration = 1 * 60.0 # [s] Duration of (linear) transition between set-points

# !! The error in the time-integration gets weirder as the ramp gets steeper !! #


## Pump model for producer
# Values from Licklederer et al:
pump_nominalspeed = 4100.0 # rpm
pump_ref1 = (0.0, 40221.0, pump_nominalspeed) # (massflow [kg/s], deltaP [Pa], speed [rpm])
pump_ref2 = (0.922, 0.0, pump_nominalspeed)


## Material properties, taking propylene glycol (Propane-1,3-diol) as fluid [Hirsch]
density = 1064.4 # [kg/m^3], value at 0°C, VDI Heat Atlas D3.1. Table 2 (pg 315)
fluid_T = DHG.Fluids.PropyleneGlycol


## Properties of pipe wall, polyethylene pipe [Hirsch]
wall_roughness = 8.116e-6 # [m], Rocha


## Interpolation schemes
dxs = [initial_dx / (refinement_ratio ^ r) for r in 0:n_refinement_levels]
convection_schemes = Dict("Upwind" => DHG.Discretisation.upwind,
                          "Linear Upwind" => DHG.Discretisation.linear_upwind,
                          "van Leer" => DHG.Discretisation.vanLeer,
                          "van Albada" => DHG.Discretisation.vanAlbada,
                          "MINMOD" => DHG.Discretisation.minmod,
                         )


## Transport properties
transport_models = DHG.TransportModels(friction_factor=DHG.Transport.friction_Churchill,
                                       Nusselt_number=DHG.Transport.Nusselt_ChiltonColburn)


## Prosumer functions
consumer_hydctrl = (t) -> (massflow)
consumer_hydchar = (ctrl_input, massflow) -> (ctrl_input) # Return control input unchanged

producer_hydctrl = (t) -> (pump_nominalspeed)
producer_hydchar = DHG.Miscellaneous.PumpModel(pump_ref1..., pump_ref2..., density, pump_nominalspeed)

consumer_thmctrl_sp1 = (t) -> (consumer_power_sp1)
consumer_thmctrl_sp2 = utils.piecewise_linear((-1.0, consumer_power_sp1),
                                              (sp_change_time, consumer_power_sp1),
                                              (sp_change_time + ramp_duration, consumer_power_sp2),
                                              (ramp_duration + 1.0, consumer_power_sp2)
                                             ) # Piecewise-linear ramp

producer_thmctrl_sp1 = (t) -> (-1.05 * consumer_thmctrl_sp1(t)) # Assuming heat loss in pipes = 5 to 20% of transmitted energy [Dang]
producer_thmctrl_sp2= (t) -> (-1.05 * consumer_thmctrl_sp2(t))


## Network structure
node_structs = (DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.ReferenceNode(p_ref),
               ) # Using tuples instead of Vectors: types are not uniform, network structure is constexpr

edge_structs_sp1 = ( # Used by steady-state solver
                    DHG.Pipe(1, 3, # src, dst
                             pipe_innerdiameter, pipe_outerdiameter, pipe_length,
                             wall_roughness, wall_conductivity), # hot pipe
                    DHG.PressureChange(2, 1,
                                       producer_hydctrl, producer_thmctrl_sp1, producer_hydchar), # producer
                    DHG.Massflow(3, 4,
                                 consumer_hydctrl, consumer_thmctrl_sp1, consumer_hydchar), # consumer
                    DHG.Pipe(4, 2,
                             pipe_innerdiameter, pipe_outerdiameter, pipe_length,
                             wall_roughness, wall_conductivity), # cold pipe
                   )

edge_structs_sp2 = ( # Used by dynamic solver
                    DHG.Pipe(1, 3, # src, dst
                             pipe_innerdiameter, pipe_outerdiameter, pipe_length,
                             wall_roughness, wall_conductivity), # hot pipe
                    DHG.PressureChange(2, 1,
                                       producer_hydctrl, producer_thmctrl_sp2, producer_hydchar), # producer
                    DHG.Massflow(3, 4,
                                 consumer_hydctrl, consumer_thmctrl_sp2, consumer_hydchar), # consumer
                    DHG.Pipe(4, 2,
                             pipe_innerdiameter, pipe_outerdiameter, pipe_length,
                             wall_roughness, wall_conductivity), # cold pipe
                   )
# --- end of parameters --- #


expected_velocity = massflow / (density * 0.25 * pi * pipe_innerdiameter^2)

n_savetimes = Int((time_interval[2] - time_interval[1]) / saveinterval) + 1

params = (density=density, T_ambient=T_ambient)


states_steady = Dict(name =>                                    # convection scheme name
                     Vector{Vector{Float64}}()                  # dx index => nodal temperatures
                     for (name, _) in convection_schemes
                    )::Dict{String, Vector{Vector{Float64}}}
    # Using Vector and not Dict(dx => result) to not use Floats as keys


states_dynamic = Dict(name =>
                      Vector{Array{Float64, 2}}() # dx index => Array(nodal temperature, timestep)
                      for (name, _) in convection_schemes
                     )



# Compare discretisation schemes
for (name, scheme) in convection_schemes
    for (_, dx) in enumerate(dxs)

        discretisation = DHG.Discretisation.FVM(dx=dx, convection=scheme)
        println("\nConvection scheme: $name, dx = $dx")

        ## Steady-state solution -- set-point 1
        ## Set up problem
        nd_fn_sp1, g = DHG.assemble(node_structs, edge_structs_sp1, transport_models, discretisation, fluid_T)
            # nd_fn = nd.network_dynamics(collect(node_structs), collect(edge_structs), g)

        ## Calculate indices of nodal temperatures: constant for a given (scheme, dx)
        nodeT_idxs = [DHG.PostProcessing.node_T_idxs(nd_fn_sp1.syms, node_idx)[1] # unpack 1-element vector
                      for node_idx in 1:length(node_structs)
                     ]

        ## Initialise state vector: massflow states must be nonzero, required for node temperature calculation
        #   number of massflow and pressure states are constant; == n_nodes
        #   number of temperature states depends on dx
        u0_sp1 = DHG.Miscellaneous.initialise(nd_fn_sp1,
                                                     (_) -> init_massflows, # massflows to init_massflow
                                                     p_ref,                 # pressures to p_ref
                                                     T_ambient              # temperatures to T_ambient
                                                    )

        ## Solve for steady-state
        print("Starting steady-state solution")
        sol_steady = DHG.Miscellaneous.solve_steadystate(nd_fn_sp1, u0_sp1, params, solver_steady())
        println(" --- done")

        ## Save node temperatures at steady state
        state_steady = [sol_steady[idx] for idx in nodeT_idxs]

        ## Dynamic solution
        nd_fn_sp2, g = DHG.assemble(node_structs, edge_structs_sp2, transport_models, discretisation, fluid_T)
        u0_sp2 = sol_steady.u

        print("Starting dynamic solution")
        if max_CFL !== nothing
            max_dt = max_CFL * dx / expected_velocity
            print(", max dt = $max_dt")
            sol_dynamic = DHG.Miscellaneous.solve_dynamic(nd_fn_sp2, u0_sp2, params, solver_dynamic(),
                                                          time_interval, saveat=saveinterval,
                                                          dtmax=max_dt,
                                                          dt=1.0,
                                                          tstops=[sp_change_time, sp_change_time + ramp_duration]
                                                         )
        else # no limit on dt
            sol_dynamic = DHG.Miscellaneous.solve_dynamic(nd_fn_sp2, u0_sp2, params, solver_dynamic(),
                                                          time_interval, saveat=saveinterval,
                                                          tstops=[sp_change_time, sp_change_time + ramp_duration]
                                                         )
        end
        println(" --- done")

        ## Save node temperatures at final dynamic state
        state_dynamic = [sol_dynamic[time_idx][nodeT_idx]
                         for nodeT_idx in nodeT_idxs, time_idx in 1:length(sol_dynamic)
                        ]

        push!(states_steady[name], state_steady)
        push!(states_dynamic[name], state_dynamic)

    end # loop over dxs
end # loop over convection schemes


# Post-processing

statediffs_steady = Dict(scheme => [dx_vec[i+1] .- dx_vec[i] for i in 1:(length(dx_vec)-1)]
                         for (scheme, dx_vec) in states_steady
                        )

statediffs_dynamic = Dict(scheme => [dx_mat[i+1] .- dx_mat[i] for i in 1:(length(dx_mat)-1)]
                         for (scheme, dx_mat) in states_dynamic
                        )

errors_dynamic = Dict(scheme => [abs.(mat) for mat in dx_mats] for (scheme, dx_mats) in statediffs_dynamic)

# convergence = Dict{String, Matrix{Float64}}()
# for (scheme, errors) in statediffs_dynamic # statediffs_dynamic for no abs()
#     for i in 1:2:(2 * div(length(errors), 2, RoundDown))
#         numerator = errors[i]
#         denominator = errors[i+1]
#         p = log.(numerator ./ denominator) / log(refinement_ratio)
#         convergence[scheme] = p
#     end
# end


maxerrors_dynamic = Dict(scheme_name =>
                         [[la.norm(statediffs_dynamic[scheme_name][refinement_idx][node_idx, :],
                                   Inf
                                  ) for node_idx in 1:length(node_structs)
                          ] for refinement_idx in 1:n_refinement_levels
                         ] for (scheme_name, _) in convection_schemes
                        )

# ## Order of convergence calculation from grid convergence analysis
# order_convergence = Dict{String, Vector{Float64}}()
# for (scheme, errors) in maxerrors_dynamic
#     for i in 1:3:(3 * div(length(errors), 3, RoundDown)) # Order of convergence works on sets of 3 dxs
#         order_convergence[scheme] = log.((errors[i+1] .- errors[i]) ./ (errors[i+2] .- errors[i+1])) / log(refinement_ratio)
#     end
# end


times = time_interval[1]:saveinterval:time_interval[2]

if export_results
    print("Exporting results ... ")
    HDF5.h5open(export_filename, writemode) do fid
    fid["times"] = collect(times)
    fid["dxs"] = dxs
        for (scheme_name, _) in convection_schemes
            g = HDF5.create_group(fid, scheme_name)
            for (dx_idx, dx) in enumerate(dxs)
                dx_g = HDF5.create_group(g, "dx_$dx_idx")
                dx_g["dx"] = dx
                dx_g["nodes_temperatures"] = states_dynamic[scheme_name][dx_idx]
            end
        end
    end
    println("done")
end

plots_node_Ts = [plt.plot(times, [states_dynamic["Upwind"][dx_idx][node_idx, :] for dx_idx in 1:length(dxs)],
               # label="dx: $(dxs[dx_idx])",
               legendposition=:inline,
               legendfontsize=3,
              )
      for node_idx in 1:length(node_structs)
     ]
fig_node_Ts = plt.plot(plots_node_Ts...,
              layout=4,
              title=reshape(["Node $i" for i in 1:length(plots_node_Ts)], (1, length(plots_node_Ts))),
              titlefontsize=8,
             )

# plots_maxerrors = [plt.plot([errors[refinement_idx][node_idx] for refinement_idx in 1:n_refinement_levels],
#                             label=scheme,
#                             legendposition=:inline,
#                             legendfontsize=4,
#                            ) for node_idx in 1:length(node_structs)
#                    for (scheme, errors) in maxerrors_dynamic
#                   ]

# fig_maxerrors = plt.plot(plots_maxerrors...,
#                          layout=4,
#                          title=reshape(["Node $i" for i in 1:length(plots_maxerrors)], (1, length(plots_maxerrors))),
#                          titlefontsize=8,
#                         )

# plots_convergence = [plt.plot(times, conv[node_idx, :],
#                             label=scheme,
#                             legendposition=:inline,
#                             legendfontsize=4,
#                            ) for node_idx in 1:length(node_structs)
#                    for (scheme, conv) in convergence
#                   ]

# fig_convergence = plt.plot(plots_convergence...,
#                          layout=4,
#                          title=reshape(["Node $i" for i in 1:length(plots_convergence)], (1, length(plots_convergence))),
#                          titlefontsize=8,
#                         )
# display(fig_convergence)
