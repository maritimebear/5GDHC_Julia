# Script to perturb steady-state solution, then study the time-evolution of the perturbed state
# for various discretisation schemes and grid sizings
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

import Plots as plt

include("../../src/DHG.jl")
import .DHG


# --- Parameters ---- #

Random.seed!(93851203598)
solver_steady = SciMLNLSolve.NLSolveJL
solver_dynamic = de.Rodas5

## Spatial discretisation
initial_dx = 10.0 # [m]
refinement_ratio = 2
n_refinement_levels = 3

## Temporal discretisation
time_interval = (0.0, 1 * 60 * 60.0) # [s]
saveinterval = 1.0 # [s]
max_CFL = 1.0 # limits maximum timstep, set to nothing to disable constraint


## Perturbation
syms_to_perturb = [:T_end_1]
perturbations_relative = [0.01] # relative magnitudes
syms_to_observe = [:T_3, :T_1]


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
consumer_heatrate = -2.7e3 # [W] Assuming temperature change across consumer = -4 K [Hirsch]


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
                         )


## Transport properties
transport_models = DHG.TransportModels(friction_factor=DHG.Transport.friction_Churchill,
                                       Nusselt_number=DHG.Transport.Nusselt_ChiltonColburn)


## Prosumer functions
consumer_hydctrl = (t) -> (massflow)
consumer_hydchar = (ctrl_input, massflow) -> (ctrl_input) # Return control input unchanged
consumer_thmctrl = (t) -> (consumer_heatrate)

producer_hydctrl = (t) -> (pump_nominalspeed)
producer_hydchar = DHG.Miscellaneous.PumpModel(pump_ref1..., pump_ref2..., density, pump_nominalspeed)
producer_thmctrl = (t) -> (-1.05 * consumer_thmctrl(t)) # Assuming heat loss in pipes = 5 to 20% of transmitted energy [Dang]


## Network structure
node_structs = (DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.ReferenceNode(p_ref),
               ) # Using tuples instead of Vectors: types are not uniform, network structure is constexpr

edge_structs = (
                DHG.Pipe(1, 3, # src, dst
                         pipe_innerdiameter, pipe_outerdiameter, pipe_length,
                         wall_roughness, wall_conductivity), # hot pipe
                DHG.PressureChange(2, 1,
                                   producer_hydctrl, producer_thmctrl, producer_hydchar), # producer
                DHG.Massflow(3, 4,
                             consumer_hydctrl, consumer_thmctrl, consumer_hydchar), # consumer
                DHG.Pipe(4, 2,
                         pipe_innerdiameter, pipe_outerdiameter, pipe_length,
                         wall_roughness, wall_conductivity), # cold pipe
               )

# --- end of parameters --- #


if size(syms_to_perturb) != size(perturbations_relative)
    error("Sizes of syms_to_perturb, perturbations_relative must be the same")
end

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

        ## Set up problem
        nd_fn, g = DHG.assemble(node_structs, edge_structs, transport_models, discretisation, fluid_T)
            # nd_fn = nd.network_dynamics(collect(node_structs), collect(edge_structs), g)

        ## Calculate indices of nodal temperatures: constant for a given (scheme, dx)
        nodeT_idxs = [DHG.PostProcessing.node_T_idxs(nd_fn.syms, node_idx)[1] # unpack 1-element vector
                      for node_idx in 1:length(node_structs)
                     ]

        ## Initialise state vector: massflow states must be nonzero, required for node temperature calculation
        #   number of massflow and pressure states are constant; == n_nodes
        #   number of temperature states depends on dx
        initial_guess = DHG.Miscellaneous.initialise(nd_fn,
                                                     (_) -> init_massflows, # massflows to init_massflow
                                                     p_ref,                 # pressures to p_ref
                                                     T_ambient              # temperatures to T_ambient
                                                    )

        ## Solve for steady-state
        print("Starting steady-state solution")
        sol_steady = DHG.Miscellaneous.solve_steadystate(nd_fn, initial_guess, params, solver_steady())
        println(" --- done")

        ## Save node temperatures at steady state
        state_steady = [sol_steady[idx] for idx in nodeT_idxs]

        ## Perturb steady-state solution
        idxs_to_perturb = [findfirst(==(sym), nd_fn.syms) for sym in syms_to_perturb]
        u0_dynamic = Vector{Float64}(sol_steady.u)

        for (idx, rel_ptbn) in zip(idxs_to_perturb, perturbations_relative)
            u0_dynamic[idx] *= (1.0 + rel_ptbn) # Apply perturbation
        end

        ## Dynamic solution
        print("Starting dynamic solution")
        if max_CFL !== nothing
            max_dt = max_CFL * dx / expected_velocity
            print(", max dt = $max_dt")
            sol_dynamic = DHG.Miscellaneous.solve_dynamic(nd_fn, u0_dynamic, params, solver_dynamic(),
                                                          time_interval, saveat=saveinterval,
                                                          dtmax=max_dt)
        else # no limit on dt
            sol_dynamic = DHG.Miscellaneous.solve_dynamic(nd_fn, u0_dynamic, params, solver_dynamic(),
                                                          time_interval, saveat=saveinterval)
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

# node_Es = [abs.(statediffs_dynamic["Upwind"][i][3, :]) for i in 1:n_refinement_levels]
# times = time_interval[1]:saveinterval:time_interval[2]
# p = plt.plot(times, node_Es[1])
# for v in node_Es[2:end]
#     p = plt.plot!(times, v)
# end
# display(p)

node_Ts = [states_dynamic["Upwind"][i][3, :] for i in 1:length(dxs)]
times = time_interval[1]:saveinterval:time_interval[2]
p = plt.plot(times, node_Ts[1])
for v in node_Ts[2:end]
    plt.plot!(times, v)
end
display(p)
