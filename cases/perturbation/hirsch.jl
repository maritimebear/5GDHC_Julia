# Script to perturb steady-state solution, then study the time-evolution of the perturbed state
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


# import Graphs as gr
# import NetworkDynamics as nd
import DifferentialEquations as de
import SciMLNLSolve
import Random

import GLMakie as mk

include("../../src/DHG.jl")
import .DHG


# --- Parameters ---- #

Random.seed!(93851203598)
solver_steady = SciMLNLSolve.NLSolveJL
solver_dynamic = de.Rodas5
time_interval = (0.0, 1 * 60 * 60.0) # [s]
saveinterval = 1.0 # [s]
dx = 10.0 # [m]


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
producer_heatrate = -(1.05* consumer_heatrate) # [W] Assuming heat loss in pipes = 5 to 20% of transmitted energy [Dang]


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
producer_thmctrl = (t) -> (producer_heatrate)


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


params = (density=density, T_ambient=T_ambient)

dynsol_tuple_T = NamedTuple{(:t, :u),
                            Tuple{Vector{Float64}, Vector{Vector{Float64}}}
                           } # DataType: (t = timesteps, u = state values)


Base.@kwdef struct Result # struct to hold results
    syms::Vector{Symbol}
    steady::Vector{Float64}
    dynamic::dynsol_tuple_T
end


results = Dict{String, Result}() # Interpolation scheme name => Result struct


# Compare discretisation schemes
for (name, scheme) in convection_schemes

    println("Convection scheme: $name, dx: $dx")

    discretisation = DHG.Discretisation.FVM(dx=dx, convection=scheme)

    ## Set up problem
    nd_fn, g = DHG.assemble(node_structs, edge_structs, transport_models, discretisation, fluid_T)
        # nd_fn = nd.network_dynamics(collect(node_structs), collect(edge_structs), g)

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

    ## Perturb steady-state solution
    idxs_to_perturb = [findfirst(==(sym), nd_fn.syms) for sym in syms_to_perturb]
    u0_dynamic = Vector{Float64}(sol_steady.u)

    for (idx, rel_ptbn) in zip(idxs_to_perturb, perturbations_relative)
        u0_dynamic[idx] *= (1.0 + rel_ptbn) # Apply perturbation
        # @show (u0_dynamic[idx] - sol_steady.u[idx]) / sol_steady.u[idx] # TODO: Cleanup
    end


    ## Dynamic solution
    print("Starting dynamic solution")
    sol_dynamic = DHG.Miscellaneous.solve_dynamic(nd_fn, u0_dynamic, params, solver_dynamic(),
                                                  time_interval, saveinterval)
    println(" --- done")

    # TODO: Cleanup
    for (idx, rel_ptbn) in zip(idxs_to_perturb, perturbations_relative)
        @show sol_dynamic.t[1]
        @show sol_dynamic.u[1] == u0_dynamic
        @show (sol_dynamic.u[1][idx] - sol_steady.u[idx]) / sol_steady.u[idx]
    end

    results[name] = Result(syms=nd_fn.syms,
                           steady=sol_steady.u,
                           dynamic=(t=sol_dynamic.t, u=sol_dynamic.u)
                          )
end


# Post-processing

## Set up figures and axes
fig_trajectories = mk.Figure()
axes_trajectories = [mk.Axis(fig_trajectories[row, 1],
                            # yticks=mk.LinearTicks(3),
                            # yminorticks=mk.IntervalsBetween(5), yminorticksvisible=true, yminorgridvisible=true,
                            # yticks = mk.MultiplesTicks(5, 1e-3, "a"),
                           ) for (row, _) in enumerate(syms_to_observe)
                    ]

for (i, ax) in enumerate(axes_trajectories)
    # Inset titles for axes
    mk.text!(ax,
             1, 1, text=String(syms_to_observe[i]), font=:bold, #fontsize=11,
             align=(:right, :top), space=:relative, offset=(-8, -14), justification=:right,
            )
end


## Calculate values of interest and plot
for (scheme, result) in results

    idxs_to_observe = [findfirst(==(sym), result.syms) for sym in syms_to_observe]

    states_steady = [result.steady[idx] for idx in idxs_to_observe]::Vector{Float64}
    states_dynamic = [[_u[idx] for _u in result.dynamic.u] for idx in idxs_to_observe]::Vector{Vector{Float64}}

    errors_states = [(u .- u0) ./ u0 for (u, u0) in zip(states_dynamic, states_steady)]::Vector{Vector{Float64}}

    # errors_states = [(u .- result.steady[idx]) ./ result.steady[idx]
    #                  for (u, idx) in zip(states_dynamic, idxs_to_observe)
    #                 ]::Vector{Vector{Float64}}

    for (i, _) in enumerate(idxs_to_observe)
        _ = mk.lines!(axes_trajectories[i],
                      result.dynamic.t, # x-axis
                      errors_states[i], # y-axis
                     ) 
    end

end
