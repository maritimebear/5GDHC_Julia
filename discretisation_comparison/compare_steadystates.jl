# Script to compare steady-state solution for various discretisation schemes and grid sizings

import Graphs as gr
import NetworkDynamics as nd
import DifferentialEquations as de
import SciMLNLSolve
import Random

import Plots as plt
import GLMakie, GraphMakie

include("../src/DHG.jl")
import .DHG

# TODO: actual values for pipe geometry, pump model, thermal controls

# --- Parameters ---- #

Random.seed!(93851203598)
solver = SciMLNLSolve.NLSolveJL

## Fixed/reference values
n_nodes = 4
init_massflows = rand(Float64, (n_nodes, )) # Initial values for massflow states
T_ambient = 273.15 + 10.0 # [K]
p_ref = 101325.0 # [Pa], pressure at reference node

## Material properties, water as fluid
density = 1e3 # [kg/m^3]
fluid_T = DHG.Fluids.Water

## Properties of pipe wall, taking steel as wall material
wall_roughness = 0.045e-3 # [m]
    # Cengel and Cimbala, "Fluid Mechanics: Fundamentals and Applications", 4th ed., table 8-2 (pg.371)

## Pipe parameters, same for all pipes
pipe_diameter = 1.0
pipe_length = 1.0

## Discretisation
# dxs = [1.0 / 2^r for r in 0:3]
dxs = [0.1]
convection_schemes = [DHG.Discretisation.upwind]
discretisations = [DHG.Discretisation.FVM(dx=dx, convection=scheme)
                   for scheme in convection_schemes # outer loop
                   for dx in dxs
                  ]

## Transport properties
transport_models = DHG.TransportModels(friction_factor=DHG.Transport.friction_Churchill,
                                       Nusselt_number=DHG.Transport.Nusselt_ChiltonCoburn)

## Prosumer parameters
pump_nominalspeed = 4100.0 # rpm
pump_ref1 = (0.0, 1.0, pump_nominalspeed) # (massflow [kg/s], deltaP [Pa], speed [rpm])
pump_ref2 = (1.0, 0.0, pump_nominalspeed)

function pumpspeed(t)
    # t in seconds, returns pump speed in rpm
    return 1.0 * pump_nominalspeed
end

function heatinput(t)
    # t in seconds, returns thermal power in W
    return 1e3
end

function massflow_valve(t) # massflow(t)
    # t in seconds, returns massflow in kg/s
    return 0.5
end

consumer_hydctrl = massflow_valve
consumer_hydchar = (ctrl_input, massflow) -> (ctrl_input) # Return control input unchanged
consumer_thmctrl = (t) -> (-0.9 * heatinput(t))


producer_hydctrl = pumpspeed
producer_hydchar = DHG.Miscellaneous.PumpModel(pump_ref1..., pump_ref2...,
                                                    density, pump_nominalspeed)
producer_thmctrl = heatinput

## Network structure
node_structs = (DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.ReferenceNode(p_ref),
               ) # Using tuples instead of Vectors: types are not uniform, network structure is constexpr

edge_structs = (
                DHG.Pipe(1, 3, # src, dst
                         pipe_diameter, pipe_length,  wall_roughness), # hot pipe
                DHG.PressureChange(2, 1,
                                   producer_hydctrl, producer_thmctrl, producer_hydchar), # producer
                DHG.Massflow(3, 4,
                             consumer_hydctrl, consumer_thmctrl, consumer_hydchar), # consumer
                DHG.Pipe(4, 2,
                         pipe_diameter, pipe_length,  wall_roughness), # cold pipe
               )

# --- end of parameters --- #


params = (density=density, T_ambient=T_ambient)

results_dict = Dict{Int,
                    NamedTuple{(:syms, :sol), Tuple{Vector{Symbol}, Vector{Float64}}}
                   }()

# Compare discretisation schemes
for (i, discn) in enumerate(discretisations)
    ## Set up problem
    nd_fn, g = DHG.assemble(node_structs, edge_structs, transport_models, discn, fluid_T)
        # nd_fn = nd.network_dynamics(collect(node_structs), collect(edge_structs), g)

    ## Initialise state vector: massflow states must be nonzero, required for node temperature calculation
    #   number of massflow, pressure states constant; == n_nodes
    #   number of temperature states depends on dx
    initial_guess = DHG.Miscellaneous.initialise(nd_fn,
                                                 (_) -> init_massflows, # massflows to init_massflow
                                                 p_ref,                 # all pressures to p_ref
                                                 T_ambient              # all temperatures to T_ambient
                                                )

    ## Solve for steady-state
    results_dict[i] = (syms=nd_fn.syms,
                       sol=DHG.Miscellaneous.solve_steadystate(nd_fn, initial_guess, params, solver())
                      )

end
