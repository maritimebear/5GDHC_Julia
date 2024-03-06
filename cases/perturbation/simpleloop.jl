# Script to perturb steady-state solution, then study the time-evolution of the perturbed state

import Graphs as gr
import NetworkDynamics as nd
import DifferentialEquations as de
import SciMLNLSolve

import Plots as plt
import GLMakie, GraphMakie

include("../../src/DHG.jl")
import .DHG

include("./utils_perturb.jl")
import .Utils_Perturb as utils_p


# Parameters

# Perturbation
syms_to_perturb = :T_end_1
perturbation = 1.0
time_interval = (0.0, 24 * 60 * 60.0) # seconds
save_interval = 5 * 60.0 # seconds
save_times = [t for t in time_interval[1] : save_interval : time_interval[end]]

# Graph
g = gr.SimpleDiGraph(4)
edges_g = [(1 => 3), # hot pipe
           (2 => 1), # producer
           (3 => 4), # consumer
           (4 => 2), # cold pipe
          ]
for e in edges_g
    # gr.add_edge!() -> bool to indicate success
    gr.add_edge!(g, e.first, e.second) ? nothing : throw("Failed to add edge $e")
end

fig_graph = GraphMakie.graphplot(g; ilabels=repr.(1:gr.nv(g)), elabels=repr.(1:gr.ne(g)))

# Material properties
# Using material properties of water
const density = 1e3 # kg/m^3
T_ambient::Float64 = 273.15 + 10.0 # kelvin
const fluid_T = DHG.Fluids.Water

# Properties of pipe wall
# Assuming steel as wall material
wall_roughness = 0.045e-3 # m, Cengel table 8-2 (pg.371)

transport_models = DHG.TransportModels(friction_factor=DHG.Transport.friction_Churchill,
                                       Nusselt_number=DHG.Transport.Nusselt_ChiltonCoburn)

## Pipe parameters
pipe_diameter = 1.0
pipe_length = 1.0
pipe_dx = 0.1

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

# Reference node
p_ref = 101325.0 # Pa

# Network components
## Using tuples instead of Vectors: types are not uniform, network structure is constexpr
node_structs = (DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.ReferenceNode(p_ref),
               )

edge_structs = (
                DHG.Pipe(1, 3, pipe_diameter, pipe_length, pipe_dx, wall_roughness),
                DHG.PressureChange(2, 1, producer_hydctrl, producer_thmctrl, producer_hydchar),
                DHG.Massflow(3, 4, consumer_hydctrl, consumer_thmctrl, consumer_hydchar),
                DHG.Pipe(4, 2, pipe_diameter, pipe_length, pipe_dx, wall_roughness),
               )

nodes = (DHG.node(x) for x in node_structs) # Generator
edges = (DHG.edge(x, transport_models, fluid_T) for x in edge_structs) # Generator
params = (density=density, T_ambient=T_ambient)

# Set up problem and solve
nd_fn = nd.network_dynamics(collect(nodes), collect(edges), g)

# Initialise state vector: massflow states must be nonzero, required for node temperature calculation
n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
initial_guess = Vector{Float64}(undef, n_states)

for (k, v) in [(:m => (n) -> rand(Float64, (n, ))), # massflows to random floats
               (:p => p_ref),                       # all pressures to p_ref
               (:T => T_ambient),                   # all temperatures to T_ambient
              ]
    DHG.Miscellaneous.set_idxs(initial_guess, (k => v), nd_fn)
end

# Steady-state solution
sol_steady = utils_p.solve_steadystate(nd_fn, initial_guess, params, SciMLNLSolve.NLSolveJL())

# Dynamic solution, for comparison with steady-state solution
sol_dynamic = utils_p.solve_dynamic(nd_fn, initial_guess, time_interval, params, de.Rodas5())

# Calculate and plot error between steady-state and dynamic solutions
errornorms_initialguess = utils_p.error_norms(sol_dynamic.u, sol_steady.u)
plot_enorms_initialguess = utils_p.plot_errornorms(sol_dynamic.t, errornorms_initialguess,
                                                   "Steady-state solver vs. dynamic solver, from initial_guess")


# Perturb solution
perturbation_vec = zeros(eltype(sol_steady.u), size(sol_steady.u))
idxs_to_perturb = nd.idx_containing(nd_fn, syms_to_perturb)::Vector{<:Integer}
perturbation_vec[idxs_to_perturb] .= perturbation

perturbed_steady = sol_steady.u + perturbation_vec
perturbed_dynamic = sol_dynamic.u[end] + perturbation_vec

# Time-evolution of perturbed state
sol_pertsteady = utils_p.solve_dynamic(nd_fn, perturbed_steady, time_interval, params,
                                       de.Rodas5(), save_times)

sol_pertdynamic = utils_p.solve_dynamic(nd_fn, perturbed_dynamic, time_interval, params,
                                        de.Rodas5(), save_times)


# Calculate and plot variation of error
errornorms_pertsteady = utils_p.error_norms(sol_pertsteady.u, sol_steady.u)
plot_enorms_pertsteady = utils_p.plot_errornorms(sol_pertsteady.t, errornorms_pertsteady,
                                           "Time-evolution of perturbed steady-state solution")



display(plot_enorms_pertsteady)

