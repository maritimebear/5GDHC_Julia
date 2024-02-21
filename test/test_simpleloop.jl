import Graphs as gr
import NetworkDynamics as nd
import DifferentialEquations as de
import SparseArrays as sp

import GLMakie, GraphMakie

include("../src/DHG.jl")
import .DHG


# Parameters

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
diameter = 1.0
length = 1.0
dx = 0.1

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
                DHG.Pipe(1, 3, diameter, length, dx, wall_roughness),
                DHG.PressureChange(2, 1, producer_hydctrl, producer_thmctrl, producer_hydchar),
                DHG.Massflow(3, 4, consumer_hydctrl, consumer_thmctrl, consumer_hydchar),
                DHG.Pipe(4, 2, diameter, length, dx, wall_roughness),
               )

nodes = (DHG.node(x) for x in node_structs) # Generator
edges = (DHG.edge(x, transport_models, fluid_T) for x in edge_structs) # Generator
params = (density=density, T_ambient=T_ambient)

# Set up problem and solve
nd_fn = nd.network_dynamics(collect(nodes), collect(edges), g)

# Initialise state vector: massflow states must be nonzero, required for node temperature calculation
n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
initial_guess = ones(n_states)

for (k, v) in [(:m => (n) -> rand(Float64, (n, ))), # massflows to random floats
               (:p => p_ref),                       # all pressures to p_ref
               (:T => T_ambient),                   # all temperatures to T_ambient
              ]
    DHG.Miscellaneous.set_idxs(initial_guess, (k => v), nd_fn)
end

prob = de.ODEProblem(nd_fn, initial_guess, (0.0, 24 * 60 * 60), params)
sol = de.solve(prob, de.Rodas5())

if sol.retcode !== de.ReturnCode.Success
    throw("Unsuccessful retcode from solver")
end
