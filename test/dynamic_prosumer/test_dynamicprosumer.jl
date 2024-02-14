import Graphs as gr
import NetworkDynamics as nd
import DifferentialEquations as de
import SparseArrays as sp

import GLMakie, GraphMakie

include("../../src/DHG.jl")
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
heat_capacity::Float64 = 4184.0 # J/kg-K, used only in heat transfer coefficient calculation
dyn_visc::Float64 = 8.9e-4 # Pa-s

# Properties of pipe wall
# Assuming steel as wall material
wall_conductivity = 30.0 # W/m-K, NOT A RELIABLE VALUE !!!
wall_thickness = 10e-3 # m

friction = (Re) -> (-1.0)
h_wall::Float64 = 2 * wall_conductivity / wall_thickness # Wikipedia: Heat transfer coefficient -- Heat transfer coefficient of pipe wall
heat_transfer::Float64 = -h_wall / (density * heat_capacity) # TODO: Why is this -ve?

transport_coeffs = DHG.TransportProperties(dynamic_viscosity=dyn_visc, wall_friction=friction,
                                             heat_capacity=heat_capacity, heat_transfer=heat_transfer)

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
    return 1.0
end

consumer_hydctrl = massflow_valve
consumer_hydchar = (ctrl_input, massflow) -> (ctrl_input) # Return control input unchanged
consumer_thmctrl = (t) -> (-0.9 * heatinput(t))


producer_hydctrl = pumpspeed
producer_hydchar = DHG.ControlFunctions.PumpModel(pump_ref1..., pump_ref2...,
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

# edge_structs = (DHG.Pipe(1, 2, diameter, length, dx),
edge_structs = (
                DHG.Pipe(1, 3, diameter, length, dx),
                DHG.Massflow(2, 1, consumer_hydctrl, consumer_thmctrl, consumer_hydchar),
                DHG.Pipe(3, 4, diameter, length, dx),
                DHG.Pipe(4, 2, diameter, length, dx),
               )

nodes = (DHG.node(x) for x in node_structs) # Generator
edges = (DHG.edge(x, transport_coeffs) for x in edge_structs) # Generator
params = (density=density, T_ambient=T_ambient)

# Set up problem and solve
nd_fn = nd.network_dynamics(collect(nodes), collect(edges), g)

n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
initial_guess = ones(n_states)

prob = de.ODEProblem(nd_fn, initial_guess, (0.0, 24 * 60 * 60), params)
sol = de.solve(prob, de.Rodas5())
