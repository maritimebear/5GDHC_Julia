import Graphs as gr
import NetworkDynamics as nd
import DifferentialEquations as de

import GLMakie, GraphMakie

include("../../src/DHG.jl")
import .DHG


# Parameters
# graph_definition = "./single_prosumer.gml"
p_ref = 101325.0 # Pa
T_fixed = 298.15 # K TODO: Remove T_fixed from ref_pressure node, calculate nodal temperature like junction nodes

## Pipe parameters
diameter = 1.0
length = 1.0
dx = 0.1

## Prosumer parameters
function pressure_control(t)
    # t in seconds, returns deltaP in Pa
    if t < (9 * 60 * 60); return -5;
    elseif t < (18 * 60 * 60); return -10;
    else; return -5;
    end
end

function heatrate_control(t)
    # t in seconds, returns heat demand in W
    if t < (9 * 60 * 60); return -100;
    elseif t < (18 * 60 * 60); return -500;
    else; return -100;
    end
end

hydraulic_characteristic = (deltaP) -> (-deltaP)

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

# Graph
g = gr.SimpleDiGraph(4)
edges_g = [(1 => 2), # producer
           (1 => 3), # hot pipe
           (2 => 4), # cold pipe
           (3 => 4), # consumer
          ]
for e in edges_g
    gr.add_edge!(g, e.first, e.second) ? nothing : throw("Failed to add edge $e")
end

fig_graph = GraphMakie.graphplot(g; ilabels=repr.(1:gr.nv(g)), elabels=repr.(1:gr.ne(g)))

# params = DHG.Parameters(density=density, T_ambient=T_ambient, p_ref=p_ref, T_fixed=T_fixed)

# nodes::Vector{nd.DirectedODEVertex} = [DHG.WrapperFunctions.junction_node() for _ in 1:3]
# push!(nodes, DHG.WrapperFunctions.fixed_node())

# edges::Vector{nd.ODEEdge} = [DHG.WrapperFunctions.pipe_edge(diameter, length, dx, transport_coeffs) for _ in 1:3]
# pushfirst!(edges, DHG.WrapperFunctions.prosumer(pressure_control, heatrate_control,
#                                                 hydraulic_characteristic, transport_coeffs))

# nd_fn = nd.network_dynamics(nodes, edges, g)

# n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
# initial_guess = ones(n_states)

# prob = de.ODEProblem(nd_fn, initial_guess, (0.0, 24 * 60 * 60), params)
# sol = de.solve(prob, de.Rodas5())
