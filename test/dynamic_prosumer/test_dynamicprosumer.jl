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
edges_g = [(1 => 2), # producer
           (1 => 3), # hot pipe
           (2 => 4), # cold pipe
           (3 => 4), # consumer
          ]
for e in edges_g
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
pump_ref1 = (0.0, 40221.0, pump_nominalspeed) # (massflow [kg/s], deltaP [Pa], speed [rpm])
pump_ref2 = (55.33 * density, 0.0, pump_nominalspeed)

function hydctrl_pump(nominal_speed)
    function pumpspeed(t)
        let n = nominal_speed
            # t in seconds, returns pump speed in rpm
            if t < (9 * 60 * 60); return 1 * n;
            elseif t < (18 * 60 * 60); return 2 * n;
            else; return 1 * n;
            end
        end
    end
    return pumpspeed
end

function producer_thmpwr(t)
    # t in seconds, returns heat demand in W
    if t < (9 * 60 * 60); return 1e3;
    elseif t < (18 * 60 * 60); return 5e3;
    else; return 1e3;
    end
end

# function consumer_massflow(t) # massflow(t)
#     if t < (9 * 60 * 60); return 1;
#     elseif t < (18 * 60 * 60); return 2;
#     else; return 1;
#     end
# end
# consumer_hydchar = (x, _) -> (x) # Forward control input


producer_hydctrl = hydctrl_pump(pump_nominalspeed)
producer_hydchar = DHG.DynamicalFunctions.PumpModel(pump_ref1..., pump_ref2...,
                                                    density, pump_nominalspeed)

hydraulic_controls = sp.SparseVector(4, [1], [producer_hydctrl]) # (size, index, value)
thermal_controls = sp.SparseVector(4, [1], [producer_thmpwr])
hydraulic_characteristics = sp.SparseVector(4, [1], [producer_hydchar])

# Reference node
p_ref = 101325.0 # Pa
T_fixed = 298.15 # K TODO: Remove T_fixed from ref_pressure node, calculate nodal temperature like junction nodes


prosumer_params = DHG.ParameterStruct.ProsumerParameters(hydraulic_controls, thermal_controls, hydraulic_characteristics)
params = DHG.ParameterStruct.Parameters(density, T_ambient, p_ref, T_fixed, prosumer_params)


nodes::Vector{nd.DirectedODEVertex} = [DHG.WrapperFunctions.junction_node() for _ in 1:3]
push!(nodes, DHG.WrapperFunctions.fixed_node())

edges::Vector{nd.ODEEdge} = [DHG.WrapperFunctions.pipe_edge(diameter, length, dx, transport_coeffs) for _ in 1:3]
pushfirst!(edges, DHG.WrapperFunctions.prosumer_deltaP(1, transport_coeffs)) # (index, coeff_fns)

nd_fn = nd.network_dynamics(nodes, edges, g)

n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
initial_guess = ones(n_states)

prob = de.ODEProblem(nd_fn, initial_guess, (0.0, 24 * 60 * 60), params)
sol = de.solve(prob, de.Rodas5())
