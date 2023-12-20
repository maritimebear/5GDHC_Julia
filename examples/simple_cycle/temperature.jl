import Graphs as gr
import NetworkDynamics as nd
import LinearAlgebra as la
import DifferentialEquations as de
import GLMakie, GraphMakie

# Only temperature transport, no hydrodynamics

include("../../src/InterpolationSchemes.jl") # provides module FVM

# Graph
g = gr.cycle_digraph(4)
fig_graph = GraphMakie.graphplot(g; ilabels=repr.(1:gr.nv(g)), elabels=repr.(1:gr.ne(g)))

# Parameters
# Using material properties of water
density = 1e3 # kg/m^3
heat_capacity = 4184.0 # J/kg-K, used only in heat transfer coefficient calculation
T_ambient = 273.15 + 10.0 # kelvin

# Properties of pipe wall
# Assuming steel as wall material
wall_conductivity = 30.0 # W/m-K, NOT A RELIABLE VALUE !!!
wall_thickness = 10e-3 # m

# Parameters for each edge
diameters = [1.0 for _ in 1:gr.ne(g)]
massflows = [1.0 for _ in 1:gr.ne(g)]
pipelengths = [1.0 for _ in 1:gr.ne(g)]
n_cells = [10 for _ in 1:gr.ne(g)] # number of cells in each edge
dx = pipelengths ./ n_cells # cell widths of each edge, not cell-centre locations

h_wall = () -> (2 * wall_conductivity / wall_thickness) # Wikipedia: Heat transfer coefficient -- Heat transfer coefficient of pipe wall
htrans_coeff = () -> (-h_wall() / (density * heat_capacity))

# Named tuple to be passed to dynamical functions
params = (density = density,
          diameter = diameters[1],
          massflow = massflows[1],
          dx = dx[1],
          delta_T = 0.0,
          T_fixed = 273.15 + 25.0,
          T_ambient = T_ambient,
          htrans_coeff = htrans_coeff,
         )


# TODO: Function to create AoS edge_params and node_params, mapping edge/node idx to type and parameters
#       eg. prosumer_edges don't need dx, diameter, heat_capacity but do need delta_T; pipe_edges are the other way around
# edge_params = (density = 1e3,
#           diameter = 1.0,
#           heat_capacity = 4184.0, # J/kg-K
#           delta_T = 5.0,
#           htrans_coeff = (_) -> (-1.0),
#           dx = pipelengths ./ n_cells # cell widths of each edge, not cell-centre locations
#          )

# node_params = (density = 1e3,
#           T_fixed = 273.15 + 25.0, # kelvin
#           massflow = [1.0 for _ in 1:gr.ne(g)], # Same massflow for all edges
#          )



function pipe_edge!(de, e, v_s, v_d, p, _)
    # Calculate local variables
    area= 0.25 * pi * (p.diameter ^ 2)
    velocity= p.massflow / (p.density * area)

    convection_term = -(1 / p.dx) .* FVM.upwind(e, v_s[1], v_d[1], velocity)
    source_term = p.htrans_coeff() .* (e[1:end] .- p.T_ambient)

    @views de .= (convection_term .+ source_term)

    return nothing
end

function prosumer_edge!(de, e, v_s, v_d, p, _)
    # Prosumer edges must always have n_cells == 1
    # de[1] == 0, algebraic constraint

    de[1] = p.delta_T - (e[2] - e[1]) # Fixed temperature change across edge
    de[2] = 0.0 # Dummy

    return nothing
end

function junction_node!(dv, v, edges_in, edges_out, p, _)
    # DirectedODEVertex
    # dv[1] = 0.0

    # Calculate node temperature from incoming and outgoing edges
    # Assumption: edge state 1 => mass flow, edge states [2:end] => temperatures in finite-volume cells
    enthalpy_in = 0.0
    massflow_out = 0.0

    # enthalpy_in += sum(map(e -> e[1] * e[end], filter(e -> e[1] > 0, edges_in)))
    #     # for each edge in, if massflow is +ve (ie. massflow into node),
    #     #   enthalpy_in += massflow * temperature at edge-node interface
    # enthalpy_in += sum(map(e -> -e[1] * e[2], filter(e -> e[1] < 0, edges_out)))
    #     # for each edge out, if massflow is -ve (ie. massflow into node, since massflow direction is defined wrt. edge direction),
    #     #   enthalpy_in += (-massflow) * temperature at edge-node interface, - since massflow is -ve

    # massflow_out += sum(map(e -> e[1], filter(e -> e[1] > 0, edges_out)))
    # massflow_out += sum(map(e -> -e[1], filter(e -> e[1] < 0, edges_in)))
    
    # Remove after testing temperature-only case
    massflow = p.massflow
    for e in edges_in
        enthalpy_in += massflow * e[end]
    end
    for e in edges_out
        massflow_out += massflow
    end

    dv[1] = v[1] - (enthalpy_in / massflow_out) # node_temp = enthalpy_in / massflow_out

    return nothing
end

function fixed_node!(dv, v, _, _, p, _)
    # DirectedODEVertex
    # dv[1] = 0.0
    dv[1] = v[1] - p.T_fixed
    return nothing
end

edges::Vector{nd.ODEEdge} = [nd.ODEEdge(f=pipe_edge!, dim=n_cells[i], coupling=:directed,
                                        sym=[Symbol("T$j") for j in 1:n_cells[i]])
                            for i in 1:4]

# pushfirst!(edges, nd.ODEEdge(f=prosumer_edge!, dim=2, coupling=:directed,
#                              sym=[Symbol("T$j") for j in 1:2]))

nodes::Vector{nd.DirectedODEVertex} = [nd.DirectedODEVertex(f=junction_node!, dim=1,
                                                            mass_matrix=zeros(1,1), sym=[:T])
                                      for i in 2:4]

pushfirst!(nodes, nd.DirectedODEVertex(f=fixed_node!, dim=1, mass_matrix=zeros(1,1), sym=[:T]))


nd_fn = nd.network_dynamics(nodes, edges, g)

# Initialise solution
n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
initial_guess = zeros(n_states)
init_prob = de.SteadyStateProblem(nd_fn, initial_guess, params)
init_sol = de.solve(init_prob, de.DynamicSS(de.Rodas5()))


# display(fig_graph)
