import Graphs as gr
import NetworkDynamics as nd
import DifferentialEquations as de
import GLMakie, GraphMakie

include("../../src/InterpolationSchemes.jl") # provides module FVM
include("../../src/DHG.jl")
import .DHG, .FVM

# Graph
g = gr.cycle_digraph(4)
fig_graph = GraphMakie.graphplot(g; ilabels=repr.(1:gr.nv(g)), elabels=repr.(1:gr.ne(g)))

# Parameters
# Using material properties of water
density = 1e3 # kg/m^3
heat_capacity = 4184.0, # J/kg-K, used only in heat transfer coefficient calculation
dyn_visc = 8.9e-4 # Pa-s

# Parameters for each edge
diameters = [1.0 for _ in 1:gr.ne(g)]
massflows = [1.0 for _ in 1:gr.ne(g)]
pipelengths = [1.0 for _ in 1:gr.ne(g)]
n_cells = [10 for _ in 1:gr.ne(g)] # number of cells in each edge
dx = pipelengths ./ n_cells # cell widths of each edge, not cell-centre locations

htrans_coeff = () -> (-1.0)
friction_fn = (Re) -> (-1.0)

# Named tuple to be passed to dynamical functions
params = (density = density,
          dyn_visc = dyn_visc,
          diameter = diameters[1],
          massflow = massflows[1],
          dx = dx[1],
          delta_p = 1.0,
          delta_T = 0.0,
          p_ref = 101325.0,
          T_fixed = 273.15 + 25.0,
          friction_fn = friction_fn,
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
    area = 0.25 * pi * (p.diameter ^ 2)
    velocity = e[1] / (p.density * area)
    Re = p.density * velocity * p.diameter / p.dyn_visc

    # e[1] : mass flow rate, algebraic constraint
    #   => de[1] == 0, used to calculate pressure drop across pipe due to friction
    # TODO: Implement pressure loss according to Cengel eqn. 8-21,
    #       Churchill or Swamee-Jain approximation for Darcy-Weisbach f
    de[1] = (v_s[1] - v_d[1]) + (p.friction_fn(Re) * velocity * abs(velocity)) # Pressure drop due to friction

    @views de[2:end] .= -(1 / p.dx) .* FVM.upwind(e[2:end], v_s[2], v_d[2], velocity)

    return nothing
end

function prosumer_edge!(de, e, v_s, v_d, p, _)
    # Prosumer edges must always have dims == 2
    # de[1:2] == 0, algebraic constraint

    de[1] = p.massflow - e[1] # Fixed mass flow rate
    de[2] = p.delta_T - e[2] # Fixed temperature change across edge

    return nothing
end

function junction_node!(dv, v, edges_in, edges_out, _, _)
    # DirectedODEVertex, dims == 2
    # dv[1:2] = 0.0
    dv[1] = sum(map(e -> e[1], edges_in)) - sum(map(e -> e[1], edges_out)) # Mass conservation

    # Calculate node temperature from incoming and outgoing edges
    # Assumption: edge state 1 => mass flow, edge states [2:end] => temperatures in finite-volume cells
    enthalpy_in = 0.0
    massflow_out = 0.0

    enthalpy_in += sum(map(e -> e[1] * e[end], filter(e -> e[1] > 0, edges_in)))
        # for each edge in, if massflow is +ve (ie. massflow into node),
        #   enthalpy_in += massflow * temperature at edge-node interface
    enthalpy_in += sum(map(e -> -e[1] * e[2], filter(e -> e[1] < 0, edges_out)))
        # for each edge out, if massflow is -ve (ie. massflow into node, since massflow direction is defined wrt. edge direction),
        #   enthalpy_in += (-massflow) * temperature at edge-node interface, - since massflow is -ve

    massflow_out += sum(map(e -> e[1], filter(e -> e[1] > 0, edges_out)))
    massflow_out += sum(map(e -> -e[1], filter(e -> e[1] < 0, edges_in)))

    dv[2] = v[2] - (enthalpy_in / massflow_out) # node_temp = enthalpy_in / massflow_out

    return nothing
end

function fixed_node!(dv, v, _, _, p, _)
    # DirectedODEVertex, dims == 2
    # dv[1:2] = 0.0
    dv[1] = v[1] - p.p_ref
    dv[2] = v[2] - p.T_ref
    return nothing
end

edges::Vector{nd.ODEEdge} = [DHG.pipe_edge(pipelengths[1], dx[1]) for _ in 1:3]
pushfirst!(edges, DHG.prosumer_edge())

nodes::Vector{nd.DirectedODEVertex} = [DHG.junction_node() for _ in 1:3]
pushfirst!(nodes, DHG.fixed_node())

nd_fn = nd.network_dynamics(nodes, edges, g)

# Initialise solution
n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
initial_guess = zeros(n_states)
init_prob = de.SteadyStateProblem(nd_fn, initial_guess, params)
init_sol = de.solve(init_prob, de.DynamicSS(de.Rodas5()))


# display(fig_graph)
