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
heat_capacity = 4184.0 # J/kg-K, used only in heat transfer coefficient calculation
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


edges::Vector{nd.ODEEdge} = [DHG.pipe_edge(pipelengths[1], dx[1]) for _ in 1:3]
pushfirst!(edges, DHG.prosumer_edge())

nodes::Vector{nd.DirectedODEVertex} = [DHG.junction_node() for _ in 1:3]
pushfirst!(nodes, DHG.fixed_node())

nd_fn = nd.network_dynamics(nodes, edges, g)

# Initialise solution
n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
initial_guess = ones(n_states)
init_prob = de.SteadyStateProblem(nd_fn, initial_guess, params)
init_sol = de.solve(init_prob, de.DynamicSS(de.Rodas5()))


# display(fig_graph)
