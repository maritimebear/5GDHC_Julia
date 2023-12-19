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
heat_capacity = 4184.0, # J/kg-K

# Parameters for each edge
diameters = [1.0 for _ in 1:gr.ne(g)]
massflows = [1.0 for _ in 1:gr.ne(g)]
pipelengths = [1.0 for _ in 1:gr.ne(g)]
n_cells = [10 for _ in 1:gr.ne(g)] # number of cells in each edge
dx = pipelengths ./ n_cells # cell widths of each edge, not cell-centre locations

htrans_coeff = (_) -> (-1.0)

# Named tuple to be passed to dynamical functions
params = (density = density,
          heat_capacity = heat_capacity,
          diameter = diameters[1],
          massflow = massflows[1],
          dx = dx[1],
          delta_T = 5.0,
          T_fixed = 273.15 + 25.0,
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

    de .= -(1 / p.dx) .* FVM.upwind(e, v_s[1], v_d[1], velocity)

    return nothing
end

function prosumer_edge!(de, e, v_s, v_d, p, _)
    # Prosumer edges must always have dim == 1
    # de[1] == 0, algebraic constraint

    de[1] = p.delta_T - e[1] # Fixed temperature change across edge

    return nothing
end

function junction_node!(dv, _, edges_in, edges_out, _, _)
    # DirectedODEVertex
    # dv[1] = 0.0
    dv[1] = sum(map(e -> e[1], edges_in)) - sum(map(e -> e[1], edges_out)) # Mass conservation
    return nothing
end

function fixed_pressure_node!(dv, v, _, _, p, _)
    # DirectedODEVertex
    dv[1] = v[1] - p.p_ref
    return nothing
end



nodes::Vector{nd.DirectedODEVertex} = [nd.DirectedODEVertex(f=junction_node!, dim=1,
                                                          mass_matrix=zeros(1,1), sym=[:p])
                                      for _ in 1:3]
pushfirst!(nodes, nd.DirectedODEVertex(f=fixed_pressure_node!, dim=1, mass_matrix=zeros(1,1),
                                      sym=[:p]))

edges::Vector{nd.ODEEdge} = [nd.ODEEdge(f=pipe_edge!, dim=1, coupling=:directed, mass_matrix=zeros(1,1), sym=[:m])
                             for _ in 1:3]
pushfirst!(edges, nd.ODEEdge(f=prosumer_edge!, dim=1, coupling=:directed,
                             mass_matrix=zeros(1,1), sym=[:m]))

nd_fn = nd.network_dynamics(nodes, edges, g)

# Initialise solution
n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
initial_guess = zeros(n_states)
init_prob = de.SteadyStateProblem(nd_fn, initial_guess, params)
init_sol = de.solve(init_prob, de.DynamicSS(de.Rodas5()))


# display(fig_graph)
