import Graphs as gr
import NetworkDynamics as nd
import LinearAlgebra as la
import DifferentialEquations as de
import GLMakie, GraphMakie
import ModelingToolkit as mtk

# Parameters
params = (density = 1e3,
          dyn_visc = 8.9e-4,
          diameter = 1.0,
          p_ref = 101325.0,
          delta_p = 1.0,
          massflow = 1.0
         )

friction_fn = (Re) -> (-1.0)


function pipe_edge!(de, e, v_s, v_d, p, _)
    # Calculate local variables
    area= 0.25 * pi * (p.diameter ^ 2)
    velocity= e[1] / (p.density * area)
    Re= p.density * velocity * p.diameter / p.dyn_visc

    # @show Re
    # @show velocity

    de[1] = p.density * area * ( (v_s[1] - v_d[1]) + (friction_fn(Re) * velocity * abs(velocity)) ) # Pressure drop due to friction
    # de[2] = e[2] + e[1]

    return nothing
end

function prosumer_edge!(de, e, v_s, v_d, p, _)
    de[1] = p.delta_p - (v_d[1] - v_s[1]) # Fixed pressure difference
    # de[1] = p.massflow - e[1] # Fixed mass flow rate
    # de[2] = e[2] + e[1]
    return nothing
end

function junction_node!(dv, _, edges_in, edges_out, _, _)
    # DirectedODEVertex
    # dv[1] = 0.0
    dv[1] = sum(map(e -> e[1], edges_in)) - sum(map(e -> e[1], edges_out)) # Mass conservation
    # dv[1] = sum(map(e -> e[1], edges_in)) - sum(map(e -> e[1], edges_out)) # Mass conservation
    return nothing
end

function fixed_pressure_node!(dv, v, _, _, p, _)
    # DirectedODEVertex
    dv[1] = v[1] - p.p_ref
    return nothing
end


# Graph
g = gr.cycle_digraph(4)

fig_graph = GraphMakie.graphplot(g; ilabels=repr.(1:gr.nv(g)), elabels=repr.(1:gr.ne(g)))

nodes::Vector{nd.DirectedODEVertex} = [nd.DirectedODEVertex(f=junction_node!, dim=1,
                                                          mass_matrix=zeros(1,1), sym=[:p])
                                      for _ in 1:3]
pushfirst!(nodes, nd.DirectedODEVertex(f=fixed_pressure_node!, dim=1, mass_matrix=zeros(1,1),
                                      sym=[:p]))

edges::Vector{nd.ODEEdge} = [nd.ODEEdge(f=pipe_edge!, dim=1, coupling=:directed, sym=[:m])
                             for _ in 1:3]
pushfirst!(edges, nd.ODEEdge(f=prosumer_edge!, dim=1, coupling=:directed,
                             mass_matrix=zeros(1,1), sym=[:m]))

nd_fn = nd.network_dynamics(nodes, edges, g)

# Initialise solution
n_states = sum([mapreduce(x -> x.dim, +, v) for v in (nodes, edges)])
sol_init = zeros(n_states)

full_prob = de.NonlinearProblem(nd_fn, sol_init, params)

# DAE index reduction
traced_sys = mtk.modelingtoolkitize(full_prob)
reduced_sys = mtk.structural_simplify(mtk.dae_index_lowering(traced_sys))
prob = de.SteadyStateProblem(reduced_sys, Pair[])

sol = de.solve(prob, de.DynamicSS(de.Rodas5()))


display(fig_graph)
