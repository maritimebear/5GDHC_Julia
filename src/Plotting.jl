module Plotting

import Graphs as gr
import GraphMakie, GLMakie
import ColorSchemes as cs

export plot_hydrodynamics

# Functions to plot solution on graph, following "Stress on Truss" example:
# https://graph.makie.org/stable/generated/truss/

function plot_hydrodynamics(solution::Vector{Float64}, graph::gr.SimpleDiGraph)
    cscheme = cs.colorschemes[:delta] # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/
    ilabels = repr.(gr.nv(graph)) # Node-internal labels: node number
    elabels = ["Edge $idx\nMass flow rate: $m"
               for (idx, m) in zip(repr.(gr.ne(graph)), massflows)]


    plot = GraphMakie.graphplot(graph,
                                ilabels = ilabels, 
                                elabels = elabels,
                               )
end



end # module
