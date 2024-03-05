module PostProcessing

# Convenience functions for plotting etc.

import Graphs

import GLMakie, GraphMakie

export plot_graph


function plot_graph(g::Graphs.SimpleDiGraph)
    # -> handle to plot
    # Displays graph g with nodes and edges numbers (following implicit Graphs.jl ordering)
    return GraphMakie.graphplot(g; ilabels=repr.(1:Graphs.nv(g)), elabels=repr.(1:Graphs.ne(g)))
end





end # module
