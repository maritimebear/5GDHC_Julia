const plot_graph = true

import Graphs as gr

include("../src/GraphParsing.jl")
import .GraphParsing

@static if plot_graph; import GLMakie, GraphMakie; end

inputfile = "./cycle4.gml"

graph, node_dict, edge_dict = GraphParsing.parse_gml(inputfile)

if plot_graph
    fig_graph = GraphMakie.graphplot(graph; ilabels=repr.(1:gr.nv(graph)), elabels=repr.(1:gr.ne(graph)))
    display(fig_graph)
end



