module utils_nsize

import Graphs as gr
import GLMakie, GraphMakie

include("../../src/DHG.jl")
import .DHG

export generate_graph, create_edge_struct, create_node_struct

function generate_graph(n_prosumers::Integer)
    # -> gr.SimpleDiGraph
    # Generates ladder-like graph with the specified number of prosumers
    #   i.e. each prosumer is connected in parallel, across two "rails" of pipes

    @assert iseven(n_prosumers)
    index_array = reshape([i for i in 1:(2 * n_prosumers)], 2, :)
    @assert size(index_array) == (2, n_prosumers)
    g = gr.SimpleDiGraph(2 * n_prosumers)
    
    # Create prosumer edges
    for col_idx in 1:n_prosumers
        row_idx = isodd(col_idx) ? (2, 1) : (1, 2) # odd column => producer, else consumer
        idxs = (index_array[row_idx[1], col_idx], index_array[row_idx[2], col_idx])
        gr.add_edge!(g, idxs...) ? nothing : throw("Failed to add edge: $idxs")
            # Not using @assert in case it gets elided; add_edge!() must run
    end

    # Create pipe edges
    for col_idx in 1:(n_prosumers - 1) # Hot rail
        idxs = (index_array[1, col_idx], index_array[1, (col_idx + 1)])
        gr.add_edge!(g, idxs...) ? nothing : throw("Failed to add edge: $idxs")
    end
    for col_idx in n_prosumers:-1:2 # Cold rail
        idxs = (index_array[2, col_idx], index_array[2, (col_idx - 1)])
        gr.add_edge!(g, idxs...) ? nothing : throw("Failed to add edge: $idxs")
    end

    return g
end


function create_edge_struct(src, dst,
                            innerdia, outerdia, length, roughness, conductivity,
                            prod_hydctrl, prod_thmctr, prod_hydchar,
                            cons_hydctrl, cons_thmctr, cons_hydchar,
                           )
    # -> <:DHG.Edge
    if dst == (src + 1) # consumer, all consumers are massflow prosumers in this case
        return DHG.Massflow(src, dst, cons_hydctrl, cons_thmctr, cons_hydchar)
    elseif dst == (src - 1) # producer, all producers are pressure change prosumers in this case
        return DHG.PressureChange(src, dst, prod_hydctrl, prod_thmctr, prod_hydchar)
    else # pipe
        return DHG.Pipe(src, dst, innerdia, outerdia, length, roughness, conductivity)
    end
    throw("This line is supposed to be unreachable")
end


function create_node_struct(node_idx, refnode_idx, refpressure)
    # -> <:DHG.Node
    # Stop Julia from complaining about type imports because of DHG import?
    return node_idx == refnode_idx ? DHG.ReferenceNode(refpressure) : DHG.JunctionNode()
end


end # module
