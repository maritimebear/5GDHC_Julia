module NetworkComponents
# Structs to interface between parsed network and DHG.jl functions

export Node, Edge, ProsumerEdge, Junction, FixedNode, Pipe, PressureChangeProsumer, MassflowProsumer


abstract type Node end
abstract type Edge end
abstract type ProsumerEdge <: Edge end


struct Junction <: Node
    # No attributes for junction nodes
end


struct FixedNode <: Node
    pressure::Float64
    temperature::Float64
end


struct Pipe <: Edge
    # TODO
end


struct PressureChangeProsumer <: ProsumerEdge
    delta_p::Float64
    delta_T::Float64
end


struct MassflowProsumer <: ProsumerEdge
    massflow::Float64
end


end # module
