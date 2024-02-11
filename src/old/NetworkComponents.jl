module NetworkComponents # submodule, included in DHG.jl
# Structs to interface between parsed network and DHG.jl functions

export Node, Edge, ProsumerEdge, Junction, FixedNode, Pipe, PressureChangeProsumer, MassflowProsumer


abstract type Node end
abstract type Edge end
abstract type ProsumerEdge <: Edge end


struct Junction <: Node
    # No attributes for junction nodes
end


Base.@kwdef struct FixedNode <: Node
    pressure::Float64
    temperature::Float64
end


Base.@kwdef struct Pipe <: Edge
    diameter::Float64
    length::Float64
    dx::Float64
end


Base.@kwdef struct PressureChangeProsumer <: ProsumerEdge
    deltaP::Float64
    deltaT::Float64
end


Base.@kwdef struct MassflowProsumer <: ProsumerEdge
    massflow::Float64
    deltaT::Float64
end


end # (sub)module
