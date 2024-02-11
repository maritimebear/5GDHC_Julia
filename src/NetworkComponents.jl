module NetworkComponents
# Structs to hold parameters for each type of node and edge in the network

export Node, Edge, Prosumer, JunctionNode, ReferenceNode, Pipe, PressureChange, Massflow


abstract type Node end
abstract type Edge end
abstract type Prosumer <: Edge end


Base.@kwdef struct JunctionNode <: Node
    # No parameters
end


Base.@kwdef struct ReferenceNode <: Node
    pressure::Float64
end


Base.@kwdef struct Pipe{srcNode_T <: Node, dstNode_T <: Node} <: Edge
    src::srcNode_T
    dst::dstNode_T
    diameter::Float64
    length::Float64
    dx::Float64
end


Base.@kwdef struct PressureChange{src_T <: Node,
                                  dst_T <: Node,
                                  HydCtrl_T <: Function,
                                  ThmCtrl_T <: Function,
                                  HydChr_T <: Function
                                 } <: Prosumer
    src::src_T
    dst::dst_T
    hydraulic_control::HydCtrl_T
    thermal_control::ThmCtrl_T
    hydraulic_characteristic::HydChr_T
end


Base.@kwdef struct Massflow{src_T <: Node,
                                  dst_T <: Node,
                                  HydCtrl_T <: Function,
                                  ThmCtrl_T <: Function,
                                  HydChr_T <: Function
                                 } <: Prosumer
    src::src_T
    dst::dst_T
    hydraulic_control::HydCtrl_T
    thermal_control::ThmCtrl_T
    hydraulic_characteristic::HydChr_T
end

end # module
