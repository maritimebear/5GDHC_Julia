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


Base.@kwdef struct Pipe{IdxType <: Integer} <: Edge
    src::IdxType
    dst::IdxType
    diameter::Float64
    length::Float64
    dx::Float64
end


Base.@kwdef struct PressureChange{IdxType <: Integer,
                                  HydCtrl_T <: Function,
                                  ThmCtrl_T <: Function,
                                  HydChr_T <: Function
                                 } <: Prosumer
    src::IdxType
    dst::IdxType
    hydraulic_control::HydCtrl_T
    thermal_control::ThmCtrl_T
    hydraulic_characteristic::HydChr_T
end


Base.@kwdef struct Massflow{IdxType <: Integer,
                            HydCtrl_T <: Function,
                            ThmCtrl_T <: Function,
                            HydChr_T <: Function
                           } <: Prosumer
    src::IdxType
    dst::IdxType
    hydraulic_control::HydCtrl_T
    thermal_control::ThmCtrl_T
    hydraulic_characteristic::HydChr_T
end

end # module
