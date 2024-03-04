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
    src::IdxType                # source node index
    dst::IdxType                # destination node index
    inner_diameter::Float64     # [m]
    outer_diameter::Float64     # [m]
    length::Float64             # [m]
    roughness::Float64          # mean height of roughness [m]
    wall_conductivity::Float64  # thermal conductivity, [W/m-K]
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
