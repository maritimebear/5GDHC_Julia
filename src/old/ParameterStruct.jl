module ParameterStruct

import SparseArrays as sp

export ProsumerParameters, Parameters


Base.@kwdef struct ProsumerParameters{SparseVec1, SparseVec2, SparseVec3}
    controls_hydraulic::SparseVec1
    controls_thermal::SparseVec2
    characteristics_hydraulic::SparseVec3
end


Base.@kwdef mutable struct Parameters
    # mutable struct to allow modifications via callbacks
    density::Float64
    T_ambient::Float64
    p_ref::Float64
    T_fixed::Float64
    prosumers::ProsumerParameters
end


end # module
