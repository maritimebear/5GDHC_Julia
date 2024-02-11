module ParameterStruct

import SparseArrays as sp

export ProsumerPair, ProsumerParameters, Parameters


Base.@kwdef struct ProsumerPair{IndexType <: Integer}
    hydraulic::sp.SparseVector{Function, IndexType}
    thermal::sp.SparseVector{Function, IndexType}
end


Base.@kwdef struct ProsumerParameters{IndexType <: Integer}
    controls::ProsumerPair{IndexType}
    characteristics::ProsumerPair{IndexType}
end


Base.@kwdef mutable struct Parameters
    # mutable struct to allow modifications via callbacks
    density::Float64
    T_ambient::Float64
    p_ref::Float64
    T_fixed::Float64
    prosumers::ProsumerParameters{Int64}
end


end # module
