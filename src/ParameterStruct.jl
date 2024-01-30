module ParameterStruct

export Parameters

Base.@kwdef mutable struct Parameters
    # mutable struct to allow modifications via callbacks
    density::Float64
    T_ambient::Float64
    p_ref::Float64
    T_fixed::Float64
end

end
