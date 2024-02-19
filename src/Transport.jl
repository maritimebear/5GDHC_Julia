module Transport # submodule, included in DHG.jl

export TransportProperties


Base.@kwdef struct TransportProperties{T1, T2}
    # Members may be constants or callable types
    wall_friction::T1
    heat_transfer::T2
end


function friction_Churchill(Re::Real, relative_roughness::Real) # -> Float64
    # Re: Reynolds Number, relative_roughness = (mean pipe roughness height / pipe diameter) [m/m]
    # Approximate explicit expression for Darcy-Weisbach friction factor
    # Valid for laminar, transition and turbulent flows through circular sections of any roughness
    # Cengel and Cimbala, Fluid Mechanics: Fundamentals and Applications, 4th ed., equation 8-55
    inv_Re = 1.0 / Re
    A = (-2.457 * log(((7.0 * inv_Re)^0.9) + (0.27 * relative_roughness)))^16
    B = (37530.0 * inv_Re)^16
    return 8.0 * (((8.0 * inv_Re)^12) + ((A + B)^(-3/2)))^(1/12)
end


end # (sub)module
