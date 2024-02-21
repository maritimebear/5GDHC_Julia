module Transport # submodule, included in DHG.jl

export TransportModels
export Reynolds_number, Prandtl_number
export Nusselt_ChiltonCoburn ,friction_Churchill


Base.@kwdef struct TransportModels{T1, T2}
    # Members may be constants or callable types
    friction_factor::T1
    Nusselt_number::T2
end


@inline function Reynolds_number(velocity::Real, density::Real, diameter::Real, dynamic_viscosity::Real)
    # -> Float64
    return density * velocity * diameter / dynamic_viscosity
end


@inline function Prandtl_number(dynamic_viscosity::Real, isobaric_spec_heat::Real, thermal_conductivity::Real)
    # -> Float
    return dynamic_viscosity * isobaric_spec_heat / thermal_conductivity
end


@inline function Nusselt_ChiltonCoburn(friction_factor::Real, Reynolds_num::Real, Prandtl_num::Real)
    # -> Float64
    # Chilton-Coburn analogy: model for Nusselt number in turbulent internal flow, forced convection.
    # Cengel, "Heat Transfer: A Practical Approach", 2nd ed., section 8-6
    return 0.125 * friction_factor * Reynolds_num * (Prandtl_num^(1/3))
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
