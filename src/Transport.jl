module Transport # submodule, included in DHG.jl

export TransportModels
export Reynolds_number, Prandtl_number
export Nusselt_ChiltonColburn ,friction_Churchill


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


@inline function Nusselt_ChiltonColburn(friction_factor::Real, Reynolds_num::Real, Prandtl_num::Real)
    # -> Float64
    # Chilton-Colburn analogy: model for Nusselt number in turbulent internal flow, forced convection.
    # Cengel, "Heat Transfer: A Practical Approach", 2nd ed., section 8-6
    return 0.125 * friction_factor * Reynolds_num * (Prandtl_num^(1/3))
end


@inline function Nusselt_Gnielinsky(friction_factor::Real, Reynolds_num::Real, Prandtl_num::Real)
    # -> Float64
    # Model for Nusselt number in turbulent internal flow, forced convection.
    # Cengel, "Heat Transfer: A Practical Approach", 2nd ed., section 8-6, equation 8-70.
    f_8 = friction_factor / 8.0
    return (f_8 * (Reynolds_num - 1000.0) * Prandtl_num) / (1.0 + (12.7 * sqrt(f_8) * (Prandtl_num^(2/3) - 1.0)))
end


@inline function Nusselt_Colburn(_, Reynolds_num::Real, Prandtl_num::Real)
    # -> Float64
    # Colburn equation: model for Nusselt number in turbulent internal flow, forced convection.
    # Cengel, "Heat Transfer: A Practical Approach", 2nd ed., section 8-6, equation 8-67
    return 0.023 * (Reynolds_num^(0.8)) * (Prandtl_num^(1/3))
end


function friction_Churchill(Re::Real, relative_roughness::Real) # -> Float64
    # Re: Reynolds Number, relative_roughness = (mean pipe roughness height / pipe diameter) [m/m]
    # Approximate explicit expression for Darcy-Weisbach friction factor
    # Valid for laminar, transition and turbulent flows through circular sections of any roughness
    # Cengel and Cimbala, Fluid Mechanics: Fundamentals and Applications, 4th ed., equation 8-55
    inv_Re = 1.0 / Re # Re must be >0, not checking here for performance
    A = (-2.457 * log(((7.0 * inv_Re)^0.9) + (0.27 * relative_roughness)))^16
    B = (37530.0 * inv_Re)^16
    return 8.0 * (((8.0 * inv_Re)^12) + ((A + B)^(-3/2)))^(1/12)
end


end # (sub)module
