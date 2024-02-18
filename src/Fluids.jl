module Fluids

export Fluid, Water # Types
export specific_heat, dynamic_viscosity, thermal_conductivity # Methods dispatching on types


abstract type Fluid end

struct Water <: Fluid end


# Thermophysical properties as functions of temperature
# Methods dispatch on Type{T <: Fluid}
# Calling convention followed in DynamicalFunctions: f(Type{<:Fluid}, temperature)

@inline function specific_heat(::Type{Water}, temperature::Float64)
    # -> Float64
    # Specific heat capacity [J/kg-K] according to PPDS equation:
    # VDI Heat Atlas (2010), section D.3.1

    A = 0.2399; B = 12.8647; C = -33.6392; D = 104.7686; E = -155.4709; F = 92.3726;
        # VDI Heat Atlas (2010), D3.1. Table 5 (pg. 336)
        # Values are for specific heat in [J/g-K], but specifying R in [J/kg-K]
        # effectively converts to [J/kg-K]

    T_c = 647.096 # Critical temperature [K]
    R = 461.5 # Specific gas constant [J/kg-K]

    t = 1.0 - (temperature/T_c)
    return R * (A/t + B + C*t + D*t^2 + E*t^3 + F*t^4)
end


@inline function dynamic_viscosity(::Type{Water}, temperature::Float64)
    # -> Float64
    # Dynamic viscosity [Pa-s] according to PPDS equation:
    # VDI Heat Atlas (2010), section D.3.1

    A = 0.45047; B = 1.39753; C = 613.181; D = 63.697; E = 0.00006896;
        # VDI Heat Atlas (2010), D3.1. Table 7 (pg. 352)
        # Table heading says values are for dynamic viscosity in [mPa-s],
        # but results after testing match known data in [Pa-s]

    K = (C - temperature) / (temperature - D)
    K1 = K < 0 ? -( (temperature - C) / (temperature - D) )^(1/3) : K^(1/3)
        # as per VDI Heat Atlas suggestion
    K2 = K1 * K
        # K2 = K^(4/3)

    return E * exp(A*K1 + B*K2)
end


@inline function thermal_conductivity(::Type{Water}, temperature::Float64)
    # -> Float64
    # Thermal conductivity [W/m-K] according to 4th-degree polynomial:
    # VDI Heat Atlas (2010), section D.3.1

    A = -2.4149; B = 2.45165e-2; C = -0.73121e-4; D = 0.99492e-7; E = -0.53730e-10;
        # VDI Heat Atlas (2010), D3.1. Table 9 (pg. 367)

    return A + B*temperature + C*temperature^2 + D*temperature^3 + E*temperature^4
end


end # module
