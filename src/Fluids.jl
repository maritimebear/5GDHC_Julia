module Fluids

export Fluid, Water, PropyleneGlycol # Types
export specific_heat, dynamic_viscosity, thermal_conductivity # Methods dispatching on types


abstract type Fluid end

struct Water <: Fluid end
struct PropyleneGlycol <: Fluid end # Propane-1,3-diol


# Thermophysical properties as functions of temperature, ignoring pressure dependence
# Methods dispatch on Type{T <: Fluid}
# Calling convention followed in DynamicalFunctions: f(Type{<:Fluid}, temperature)

@inline function specific_heat(::Type{Water}, temperature)
    # -> Float64
    # Specific heat capacity [J/kg-K] according to PPDS equation:
    # VDI Heat Atlas (2010), section D.3.1

    T = clamp(temperature, 273.15, 250.0 + 273.15)
        # Clamped to limits from VDI Heat Atlas (2010), D3.1. Table 5 (pg. 336)

    A = 0.2399; B = 12.8647; C = -33.6392; D = 104.7686; E = -155.4709; F = 92.3726;
        # VDI Heat Atlas (2010), D3.1. Table 5 (pg. 336)
        # Values are for specific heat in [J/g-K], but specifying R in [J/kg-K]
        # effectively converts to [J/kg-K]

    T_c = 647.10 # Critical temperature [K], VDI Heat Atlas D3.1. Table 1 (pg 302)
    R = 461.402 # Specific gas constant [J/kg-K], R = R_universal [J/mol-K] / Molar mass [kg/mol]

    t = 1.0 - (T/T_c)
    return R * (A/t + B + C*t + D*t^2 + E*t^3 + F*t^4)
end


@inline function specific_heat(::Type{PropyleneGlycol}, temperature)
    # -> Float64
    # Specific heat capacity [J/kg-K] according to PPDS equation:
    # VDI Heat Atlas (2010), section D.3.1

    T = clamp(temperature, -25.0 + 273.15, 250.0 + 273.15)
        # Clamped to limits from VDI Heat Atlas (2010), D3.1. Table 5 (pg. 341)

    B = 41.1843; C = -33.4694; # A = D = E = F = 0
        # VDI Heat Atlas (2010), D3.1. Table 5 (pg. 341)
        # Values are for specific heat in [J/g-K], but specifying R in [J/kg-K]
        # effectively converts to [J/kg-K]

    T_c = 724.05 # Critical temperature [K], VDI Heat Atlas D3.1. Table 1 (pg 306)
    R = 109.271 # Specific gas constant [J/kg-K], R = R_universal [J/mol-K] / Molar mass [kg/mol]

    t = 1.0 - (T/T_c)
    return R * (B + C*t)
end


@inline function dynamic_viscosity(::Type{Water}, temperature)
    # -> Float64
    # Dynamic viscosity [Pa-s] according to PPDS equation:
    # VDI Heat Atlas (2010), section D.3.1

    T = clamp(temperature, 273.15, 200.0 + 273.15)
        # Lower limit set to 0°C as this is within the expected temperature range for 5GDHCs,
        # upper limit from VDI Heat Atlas (2010), D3.1. Table 7 (pg. 352)

    A = 0.45047; B = 1.39753; C = 613.181; D = 63.697; E = 0.00006896;
        # VDI Heat Atlas (2010), D3.1. Table 7 (pg. 352)
        # Table heading says values are for dynamic viscosity in [mPa-s],
        # but results after testing match known data in [Pa-s]

    K = (C - T) / (T - D)
    K1 = K < 0 ? -( (T - C) / (T - D) )^(1/3) : K^(1/3)
        # as per VDI Heat Atlas suggestion
    K2 = K1 * K
        # K2 = K^(4/3)

    return E * exp(A*K1 + B*K2)
end


@inline function dynamic_viscosity(::Type{PropyleneGlycol}, temperature)
    # -> Float64
    # Dynamic viscosity [Pa-s] according to PPDS equation:
    # VDI Heat Atlas (2010), section D.3.1

    T = clamp(temperature, 50.0 + 273.15, 150.0 + 273.15)
        # Limits from VDI Heat Atlas (2010), D3.1. Table 7 (pg. 357)

    A = 4.55501; B = 0.75893; C = 711.876; D = 132.080; E = 0.00000695;
        # VDI Heat Atlas (2010), D3.1. Table 7 (pg. 357)
        # Table heading says values are for dynamic viscosity in [mPa-s],
        # but results after testing match known data in [Pa-s]

    K = (C - T) / (T - D)
    K1 = K < 0 ? -( (T - C) / (T - D) )^(1/3) : K^(1/3)
        # as per VDI Heat Atlas suggestion
    K2 = K1 * K
        # K2 = K^(4/3)

    return E * exp(A*K1 + B*K2)
end


@inline function thermal_conductivity(::Type{Water}, temperature)
    # -> Float64
    # Thermal conductivity [W/m-K] according to 4th-degree polynomial:
    # VDI Heat Atlas (2010), section D.3.1

    T = clamp(temperature, 273.15, 200.0 + 273.15)
        # Clamped to limits from VDI Heat Atlas (2010), D3.1. Table 9 (pg. 367)

    A = -2.4149; B = 2.45165e-2; C = -0.73121e-4; D = 0.99492e-7; E = -0.53730e-10;
        # VDI Heat Atlas (2010), D3.1. Table 9 (pg. 367)

    return A + B*T + C*T^2 + D*T^3 + E*T^4
end


@inline function thermal_conductivity(::Type{PropyleneGlycol}, temperature)
    # -> Float64
    # Thermal conductivity [W/m-K] according to 4th-degree polynomial:
    # VDI Heat Atlas (2010), section D.3.1

    T = clamp(temperature, -25.0 + 273.15, 200.0 + 273.15)
        # Clamped to limits from VDI Heat Atlas (2010), D3.1. Table 9 (pg. 373)

    A = 0.0867; B = 0.0667e-2; C = -0.00281e-4; D = -0.01797e-7; E = 0.01228e-10;
        # VDI Heat Atlas (2010), D3.1. Table 9 (pg. 373)

    return A + B*T + C*T^2 + D*T^3 + E*T^4
end


end # module
