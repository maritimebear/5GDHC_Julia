module Fluids

export Fluid, Water, specific_heat

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
end


@inline function thermal_conductivity(::Type{Water}, temperature::Float64)
end


end # module
