# Script to plot variation of fluid thermophysical properties over temperature,
# to see if models behave as expected

import Plots as plt

include("../../src/Fluids.jl")

T_range = (0 + 273.15, 200 + 273.15) # [K]
dT = 1.0 # [K]

fluid = Fluids.PropyleneGlycol
density = 1e3 # [kg/m^3]


T = T_range[1] : dT : T_range[2]

struct Property
    name::String
    unit::String
    values::Vector{Float64}
end

names = ("Dynamic Viscosity", "Specific Heat Capacity", "Thermal Conductivity")
units = ("Pa-s", "J/kg-K", "W/m-K")
fns = (Fluids.dynamic_viscosity, Fluids.specific_heat, Fluids.thermal_conductivity)

# Calculate thermophysical properties over temperature range
properties = Dict(name => Property(name, unit, fn.(fluid, T))
                  for (name, unit, fn) in zip(names, units, fns)
                 )::Dict{String, Property}

plots = [plt.plot(T, prop.values,
                  label=name, legendposition=:inline, legendfontsize=8,
                 ) for (name, prop) in properties
        ]

fig = plt.plot(plots..., layout = (length(plots), 1), )

plt.xlabel!("Temperature [K]")
