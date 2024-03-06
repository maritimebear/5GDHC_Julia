# Script to compare models for Nusselt number
#
# Mass flow rate and pipe diameter_inner from Hirsch and Nicolai, "An efficient numerical solution
# method for detailed modelling of large 5th generation district heating and cooling networks",
# 2022, Section 4.1, Case 1.

import Plots as plt

include("../src/Transport.jl")
include("../src/Fluids.jl")

T_range = (0 + 273.15, 40 + 273.15) # [K]
    # Temperature range of 5GDHCs: 0°C to 35°C, Buffa et al, 2019
dT = 1.0 # [K]

fluid = Fluids.Water
density = 1e3 # [kg/m^3]
massflow = 0.3 # [kg/s], Hirsch and Nicolai

diameter_inner = 40.8e-3 # [m], Hirsch and Nicolai

wall_roughness = 0.045e-3 # [m] Cengel and Cimbala, "Fluid Mechanics: Fundamentals and Applications", 4th ed., table 8-2 (pg.371)
rel_roughness = wall_roughness / diameter_inner

velocity = massflow / (density * 0.25 * pi * diameter_inner^2)

T = T_range[1] : dT : T_range[2]

# Calculate transport properties over temperature range
mu = Fluids.dynamic_viscosity.(fluid, T) # [Pa-s]
Cp = Fluids.specific_heat.(fluid, T) # [J/kg-K]
k = Fluids.thermal_conductivity.(fluid, T) # [W/m-K]

Re = Transport.Reynolds_number.(velocity, density, diameter_inner, mu)
Pr = Transport.Prandtl_number.(mu, Cp, k)

friction = Transport.friction_Churchill.(Re, rel_roughness) # Darcy-Weisbach friction factor

# Nusselt number by different models/approximations
Nu = Dict("Chilton-Colburn" => Transport.Nusselt_ChiltonColburn.(friction, Re, Pr),
          "Colburn" => Transport.Nusselt_Colburn.(friction, Re, Pr),
          "Gnielinsky" => Transport.Nusselt_Gnielinsky.(friction, Re, Pr),
         )::Dict{String, Vector{Float64}}

for (name, _Nu) in Nu
    plt.plot!(T, _Nu, label=name,
              legendposition=:inline,
              legendtitle="Nusselt number model", legendfontsize=8,
             )
end
plt.xlabel!("Temperature [K]")
plt.ylabel!("Nusselt number")
plt.title!("Comparison of Nusselt number approximations")
