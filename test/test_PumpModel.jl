# Script to test Miscellaneous.jl::PumpModel
# Generate test domain |> Closure from PumpModel |> plot

# PumpModel and reference points from Licklederer et al,
# "Thermohydraulic model of Smart Thermal Grids with bidirectional
# power flow between prosumers", 2021, pages 4, 12.

import Plots as plt

include("../src/Miscellaneous.jl")
import .Miscellaneous: PumpModel

Base.@kwdef struct ReferencePoint
    m::Float64 # mass flow rate [kg/s]
    dP::Float64 # pressure change [Pa]
    n::Float64 # pump speed [rpm]
end

n_nominal = 4100.0 # nominal pump speed [rpm]
fluid_density = 1e3 # [kg/s]

ref_1 = ReferencePoint(m=0.0, dP=40221.0, n=1.0 * n_nominal)
ref_2 = ReferencePoint(m=0.922, dP=0.0, n=1.0 * n_nominal)

model = PumpModel(ref_1.m, ref_1.dP, ref_1.n,
                  ref_2.m, ref_2.dP, ref_2.n,
                  fluid_density, n_nominal)

dx = 1e-2
x = -ref_2.m:dx:ref_2.m # Test domain symmetric about 0, x = massflow
y = model.(n_nominal, x)

plt.plot(x, y)
plt.xlabel!("massflow [kg/s]")
plt.ylabel!("deltaP [Pa]")
