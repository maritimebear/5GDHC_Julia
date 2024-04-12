# Script to setup cases for performance benchmarks vs network size (varying the number of prosumers)
#
# References:
#
# Hirsch and Nicolai, "An efficient numerical solution method for detailed modelling of large
#   5th generation district heating and cooling networks", 2022, Section 4.1, Case 1
#
# Licklederer et al, "Thermohydraulic model of Smart Thermal Grids with bidirectional power flow
#   between prosumers", 2021
#
# Rocha et al, "Internal surface roughness of plastic pipes for irrigation", 2017
#
# Dang et al, "Fifth generation district heating and cooling: A comprehensive survey", 2024


import DifferentialEquations as de
import Random


# --- Parameters ---- #

Random.seed!(93851203598)
solver = de.Rodas5


## Time integration parameters
time_interval = (0.0, 1 * 60.0 * 60.0) # [s]
save_interval = 1 * 60.0 # [s]
save_times = time_interval[1]:save_interval:time_interval[2]
max_CFL = 1.0

solver_kwargs = Dict{Symbol, Any}()
solver_kwargs[:saveat] = collect(save_times)


## Benchmark parameters
benchmarkparameters = Dict(:samples => 10,
                           :seconds => 100.0,
                           :evals => 1,
                          )


## Fixed/reference values
# n_nodes = 4
# init_massflows = rand(Float64, (n_nodes, )) # Initial values for massflow states
T_ambient = 273.15 + 5.0 # [K], Hirsch
p_ref = 0.0 # [Pa], reference pressure


## Pipe geometry, same for all pipes
# Values from Hirsch and Nicolai:
pipe_innerdiameter = 40.8e-3 # [m]
pipe_outerdiameter = 50e-3 # [m]
pipe_length = 100.0 # [m]
wall_conductivity = 0.4 # [W/m-K]


## Pump model for producer
# Values from Licklederer et al:
pump_nominalspeed = 4100.0 # rpm
pump_ref1 = (0.0, 40221.0, pump_nominalspeed) # (massflow [kg/s], deltaP [Pa], speed [rpm])
pump_ref2 = (0.922, 0.0, pump_nominalspeed)




## Properties of pipe wall, polyethylene pipe [Hirsch]
wall_roughness = 8.116e-6 # [m], Rocha




# --- end of parameters --- #


params = (density=density, T_ambient=T_ambient)


result_tuple_T = NamedTuple{(:syms, :sol),
                            Tuple{Vector{Symbol}, Vector{Float64}}
                           } # DataType: (syms=symbols vector, sol=solution vector)

results = Dict{Int,                         # number of prosumers
               Vector{result_tuple_T}       # dx index => result
              }() # Using Vector and not Dict(dx => result) to not use Floats as keys

