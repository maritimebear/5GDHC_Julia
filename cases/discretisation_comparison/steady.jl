# Script to compare steady-state solution for various discretisation schemes and grid sizings

# Sources for reference values:
#
# Hirsch and Nicolai, "An efficient numerical solution method for detailed modelling of large
#   5th generation district heating and cooling networks", 2022, Section 4.1, Case 1
#
# Cengel and Cimbala, "Fluid Mechanics: Fundamentals and Applications", 4th ed., table 8-2 (pg.371)
#
# Licklederer et al, "Thermohydraulic model of Smart Thermal Grids with bidirectional power flow
#   between prosumers", 2021


import Graphs as gr
import NetworkDynamics as nd
import DifferentialEquations as de
import SciMLNLSolve
import Random

import GLMakie as mk

include("../../src/DHG.jl")
import .DHG

include("./utils_discncomparison.jl")
import .utils

# TODO: actual values for pipe geometry, pump model, thermal controls

# --- Parameters ---- #

Random.seed!(93851203598)
solver = SciMLNLSolve.NLSolveJL
initial_dx = 10.0 # [m]
refinement_ratio = 2
n_refinement_levels = 6


## Fixed/reference values
n_nodes = 4
init_massflows = rand(Float64, (n_nodes, )) # Initial values for massflow states
T_ambient = 273.15 + 10.0 # [K]
p_ref = 0.0 # [Pa], reference pressure


## Pipe geometry, same for all pipes
# Values from Hirsch and Nicolai:
pipe_innerdiameter = 40.8e-3 # [m]
pipe_outerdiameter = 50e-3 # [m]
pipe_length = 100.0 # [m]
wall_conductivity = 0.4 # [W/m-K]


## Prosumers: massflow, thermal power
massflow = 0.3 # [kg/s], Hirsch and Nicolai
producer_heatrate = 1e3 # [W]
consumer_heatrate = -0.5e3 # [W]


## Pump model for producer
# Values from Licklederer et al:
pump_nominalspeed = 4100.0 # rpm
pump_ref1 = (0.0, 40221.0, pump_nominalspeed) # (massflow [kg/s], deltaP [Pa], speed [rpm])
pump_ref2 = (0.922, 0.0, pump_nominalspeed)


## Material properties, water as fluid
density = 1e3 # [kg/m^3]
fluid_T = DHG.Fluids.Water


## Properties of pipe wall, taking steel as wall material
wall_roughness = 0.045e-3 # [m], Cengel and Cimbala


## Discretisation
dxs = [initial_dx / (refinement_ratio ^ r) for r in 0:n_refinement_levels]

convection_schemes = Dict("Upwind" => DHG.Discretisation.upwind,
                         )


## Transport properties
transport_models = DHG.TransportModels(friction_factor=DHG.Transport.friction_Churchill,
                                       Nusselt_number=DHG.Transport.Nusselt_ChiltonColburn)


## Prosumer functions
consumer_hydctrl = (t) -> (massflow)
consumer_hydchar = (ctrl_input, massflow) -> (ctrl_input) # Return control input unchanged
consumer_thmctrl = (t) -> (consumer_heatrate)

producer_hydctrl = (t) -> (pump_nominalspeed)
producer_hydchar = DHG.Miscellaneous.PumpModel(pump_ref1..., pump_ref2..., density, pump_nominalspeed)
producer_thmctrl = (t) -> (producer_heatrate)


## Network structure
node_structs = (DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.JunctionNode(),
                DHG.ReferenceNode(p_ref),
               ) # Using tuples instead of Vectors: types are not uniform, network structure is constexpr

edge_structs = (
                DHG.Pipe(1, 3, # src, dst
                         pipe_innerdiameter, pipe_outerdiameter, pipe_length,
                         wall_roughness, wall_conductivity), # hot pipe
                DHG.PressureChange(2, 1,
                                   producer_hydctrl, producer_thmctrl, producer_hydchar), # producer
                DHG.Massflow(3, 4,
                             consumer_hydctrl, consumer_thmctrl, consumer_hydchar), # consumer
                DHG.Pipe(4, 2,
                         pipe_innerdiameter, pipe_outerdiameter, pipe_length,
                         wall_roughness, wall_conductivity), # cold pipe
               )

# --- end of parameters --- #


params = (density=density, T_ambient=T_ambient)

result_tuple_T = NamedTuple{(:syms, :sol),
                            Tuple{Vector{Symbol}, Vector{Float64}}
                           } # DataType: (syms=symbols vector, sol=solution vector)

results = Dict{String,                      # convection scheme name
               Vector{result_tuple_T}       # dx index => result
              }() # Using Vector and not Dict(dx => result) to not use Floats as keys

# Compare discretisation schemes
for (name, scheme) in convection_schemes

    println("Convection scheme: $name")
    results[name] = Vector{result_tuple_T}()

    for (idx_dx, dx) in enumerate(dxs)
        discretisation = DHG.Discretisation.FVM(dx=dx, convection=scheme)

        ## Set up problem
        nd_fn, g = DHG.assemble(node_structs, edge_structs, transport_models, discretisation, fluid_T)
            # nd_fn = nd.network_dynamics(collect(node_structs), collect(edge_structs), g)

        ## Initialise state vector: massflow states must be nonzero, required for node temperature calculation
        #   number of massflow and pressure states are constant; == n_nodes
        #   number of temperature states depends on dx
        initial_guess = DHG.Miscellaneous.initialise(nd_fn,
                                                     (_) -> init_massflows, # massflows to init_massflow
                                                     p_ref,                 # pressures to p_ref
                                                     T_ambient              # temperatures to T_ambient
                                                    )

        ## Solve for steady-state
        print("Starting solve: dx = $dx")
        push!(results[name],
              (syms=nd_fn.syms, sol=DHG.Miscellaneous.solve_steadystate(nd_fn, initial_guess, params, solver()))
             )
        println(" --- done")
    end
end

# Post-processing

node_Ts = [utils.get_states(DHG.PostProcessing.node_T_idxs, results, node_idx) |>
           dict -> Dict(k => [v[1] for v in vs] for (k, vs) in dict) # Unpack 1-element Vector{Float64}
           for (node_idx, _) in enumerate(node_structs)
          ] # Temperature at each node across discn. schemes and dxs

plot_x = [(i => dx) for (i, dx) in enumerate(dxs)]

## Plots

fig_nodeTs = mk.Figure()
axes_nodeTs = [mk.Axis(fig_nodeTs[row, 1],
                       xtickformat = xs -> ["1/$(2^Int(x))" for x in xs],
                       yticks=mk.LinearTicks(3),
                       # yminorticks=mk.IntervalsBetween(5), yminorticksvisible=true, yminorgridvisible=true,
                       # yticks = mk.MultiplesTicks(5, 1e-3, "a"),
                      ) for row in 1:4
              ]

lines_nodeTs = Dict(scheme_name => [mk.scatterlines!(axes_nodeTs[i],
                                                     1:length(dxs), node_Ts[i][scheme_name],
                                                     linestyle=:dash, marker=:circle,
                                                    )
                                    for (i, _) in enumerate(node_structs)
                                   ]
                    for (scheme_name, _) in convection_schemes
                   )

for (i, ax) in enumerate(axes_nodeTs)
    # Inset titles for axes
    mk.text!(ax,
             1, 1, text="Node $i", font=:bold, #fontsize=11,
             align=(:right, :top), space=:relative, offset=(-8, -4), justification=:right,
            )
end

# Figure title
axes_nodeTs[1].title = "Node temperatures vs. discretisation"

# Axis labels as text to offset
# mk.Label(fig_nodeTs[4, 2, mk.Bottom()], "test", halign=:left)
# mk.text!(fig_nodeTs.scene, 0, 1, space=:relative,
#          text="(cell size / base size)",
#          align=(:right, :bottom), offset=(0, 0), justification=:right,
#         )
