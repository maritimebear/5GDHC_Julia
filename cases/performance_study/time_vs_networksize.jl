# Script to benchmark performance: computation time vs GCI for the Hirsch test case

import Graphs as gr
import Distributions as dbn

include("./setup_networksize.jl") # provides input parameters

include("./utils_networksize.jl")
import .utils_nsize

include("./utils_gci.jl")
import .utils as utils_gci

include("../../src/DHG.jl")
import .DHG

## Number of prosumers to run benchmarks on
n_prosumers = [2^i for i in 1:5]


## Discretisation
dx = 10.0 # [m]
convection_scheme = DHG.Discretisation.upwind
scheme_name = "Upwind"
discretisation = DHG.Discretisation.FVM(dx, convection_scheme)


## Fluid and Transport properties have to be defined in this script to keep Julia from complaining about type mismatches
# Material properties, taking propylene glycol (Propane-1,3-diol) as fluid [Hirsch]
density = 1064.4 # [kg/m^3], value at 0°C, VDI Heat Atlas D3.1. Table 2 (pg 315)
fluid_T = DHG.Fluids.PropyleneGlycol

# Transport properties
transport_models = DHG.TransportModels(friction_factor=DHG.Transport.friction_Churchill,
                                       Nusselt_number=DHG.Transport.Nusselt_ChiltonColburn)


## Prosumers: massflow, thermal power
massflow = 0.3 # [kg/s], Hirsch and Nicolai
consumer_heatrate = -2.7e3 # [W] Assuming temperature change across consumer = -4 K [Hirsch]

consumer_hydchar = (ctrl_input, massflow) -> (ctrl_input) # Return control input unchanged
producer_hydchar = DHG.Miscellaneous.PumpModel(pump_ref1..., pump_ref2..., density, pump_nominalspeed)

# Random variation in prosumer control functions
# Prosumer control functions are closures
consumer_hydctrl() = (t) -> (massflow * rand(dbn.Uniform(0.0, 1.0))) # massflow is scaled down to maintain CFL/maxdt limit
consumer_thmctrl() = (t) -> (consumer_heatrate * rand(dbn.Uniform(0.0, 2.0)))

producer_hydctrl() = (t) -> (pump_nominalspeed * rand(dbn.Uniform(0.0, 2.0)))
producer_thmctrl() = (t) -> (-consumer_heatrate * rand(dbn.Uniform(0.0, 2.0)))

expected_velocity = massflow / (density * 0.25 * pi * pipe_innerdiameter^2)

params = (density=density, T_ambient=T_ambient)

for n_prosumer in n_prosumers
    println("Number of prosumers: $n_prosumer")

    # Create Graphs.graph for easier edge struct creation
    g = utils_nsize.generate_graph(n_prosumer)
    println("Total pipe length: $((gr.ne(g) - n_prosumer) * pipe_length) m")

    node_structs = [DHG.ReferenceNode(p_ref); [DHG.JunctionNode() for _ in 2:(gr.nv(g))]...]
        # Concatenate vectors, 2:n_nodes since node 1 is set to reference node

    edge_structs = Vector{DHG.NetworkComponents.Edge}()
    for edge in gr.edges(g)
        if edge.dst == (edge.src + 1) # consumer, all consumers are massflow prosumers in this case
            push!(edge_structs, DHG.Massflow(edge.src, edge.dst,
                                             consumer_hydctrl(), consumer_thmctrl(), consumer_hydchar))
        elseif edge.dst == (edge.src - 1) # producer, all producers are pressure change prosumers in this case
            push!(edge_structs, DHG.PressureChange(edge.src, edge.dst,
                                                   producer_hydctrl(), producer_thmctrl(), producer_hydchar))
        else # pipe
            push!(edge_structs, DHG.Pipe(edge.src, edge.dst,
                                         pipe_innerdiameter, pipe_outerdiameter,
                                         pipe_length, wall_roughness, wall_conductivity))
        end
    end

    ## Set up problem
    nd_fn, _ = DHG.assemble(node_structs, edge_structs, transport_models, discretisation, fluid_T)
        # nd_fn = NetworkDynamics.network_dynamics(collect(node_structs), collect(edge_structs), g)

    ## Initialise state vector: massflow states must be nonzero, required for node temperature calculation
    initial_guess = DHG.Miscellaneous.initialise(nd_fn,
                                                 (_) -> rand(Float64, (gr.ne(g), )), # massflows to random #s
                                                 p_ref,                               # pressures to p_ref
                                                 T_ambient                            # temperatures to T_ambient
                                                )
    println("\nConvection scheme: $scheme_name, dx = $dx")

    if max_CFL !== nothing
        max_dt = max_CFL * dx / expected_velocity
        println("max dt = $max_dt")
        solver_kwargs[:dtmax] = max_dt
    end

    println("Starting benchmark\n")
    sol, benchmark = utils_gci.benchmark_solve(nd_fn, initial_guess, params, solver(), time_interval,
                                               solver_kwargs, benchmarkparameters
                                              )

    dump(benchmark)

    # push!(results[n_prosumer], (syms=nd_fn.syms, sol=sol.u[end]))
    
    println("\n --- done ---\n")
end
