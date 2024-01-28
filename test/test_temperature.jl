import Graphs as gr
import NetworkDynamics as nd
import DifferentialEquations as de

include("../src/DHG.jl")
import .DHG

# Script parameters
const inputfile = "./cycle4.gml"
const plot_graph = false

# Parameters
# Using material properties of water
const density = 1e3 # kg/m^3
T_ambient::Float64 = 273.15 + 10.0 # kelvin

# Properties of pipe wall
# Assuming steel as wall material
wall_conductivity = 30.0 # W/m-K, NOT A RELIABLE VALUE !!!
wall_thickness = 10e-3 # m

dyn_visc = () -> (8.9e-4) # Pa-s
friction = (Re) -> (-1.0)
heat_capacity = () -> 4184.0 # J/kg-K, used only in heat transfer coefficient calculation
h_wall = () -> (2 * wall_conductivity / wall_thickness) # Wikipedia: Heat transfer coefficient -- Heat transfer coefficient of pipe wall
heat_transfer = () -> (-h_wall() / (density * heat_capacity()))


@static if plot_graph; import GLMakie, GraphMakie; end


global_params = DHG.GlobalParameters(density=density, T_ambient=T_ambient)
transport_coeffs = DHG.TransportCoefficients(dynamic_viscosity=dyn_visc,
                                             wall_friction=friction,
                                             heat_capacity=heat_capacity,
                                             heat_transfer=heat_transfer)

dhg = DHG.DHGStruct(() -> DHG.parse_gml(inputfile),
                    global_params,
                    transport_coeffs
                   )


# Initialise solution

## Initial values must not be zero for initialiser to work
initial_guess = [1.0 for _ in 1:sum(reduce(+, v) for v in (dhg.n_states.nodes, dhg.n_states.edges))]
initialiser! = DHG.Utilities.initialiser(de.DynamicSS(de.Rodas5())) # Get closure function
init_retcode = initialiser!(dhg.f, initial_guess, dhg.parameters) # Call closure function, update initial_guess
# prob = de.ODEProblem(dhg.f, initial_guess, (0.0, 5*787.0), dhg.parameters)
# sol = de.solve(prob, de.Rodas5())


if plot_graph
    fig_graph = GraphMakie.graphplot(dhg.graph; ilabels=repr.(1:gr.nv(dhg.graph)), elabels=repr.(1:gr.ne(dhg.graph)))
    display(fig_graph)
end
