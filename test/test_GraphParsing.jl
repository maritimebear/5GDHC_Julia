import Graphs as gr

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
heat_transfer = () -> (-h_wall() / (density * heat_capacity))


@static if plot_graph; import GLMakie, GraphMakie; end


graph, node_dict, edge_dict = DHG.GraphParsing.parse_gml(inputfile)

if plot_graph
    fig_graph = GraphMakie.graphplot(graph; ilabels=repr.(1:gr.nv(graph)), elabels=repr.(1:gr.ne(graph)))
    display(fig_graph)
end


node_params = DHG.NodeParameters(node_dict)
edge_params = DHG.EdgeParameters(edge_dict)
global_params = DHG.GlobalParameters(density=density, T_ambient=T_ambient)
parameters = DHG.Parameters(global_params, node_params, edge_params)

transport_coeffs = DHG.TransportCoefficients(dynamic_viscosity=dyn_visc,
                                             wall_friction=friction,
                                             heat_capacity=heat_capacity,
                                             heat_transfer=heat_transfer)
