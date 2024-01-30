import Graphs as gr
import NetworkDynamics as nd
import DifferentialEquations as de

include("../src/DHG.jl")
import .DHG

const plot_graph = false
@static if plot_graph; import GLMakie, GraphMakie; end

# Parameters
p_ref = 101325.0 # Pa
T_fixed = 298.15 # K TODO: Remove T_fixed from ref_pressure node, calculate nodal temperature like junction nodes

# Using material properties of water
const density = 1e3 # kg/m^3
T_ambient::Float64 = 273.15 + 10.0 # kelvin
heat_capacity::Float64 = 4184.0 # J/kg-K, used only in heat transfer coefficient calculation
dyn_visc::Float64 = 8.9e-4 # Pa-s

# Properties of pipe wall
# Assuming steel as wall material
wall_conductivity = 30.0 # W/m-K, NOT A RELIABLE VALUE !!!
wall_thickness = 10e-3 # m

friction = (Re) -> (-1.0)
h_wall::Float64 = 2 * wall_conductivity / wall_thickness # Wikipedia: Heat transfer coefficient -- Heat transfer coefficient of pipe wall
heat_transfer::Float64 = -h_wall / (density * heat_capacity) # TODO: Why is this -ve?

transport_coeffs = DHG.TransportProperties(dynamic_viscosity=dyn_visc, wall_friction=friction,
                                             heat_capacity=heat_capacity, heat_transfer=heat_transfer)

params = DHG.Parameters(density=density, T_ambient=T_ambient, p_ref=p_ref, T_fixed=T_fixed)

