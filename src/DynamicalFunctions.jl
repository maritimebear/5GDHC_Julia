module DynamicalFunctions # submodule, included in DHG.jl
# Dynamical functions for NetworkDynamics ODEEdge and DirectedODEVertex

export pipe, prosumer_massflow, prosumer_deltaP, junction!, fixed_node!


import ..FVM
import ..Transport: TransportProperties


# TODO: Type annotations for local parameters? (eg. density::T where T is captured from outer function?)
#           # Probably better to annotate p::ParameterStructs.Parameters{ValueType, IndexType} in function signature
# TODO: Refactor junction!() with pipes to chain filter |> map |> sum
# TODO: Cleanup unused code


# Prosumer edges

@inline function prosumer_outlet_T(thermal_power::Real, massflow::Real, temperature_in::Real, spec_heat::Real) 
    # -> T_out = Q / (m * Cp) + T_in
    # Requires Q, m, Cp and T_in to be known
    return thermal_power / (massflow + spec_heat) + temperature_in
end


function PumpModel(massflow_ref1, deltaP_ref1, speed_ref1,
                    massflow_ref2, deltaP_ref2, speed_ref2,
                    fluid_density, nominal_speed
                )
    # Returns closure with reference values captured
    # Following Licklederer et al, "Thermohydraulic model of Smart Thermal Grids with bidirectional power flow between prosumers", 2021, pages 4, 12
    u_ref1 = speed_ref1 / nominal_speed
    u_ref2 = speed_ref2 / nominal_speed 
    vol_ref1 = massflow_ref1 / fluid_density
    vol_ref2 = massflow_ref2 / fluid_density

    c1 = (deltaP_ref1 - (u_ref1 / u_ref2)^2 * deltaP_ref2) / (vol_ref1^2 - (u_ref1 / u_ref2)^2 * vol_ref2^2)
    c2 = (1.0 / u_ref2)^2 * (deltaP_ref2 - (vol_ref2^2 * c1))

    function model(massflow, pump_speed)
        # Returns pressure increase across pump
        # Constant outer values in let-block for performance
        let c1 = c1
            c2 = c2
            nominal_speed = nominal_speed
            fluid_density = fluid_density
            # Local variables
            speed_ratio = pump_speed / nominal_speed
            vol_rate = massflow / fluid_density
            return (c1 * vol_rate^2) + (c2 * speed_ratio^2)
        end
    end

    return model
end


@inline @inbounds function prosumerstate_thermal(index, de, e, v_s, v_d, p, t)
    # -> Nothing
    # Sets thermal state in de, intended to be called by prosumer dynamic functions

    inlet_T = e[1] >= 0 ? v_s[1] : v_d[1] # Upwind convection
    control_input = p.prosumers.controls.thermal[index](t) # thermal power
    state_old = e[2] # Old outlet temperature
    state_new = prosumer_outlet_T(control_input, e[1], inlet_T, p.transport_properties.spec_heat)

    de[2] = state_new - state_old
    return nothing
end


@inline @inbounds function prosumerstate_fixdeltaP(index, de, e, v_s, v_d, p, t)
    # -> Nothing
    # Sets hydraulic state by fixing pressure change across prosumer,
    # intended to be called by prosumer dynamic functions

    control_input = p.prosumers.controls_hydraulic[index](t) # Pump speed
    state_old = v_d[1] - v_s[1]
    state_new = p.prosumers.characteristics.hydraulic[index](e[1], control_input) # new deltaP

    de[1] = state_new - state_old
    return nothing
end


@inline @inbounds function prosumerstate_fixmassflow(index, de, e, v_s, v_d, p, t)
    # -> Nothing
    # Sets hydraulic state by fixing mass flow rate across prosumer,
    # intended to be called by prosumer dynamic functions

    # control_input = p.prosumers.controls_hydraulic[index](t) # Not implemented, model valve opening?
    state_old = e[1]
    state_new = p.prosumers.characteristics.hydraulic[index](t) # new mass flow rate

    de[1] = state_new - state_old
    return nothing
end


## Each edge can have different parameters, so each edge function is a closure with an 'index' captured.
## 'index' is used to access corresponding values/functions from parameters 'p'

function prosumer_deltaP(index::Integer)
    # Returns closure with 'index' captured
    function f!(de, e, v_s, v_d, p, t)
        # Closure, implements physics for prosumer edges with fixed pressure change
        let
            index = index

            # Prosumer edges must always have dims == 3
            # de[1:3] == 0, algebraic constraint

            # Physics implementation
            prosumerstate_fixdeltaP(index, de, e, v_s, v_d, p, t) # Set pressure difference
            prosumerstate_thermal(index, de, e, v_s, v_d, p, t) # Set prosumer outlet temperature

            return nothing
        end # let block
    end # f!(...)
    return f!
end


function prosumer_massflow(index::Integer)
    # Returns closure with 'index' captured
    function f!(de, e, v_s, v_d, p, t)
        # Closure, implements physics for prosumer edges with fixed pressure change
        let
            index = index

            # Prosumer edges must always have dims == 3
            # de[1:3] == 0, algebraic constraint

            # Physics implementation
            prosumerstate_fixmassflow(index, de, e, v_s, v_d, p, t) # Set mass flow rate
            prosumerstate_thermal(index, de, e, v_s, v_d, p, t) # Set prosumer outlet temperature

            return nothing
        end # let block
    end # f!(...)
    return f!
end


# Pipe edges

function pipe(diameter::Float64, dx::Float64, coeff_fns::TransportProperties{F1, F2, F3, F4}) where {F1, F2, F3, F4}
    # Returns closure with 'index' and transport property functions captured

    function f!(de, e, v_s, v_d, p, _)
    # Closure, implements physics for pipe edges: wall friction, heat loss to environment
        let
            # Capture constants in let-block for performance
            # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured
            diameter = diameter
            dx = dx
            dyn_visc = coeff_fns.dynamic_viscosity # TODO: Performance: type annotate lhs?
            friction = coeff_fns.wall_friction
            htrans_coeff = coeff_fns.heat_transfer

            # Constant
            area = 0.25 * pi * (diameter ^ 2) # TODO: Mark const?

            # Get local parameters from Parameters struct
            density = p.density
            T_ambient = p.T_ambient
            # Calculate local variables
            velocity = e[1] / (density * area)
            Re = density * velocity * diameter / dyn_visc
            ## Momentum equation
            deltaP = friction(Re) * sign(velocity) * (velocity^2) # Pressure drop due to friction
            ## Energy equation
            @views convection = -(1 / dx) .* FVM.upwind(e[2:end], v_s[2], v_d[2], velocity)
            @views source = htrans_coeff .* (e[2:end] .- T_ambient)

            # Physics implementation
            # e[1] : mass flow rate, algebraic constraint
            # e[2:end] : temperatures in finite-volume cells
            #   => de[1] == 0, used to calculate pressure drop across pipe due to friction
            # TODO: Implement pressure loss according to Cengel eqn. 8-21,
            #       Churchill or Swamee-Jain approximation for Darcy-Weisbach f
            de[1] = deltaP - (v_d[1] - v_s[1]) # Momentum equation
            de[2:end] .= convection .+ source # Energy equation

            return nothing
        end # let block
    end # f!(...)

    return f!
end

# Dynamical functions for nodes

## No closures for node functions, as 'index' is not needed for node parameters

function junction!(dv, v, edges_in, edges_out, _, _)
    # DirectedODEVertex, dims == 2
    # dv[1:2] = 0.0

    # Calculate node temperature from incoming and outgoing edges
    # Assumption: edge state 1 => mass flow, edge states [2:end] => temperatures in finite-volume cells
    enthalpy_in = 0.0
    massflow_out = 0.0

    enthalpy_in += sum(map(e -> e[1] * e[end], filter(e -> e[1] > 0, edges_in)))
        # for each edge in, if massflow is +ve (ie. massflow into node),
        #   enthalpy_in += massflow * temperature at edge-node interface
    enthalpy_in += sum(map(e -> -e[1] * e[2], filter(e -> e[1] < 0, edges_out)))
        # for each edge out, if massflow is -ve (ie. massflow into node, since massflow direction is defined wrt. edge direction),
        #   enthalpy_in += (-massflow) * temperature at edge-node interface, - since massflow is -ve

    massflow_out += sum(map(e -> e[1], filter(e -> e[1] > 0, edges_out)))
    massflow_out += sum(map(e -> -e[1], filter(e -> e[1] < 0, edges_in)))

    # Physics implementation
    dv[1] = sum(map(e -> e[1], edges_in)) - sum(map(e -> e[1], edges_out)) # Mass conservation
    dv[2] = v[2] - (enthalpy_in / massflow_out) # node_temp = enthalpy_in / massflow_out

    return nothing
end

function fixed_node!(dv, v, _, _, p, _)
    # DirectedODEVertex, dims == 2
    # dv[1:2] = 0.0

    # Physics implementation
    dv[1] = v[1] - p.p_ref
    dv[2] = v[2] - p.T_fixed # TODO: Remove T_fixed, calculate nodal temperature like junction nodes
    return nothing
end

# function prosumer_massflow(index::Integer)
#     # Returns closure with 'index' captured
#     function f!(de, e, v_s, _, p, _)
#         # Closure, implements physics for prosumer edges with fixed massflow
#         let
#             index = index

#             # Prosumer edges must always have dims == 3
#             # de[1:3] == 0, algebraic constraint

#             # Get local parameters from Parameters struct
#             massflow = p.edge_parameters.massflow[index]
#             deltaT = p.edge_parameters.deltaT[index]

#             # Physics implementation
#             de[1] = massflow - e[1] # Fixed mass flow rate
#             de[2] = e[2] - v_s[2] # Upwind convection, edge start temp. == source node temp.
#             de[3] = deltaT - (e[3] - e[2]) # Fixed temperature change across edge

#             return nothing
#         end # let block
#     end # f!(...)
#     return f!
# end


# function prosumer_deltaP(index::Integer)
#     # Returns closure with 'index' captured
#     function f!(de, e, v_s, v_d, p, _)
#         # Closure, implements physics for prosumer edges with fixed pressure change
#         let
#             index = index

#             # Prosumer edges must always have dims == 3
#             # de[1:3] == 0, algebraic constraint

#             # Get local parameters from Parameters struct
#             deltaP = p.edge_parameters.deltaP[index]
#             deltaT = p.edge_parameters.deltaT[index]

#             # Physics implementation
#             de[1] = deltaP - (v_d[1] - v_s[1]) # Fixed pressure difference
#             de[2] = e[2] - v_s[2] # Upwind convection, edge start temp. == source node temp.
#             de[3] = deltaT - (e[3] - e[2]) # Fixed temperature change across edge

#             return nothing
#         end # let block
#     end # f!(...)
#     return f!
# end


# function prosumer(pressure_control::Function,
#         heatrate_control::Function,
#         hydraulic_characteristic::Function,
#         coeff_fns::TransportProperties)
#     # Returns closure
#     function f!(de, e, v_s, v_d, _, t)
#         # Closure, implements physics for prosumer edges
#         let
#             deltaP = pressure_control
#             heatrate = heatrate_control
#             characteristic = hydraulic_characteristic
#             specific_heat::Float64 = coeff_fns.heat_capacity # Assuming spec. heat is constant

#             # Calculate local variables
#             massflow::Float64 = deltaP(t) |> characteristic # Assuming massflow is always in the direction of decreasing pressure
#             inlet_T = massflow >= 0 ? v_s[2] : v_d[2] # Upwind convection
#             outlet_T = heatrate(t) / (massflow * specific_heat) + inlet_T

#             # Physics implementation
#             de[1] = massflow - e[1] # state 1 : massflow
#             de[2] = inlet_T - e[2] # state 2 : inlet temperature
#             de[3] = outlet_T - e[3] # state 3 : outlet temperature

#             return nothing
#         end # let block
#     end # f!(...)

#     return f!
# end

end # module
