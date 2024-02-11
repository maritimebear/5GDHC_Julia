module DynamicalFunctions

# Dynamical functions for NetworkDynamics ODEEdge and DirectedODEVertex

import ..Transport: TransportProperties
import ..NetworkComponents as nc
import ..FVM

export prosumer_outlet_T, prosumer_deltaP, prosumer_massflow, pipe, node_temperature, junction!, reference_node


@inline function prosumer_outlet_T(thermal_power::Real, massflow::Real, temperature_in::Real, spec_heat::Real) 
    # -> T_out = Q / (m * Cp) + T_in
    # Requires Q, m, Cp and T_in to be known
    return thermal_power / (massflow + spec_heat) + temperature_in
end


function prosumer_deltaP(prosumerstruct::nc.PressureChange, transport_coeffs::TransportProperties)
    # Returns closure with constants captured
    function f!(de, e, v_s, v_d, _, t)
        # Closure, implements physics for pressure change prosumer
        # Prosumer edges must always have dims == 2
        # de[1:2] == 0, algebraic constraints
        let
            # Capture constants in let-block for performance
            # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured
            spec_heat = transport_coeffs.heat_capacity
            hyd_ctrl = prosumerstruct.hydraulic_control
            thm_ctrl = prosumerstruct.thermal_control
            hyd_chr = prosumerstruct.hydraulic_characteristic

            # Local variables
            deltaP = hyd_chr(hyd_ctrl(t), e[1]) # Calculate pressure change from dynamic control input and massflow
            inlet_T = e[1] >= 0 ? v_s[1] : v_d[1] # Upwind convection: decide upstream direction based on sign of massflow

            # Physics implementation
            # e[1] : mass flow rate, algebraic constraint
            # e[2] : outflow temperature, algebraic constraint
            de[1] = deltaP - (v_d[1] - v_s[1])
            de[2] = prosumer_outlet_T(thm_ctrl(t), e[1], inlet_T, spec_heat) - e[2]

            return nothing
        end # let block
    end # f!(...)
    return f!
end


function prosumer_massflow(prosumerstruct::nc.Massflow, transport_coeffs::TransportProperties)
    # Returns closure with constants captured
    function f!(de, e, v_s, v_d, _, t)
        # Closure, implements physics for massflow prosumer
        let
            # Capture constants
            spec_heat = transport_coeffs.heat_capacity
            hyd_ctrl = prosumerstruct.hydraulic_control
            thm_ctrl = prosumerstruct.thermal_control
            hyd_chr = prosumerstruct.hydraulic_characteristic

            # Local variables
            massflow = hyd_chr(hyd_ctrl(t), e[1])
            inlet_T = e[1] >= 0 ? v_s[1] : v_d[1]

            # Physics implementation
            de[1] = massflow - e[1]
            de[2] = prosumer_outlet_T(thm_ctrl(t), e[1], inlet_T, spec_heat) - e[2]

            return nothing
        end # let block
    end # f!(...)
    return f!
end


function pipe(pipestruct::nc.Pipe, transport_coeffs::TransportProperties)
    # Returns closure with constants captured
    function f!(de, e, v_s, v_d, p, _)
        # Closure, implements physics for pipe edges: wall friction, heat loss to environment
        let
            diameter = pipestruct.diameter
            dx = pipestruct.dx
            dyn_visc = transport_coeffs.dynamic_viscosity # TODO: Performance: type annotate lhs?
            friction = transport_coeffs.wall_friction
            htrans_coeff = transport_coeffs.heat_transfer

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
        end
    end
end


function node_temperature(edges_in, edges_out)
    # Calculate node temperature after mixing of incoming flows
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

    return enthalpy_in / massflow_out
end


function junction!(dv, v, edges_in, edges_out, _, _)
    # DirectedODEVertex, dims == 2
    # dv[1:2] = 0.0

    # Physics implementation
    dv[1] = sum(map(e -> e[1], edges_in)) - sum(map(e -> e[1], edges_out)) # Mass conservation
    dv[2] = v[2] - node_temperature(edges_in, edges_out)

    return nothing
end


function reference_node(node_struct::nc.ReferenceNode)
    function f!(dv, v, edges_in, edges_out, _, _)
        let p_ref = node_struct.pressure
            # DirectedODEVertex, dims == 2
            # dv[1:2] = 0.0

            # Physics implementation
            dv[1] = v[1] - p_ref
            dv[2] = v[2] - node_temperature(edges_in, edges_out)
            # dv[2] = v[2] - p.T_fixed # TODO: Remove T_fixed, calculate nodal temperature like junction nodes
            return nothing
        end
    end
    return f!
end
end # module
