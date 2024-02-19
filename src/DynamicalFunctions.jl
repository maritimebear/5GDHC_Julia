module DynamicalFunctions

# Dynamical functions for NetworkDynamics ODEEdge and DirectedODEVertex

import ..Transport: TransportProperties
import ..NetworkComponents as nc
import ..FVM
import ..Fluids as fl

export prosumer_outlet_T, prosumer_deltaP, prosumer_massflow, pipe, node_temperature, junction!, reference_node


@inline function prosumer_outlet_T(thermal_power::Real, massflow::Real, temperature_in::Real, spec_heat::Real) 
    # -> T_out = Q / (m * Cp) + T_in
    # Requires Q, m, Cp and T_in to be known
    # m always positive, heat transfer direction specified by sign of Q
    return thermal_power / (abs(massflow) * spec_heat) + temperature_in
end


function prosumer_deltaP(prosumerstruct::nc.PressureChange, ::Type{fluid_T}) where {fluid_T <: fl.Fluid}
    # Returns closure with constants captured
    function f!(de, e, v_s, v_d, _, t)
        # Closure, implements physics for pressure change prosumer
        # Prosumer edges must always have dims == 3
        # de[1:3] == 0, algebraic constraints
        let
            # Capture constants in let-block for performance
            # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured
            hyd_ctrl = prosumerstruct.hydraulic_control
            thm_ctrl = prosumerstruct.thermal_control
            hyd_chr = prosumerstruct.hydraulic_characteristic

            # Local variables
            T_mean = 0.5 * (e[2] + e[3]) # TODO: switch to LMTD
            spec_heat = fl.specific_heat(fluid_T, T_mean)
            deltaP = hyd_chr(hyd_ctrl(t), e[1]) # Calculate pressure change from dynamic control input and massflow
            m_aligned = e[1] >= 0 # true <=> massflow is along edge direction or zero
            ## inlet_T, outlet_T defined wrt massflow direction: inlet_T always upstream of massflow
            inlet_T = m_aligned ? v_s[2] : v_d[2] # Upwind convection
            outlet_T = prosumer_outlet_T(thm_ctrl(t), e[1], inlet_T, spec_heat)

            ## T_src, T_dst wrt edge directivity, temperatures at edge-node interfaces
            T_src = m_aligned ? inlet_T : outlet_T
            T_dst = m_aligned ? outlet_T : inlet_T

            # Physics implementation
            # e[1] : mass flow rate, algebraic constraint
            # e[2] : temperature at interface with source node, algebraic constraint
            # e[3] : temperature at interface with destination node, algebraic constraint
            de[1] = deltaP - (v_d[1] - v_s[1])
            de[2] = T_src - e[2]
            de[3] = T_dst - e[3]

            return nothing
        end # let block
    end # f!(...)
    return f!
end


function prosumer_massflow(prosumerstruct::nc.Massflow, ::Type{fluid_T}) where {fluid_T <: fl.Fluid}
    # Returns closure with constants captured
    function f!(de, e, v_s, v_d, _, t)
        # Closure, implements physics for massflow prosumer
        let
            # Capture constants
            hyd_ctrl = prosumerstruct.hydraulic_control
            thm_ctrl = prosumerstruct.thermal_control
            hyd_chr = prosumerstruct.hydraulic_characteristic

            # Local variables
            T_mean = 0.5 * (e[2] + e[3]) # TODO: switch to LMTD
            spec_heat = fl.specific_heat(fluid_T, T_mean)
            massflow = hyd_chr(hyd_ctrl(t), e[1])
            m_aligned = e[1] >= 0 # true <=> massflow is along edge direction or zero
            ## inlet_T, outlet_T defined wrt massflow direction: inlet_T always upstream of massflow
            inlet_T = m_aligned ? v_s[2] : v_d[2] # Upwind convection
            outlet_T = prosumer_outlet_T(thm_ctrl(t), e[1], inlet_T, spec_heat)

            # TODO: verify, repeat for deltaP prosumer, correct syms, comments
            ## T_src, T_dst wrt edge directivity, temperatures at edge-node interfaces
            T_src = m_aligned ? inlet_T : outlet_T
            T_dst = m_aligned ? outlet_T : inlet_T

            # Physics implementation
            de[1] = massflow - e[1] # massflow through edge
            de[2] = T_src - e[2]    # source node interface temperature
            de[3] = T_dst - e[3]    # destination node interface temperature

            return nothing
        end # let block
    end # f!(...)
    return f!
end


function pipe(pipestruct::nc.Pipe, transport_coeffs::TransportProperties, ::Type{fluid_T}) where {fluid_T <: fl.Fluid}
    # Returns closure with constants captured
    function f!(de, e, v_s, v_d, p, _)
        # Closure, implements physics for pipe edges: wall friction, heat loss to environment
        let
            diameter = pipestruct.diameter
            dx = pipestruct.dx
            friction = transport_coeffs.wall_friction
            htrans_coeff = transport_coeffs.heat_transfer

            # Constants TODO: Mark const?
            area = 0.25 * pi * (pipestruct.diameter ^ 2) # cross-sectional area, velocity calculation
            area_curved = pi * pipestruct.diameter * pipestruct.dx # heat transfer area
            rel_roughness = pipestruct.roughness / pipestruct.diameter
            aspect_ratio = pipestruct.length / pipestruct.diameter

            # Get local parameters from Parameters struct
            density = p.density # density in parameters (not based on fluid_T) since constant, not a function of temperature
            T_ambient = p.T_ambient
            # Calculate local variables
            T_mean = 0.5 * (e[2] + e[end])
            dyn_visc = fl.dynamic_viscosity(fluid_T, T_mean)

            velocity = e[1] / (density * area)
            Re = density * velocity * diameter / dyn_visc
            dynamic_pressure = 0.5 * density * (velocity^2)

            # Momentum equation
            deltaP = -sign(velocity) * friction(Re, rel_roughness) * aspect_ratio * dynamic_pressure
            # Pressure drop due to friction, Cengel & Cimbala, Fluid Mechanics: Fundamentals and Applications, 4th ed., equation 8-21
            ##   -sign(velocity): massflow > 0 => deltaP < 0, massflow < 0 => deltaP > 0; deltaP = pressure(dst) - pressure(src)

            # Energy equation
            @views convection = -(1 / dx) .* FVM.upwind(e[2:end], v_s[2], v_d[2], velocity)
            @views source = (htrans_coeff * area_curved) .* (e[2:end] .- T_ambient)

            # Physics implementation
            # e[1] : mass flow rate, algebraic constraint
            #   => de[1] == 0, used to enforce pressure drop across pipe
            # e[2:end] : temperatures in finite-volume cells
            de[1] = deltaP - (v_d[1] - v_s[1])
            de[2:end] .= convection .+ source

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
            return nothing
        end
    end
    return f!
end
end # module
