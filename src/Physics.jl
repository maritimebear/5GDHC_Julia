# Dynamical functions for NetworkDynamics ODEEdge and DirectedODEVertex

module Physics

include("./InterpolationSchemes.jl") # provides module FVM
import .FVM

export pipe!, prosumer!, junction!, reference_node!

function pipe!(de, e, v_s, v_d, p, _)
    # Calculate local variables
    area = 0.25 * pi * (p.diameter ^ 2)
    velocity = e[1] / (p.density * area)
    Re = p.density * velocity * p.diameter / p.dyn_visc

    # e[1] : mass flow rate, algebraic constraint
    #   => de[1] == 0, used to calculate pressure drop across pipe due to friction
    # TODO: Implement pressure loss according to Cengel eqn. 8-21,
    #       Churchill or Swamee-Jain approximation for Darcy-Weisbach f
    de[1] = (v_s[1] - v_d[1]) + (p.friction_fn(Re) * velocity * abs(velocity)) # Pressure drop due to friction

    # Energy equation, e[2:end] => temperatures in finite-volume cells
    @views convection = -(1 / p.dx) .* FVM.upwind(e[2:end], v_s[2], v_d[2], velocity)
    @views source = p.htrans_coeff() .* (e[2:end] .- p.T_ambient)
    @views de[2:end] .= convection .+ source
    return nothing
end

function prosumer!(de, e, v_s, v_d, p, _)
    # Prosumer edges must always have dims == 3
    # de[1:3] == 0, algebraic constraint

    de[1] = p.massflow - e[1] # Fixed mass flow rate
    de[2] = e[2] - v_s[2] # Upwind convection, edge start temp. == source node temp.
    de[3] = p.delta_T - (e[3] - e[2]) # Fixed temperature change across edge

    return nothing
end

function junction!(dv, v, edges_in, edges_out, _, _)
    # DirectedODEVertex, dims == 2
    # dv[1:2] = 0.0
    dv[1] = sum(map(e -> e[1], edges_in)) - sum(map(e -> e[1], edges_out)) # Mass conservation

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

    dv[2] = v[2] - (enthalpy_in / massflow_out) # node_temp = enthalpy_in / massflow_out

    return nothing
end

function reference_node!(dv, v, _, _, p, _)
    # DirectedODEVertex, dims == 2
    # dv[1:2] = 0.0
    dv[1] = v[1] - p.p_ref
    dv[2] = v[2] - p.T_fixed
    return nothing
end

end # module
