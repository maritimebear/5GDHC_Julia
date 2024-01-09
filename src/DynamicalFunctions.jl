module DynamicalFunctions # submodule, included in DHG.jl
# Dynamical functions for NetworkDynamics ODEEdge and DirectedODEVertex

export pipe, prosumer_massflow, prosumer_deltaP, junction!, fixed_node!


import ..FVM
import ..TransportProperties: TransportCoefficients


# TODO: 'let' variables in closures for performance
# TODO: Type annotations for local parameters? (eg. density::T where T is captured from outer function?)
#           # Probably better to annotate p::ParameterStructs.Parameters{ValueType, IndexType} in function signature
# TODO: Refactor junction!() with pipes to chain filter |> map |> sum


# Dynamical functions for edges

## Each edge can have different parameters, so each edge function is a closure with an 'index' captured.
## 'index' is used to access corresponding values from ParameterStructs.EdgeParameters.

function pipe(diameter::Float64, dx::Float64, coeff_fns::TransportCoefficients{F1, F2, F3, F4}) where {F1, F2, F3, F4}
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
            density = p.global_parameters.density
            T_ambient = p.global_parameters.T_ambient
            # Calculate local variables
            velocity = e[1] / (density * area)
            Re = density * velocity * diameter / dyn_visc()

            # Physics implementation
            # e[1] : mass flow rate, algebraic constraint
            #   => de[1] == 0, used to calculate pressure drop across pipe due to friction
            # TODO: Implement pressure loss according to Cengel eqn. 8-21,
            #       Churchill or Swamee-Jain approximation for Darcy-Weisbach f
            de[1] = (v_s[1] - v_d[1]) + (friction(Re) * velocity * abs(velocity)) # Pressure drop due to friction

            # Energy equation, e[2:end] => temperatures in finite-volume cells
            @views convection = -(1 / dx) .* FVM.upwind(e[2:end], v_s[2], v_d[2], velocity)
            @views source = htrans_coeff() .* (e[2:end] .- T_ambient)
            @views de[2:end] .= convection .+ source
            return nothing
        end # let block
    end # f!(...)

    return f!
end


function prosumer_massflow(index::Integer)
    # Returns closure with 'index' captured
    function f!(de, e, v_s, _, p, _)
        # Closure, implements physics for prosumer edges with fixed massflow
        let
            index = index

            # Prosumer edges must always have dims == 3
            # de[1:3] == 0, algebraic constraint

            # Get local parameters from Parameters struct
            massflow = p.edge_parameters.massflow[index]
            deltaT = p.edge_parameters.deltaT[index]

            # Physics implementation
            de[1] = massflow - e[1] # Fixed mass flow rate
            de[2] = e[2] - v_s[2] # Upwind convection, edge start temp. == source node temp.
            de[3] = deltaT - (e[3] - e[2]) # Fixed temperature change across edge

            return nothing
        end # let block
    end # f!(...)
    return f!
end


function prosumer_deltaP(index::Integer)
    # Returns closure with 'index' captured
    function f!(de, e, v_s, v_d, p, _)
        # Closure, implements physics for prosumer edges with fixed pressure change
        let
            index = index

            # Prosumer edges must always have dims == 3
            # de[1:3] == 0, algebraic constraint

            # Get local parameters from Parameters struct
            deltaP = p.edge_parameters.deltaP[index]
            deltaT = p.edge_parameters.deltaT[index]

            # Physics implementation
            de[1] = deltaP - (v_d[1] - v_s[1]) # Fixed pressure difference
            de[2] = e[2] - v_s[2] # Upwind convection, edge start temp. == source node temp.
            de[3] = deltaT - (e[3] - e[2]) # Fixed temperature change across edge

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
    dv[1] = v[1] - p.node_parameters.p_ref
    dv[2] = v[2] - p.node_parameters.T_fixed
    return nothing
end

end # (sub)module
