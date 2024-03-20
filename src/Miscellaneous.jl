module Miscellaneous

import NetworkDynamics as nd
import DifferentialEquations as de

export solve_steadystate, solve_dynamic
export PumpModel
export set_idxs, initialise
export solve_steadystate, solve_dynamic
export adjacent_find


function PumpModel(massflow_ref1, deltaP_ref1, speed_ref1,
                    massflow_ref2, deltaP_ref2, speed_ref2,
                    fluid_density, nominal_speed
                )
    # Returns closure with reference values captured
    # Following Licklederer et al, "Thermohydraulic model of Smart Thermal Grids with bidirectional power flow between prosumers", 2021, pages 4, 12
    
    # Unit conversions to conform to Licklederer et al
    # non-dimensionalise [rpm]:
    u_1 = speed_ref1 / nominal_speed
    u_2 = speed_ref2 / nominal_speed 
    # [kg/s] -> [l/min]:
    V_1 = massflow_ref1 / fluid_density * 60 * 1e3
    V_2 = massflow_ref2 / fluid_density * 60 * 1e3
    # [Pa] -> [hPa]:
    dP_1 = deltaP_ref1 * 1e-2
    dP_2 = deltaP_ref2 * 1e-2

    c1 = (dP_1 - (u_1 / u_2)^2 * dP_2) / (V_1^2 - (u_1 / u_2)^2 * V_2^2)
    c2 = (1.0 / u_2)^2 * (dP_2 - (V_2^2 * c1))

    function model(pump_speed, massflow)
        # Returns pressure increase across pump
        # Constant outer values in let-block for performance
        let c1 = c1
            c2 = c2
            nominal_speed = nominal_speed
            fluid_density = fluid_density
            # Local variables
            speed_ratio = pump_speed / nominal_speed
            vol_rate = massflow / fluid_density * 60 * 1e3 # [kg/s] -> [l/min]
            return ((c1 * vol_rate^2) + (c2 * speed_ratio^2)) * 1e2 # [hPa] -> [Pa]
        end
    end

    return model
end


function set_idxs(state_vector::Vector{Float64}, key_value::Pair{K, V}, nd_fn) where {
                    K <: Union{Symbol, AbstractString, Regex}, V <: Union{Real, Function}
                   }
    # -> Nothing
    # nd_fn: return from network_dynamics() , <: SciMLBase.ODEFunction
    # Wrapper around NetworkDynamics::idx_containing(),
    #   sets values of state_vector at indices matching keys
    # If V <: Function, then V(n) must return Real or Vector{Real},
    #   n::Integer = number of indices matching keys

    idxs = nd.idx_containing(nd_fn, key_value.first)::Vector{<:Integer}

    if V <: Real
        state_vector[idxs] .= key_value.second
            # all elements at state_vector[idxs] get the same value,
            # implicit conversion: V -> Float64
        return nothing

    else # V <: Function
        val = key_value.second(length(idxs))::Union{Real, Vector{<:Real}}
        state_vector[idxs] .= val
            # val::Real => all elements get the same value,
            # val::Vector => val gets copied into state_vector[idxs], implicit dimension check
        return nothing
    end

end


function initialise(nd_fn,
                    massflow_init::T1,
                    pressure_init::T2,
                    temperature_init::T3
                ) where {T1, T2, T3 <: Union{Real, Function}}
    # -> Vector{Float64}
    # Calls/assigns each _init to the appropriate states in initial state vector
    # nd_fn: return from network_dynamics() , <: SciMLBase.ODEFunction
    #
    # Symbols for massflow, pressure and temperature states in nd_fn set in WrapperFunctions.jl,
    # assuming these to be unchanged

    n_states = length(nd_fn.syms)
    initial_state = Vector{Float64}(undef, n_states)

    for (k, v) in [(:m => massflow_init),
                   (:p => pressure_init),
                   (:T => temperature_init),
                  ]
        set_idxs(initial_state, (k => v), nd_fn) # set_idxs() defined in this file
    end

    return initial_state
end


function check_retcode(solution, solvername)
    # Checks return code from de.solve() result, raise warning if not successful
    if solution.retcode !== de.ReturnCode.Success
        @warn "Unsuccessful retcode from $(solvername): $(solution.retcode)"
    end
end


function solve_steadystate(f, x0, p, solver; kwargs...)
    # -> de.solve(de.SteadyStateProblem, solver) type
    # Wrapper to make main script more readable
    prob = de.SteadyStateProblem(f, x0, p)
    sol = de.solve(prob, solver; kwargs...)
    check_retcode(sol, "steady-state solver")
    return sol
end


function solve_dynamic(f, x0, p, solver, time_interval; kwargs...)
    # -> de.solve(de.ODEProblem, solver, saveat=save_times) type
    # Wrapper to make main script more readable
    prob = de.ODEProblem(f, x0, time_interval, p)
    sol = de.solve(prob, solver; kwargs... #saveat=save_times
                  )
    check_retcode(sol, "dynamic solver")
    return sol
end


function adjacent_find(binary_predicate, array)
    # -> Int
    # Compare adjacent elements of 'array' using 'binary_predicate'.
    # Returns array index of first 'true' from 'binary_predicate', or last index if no match found
    # Following C++ std::adjacent_find(), https://en.cppreference.com/w/cpp/algorithm/adjacent_find
    # No find-like function in Julia that takes a binary predicate?

    if length(array) == 0
        throw(DomainError("invalid argument: empty array"))
    end

    for i in 1:(length(array) - 1)
        @inbounds if binary_predicate(array[i], array[i+1])
            return i
        end
    end
    return length(array) # No match found
end


end # module
