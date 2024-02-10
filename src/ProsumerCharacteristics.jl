module ProsumerCharacteristics
# Functions to implement prosumer dynamics by modelling pumps, boilers/condensers

@inline function temperature_out(thermal_power::Real, massflow::Real, temperature_in::Real, spec_heat::Real) 
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


end # module
