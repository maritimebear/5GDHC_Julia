module Transport # submodule, included in DHG.jl

export Transport


Base.@kwdef struct TransportProperties{T1, T2, T3, T4}
    # Members may be constants or callable types
    dynamic_viscosity::T1
    wall_friction::T2
    heat_capacity::T3
    heat_transfer::T4
end


function friction_Churchill(Re::Real) # -> Float64
    # TODO
    return # float
end


end # (sub)module
