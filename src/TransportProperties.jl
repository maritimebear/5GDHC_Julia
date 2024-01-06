module TransportProperties # submodule, included in DHG.jl
# Dynamical functions for NetworkDynamics ODEEdge and DirectedODEVertex

export TransportProperties


Base.@kwdef struct TransportCoefficients{Functor1, Functor2, Functor3, Functor4}
    dynamic_viscosity::Functor1
    wall_friction::Functor2
    heat_capacity::Functor3
    heat_transfer::Functor4
end


function friction_Churchill(Re::Real) # -> Float64
    # TODO
    return # float
end


end # (sub)module
