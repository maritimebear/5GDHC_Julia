module ParameterStructs


export EdgeParameters, NodeParameters, GlobalParameters, Parameters


Base.@kwdef struct EdgeParameters{SparseVectorType}
    # Immutable struct, attribute arrays can still be modified:
    #   eg. edge_parameters.diameter[i] = new_value

    # Pipe edge parameters
    diameter::SparseVectorType
    length::SparseVectorType
    dx::SparseVectorType
    # Prosumer edge parameters
    massflow::SparseVectorType
    delta_p::SparseVectorType
    delta_T::SparseVectorType
end


Base.@kwdef struct NodeParameters
    p_ref::Float64
    T_fixed::Float64
end


Base.@kwdef struct GlobalParameters{Functor1, Functor2, Functor3, Functor4}
    # TODO: dynamic viscosity, heat capacity as constants or functions?
    density::Float64
    T_ambient::Float64
    dynamic_viscosity::Functor1
    friction_coeff::Functor2
    heat_capacity::Functor3
    heat_loss_coeff::Functor4
end


Base.@kwdef struct Parameters
    # Structure of structures of arrays, for modification via DifferentialEquations callbacks:
    #   integrator.p.edge_parameters.whatever[index] = new_value
    global_parameters::GlobalParameters
    node_parameters::NodeParameters
    edge_parameters::EdgeParameters
end


end # module
