module ParameterStructs # submodule, included in DHG.jl

export EdgeParameters, NodeParameters, GlobalParameters, Parameters


import ..NetworkComponents as nc
import ..GraphParsing as gp

import SparseArrays as sp


# Structs to hold parameters, modifiable via DifferentialEquations callbacks

Base.@kwdef struct EdgeParameters{SparseVectorType}
    # Immutable struct, attribute arrays can still be modified:
    #   eg. edge_parameters.diameter[i] = new_value

    # Pipe edge parameters
    diameter::SparseVectorType
    length::SparseVectorType
    dx::SparseVectorType
    # Prosumer edge parameters
    massflow::SparseVectorType
    deltaP::SparseVectorType
    deltaT::SparseVectorType
end


Base.@kwdef mutable struct NodeParameters{T<:Real}
    # mutable struct <= p_ref, T_fixed can be changed
    p_ref::T
    T_fixed::T
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


# Custom constructors

function NodeParameters(node_dict::gp.ComponentDict{IdxType, nc.Node}) where {IdxType <: Integer}
    # parse_gml() -> ComponentDict |> this constructor -> NodeParameters struct
    if (length(node_dict.indices[:fixed]) != 1) # Should be false, checked already in parse_gml()
        throw(ArgumentError("multiple fixed nodes defined"))
    end
    fixed_idx::IdxType = node_dict.indices[:fixed][1]
    return NodeParameters(p_ref=node_dict.components[fixed_idx].pressure,
                          T_fixed=node_dict.components[fixed_idx].temperature
                         )
end


end # (sub)module
