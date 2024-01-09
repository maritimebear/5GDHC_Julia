module ParameterStructs # submodule, included in DHG.jl

export EdgeParameters, NodeParameters, GlobalParameters, Parameters


import ..NetworkComponents as nc
import ..GraphParsing as gp

import SparseArrays as sp


# Structs to hold parameters, modifiable via DifferentialEquations callbacks

Base.@kwdef struct EdgeParameters{ValueType <: Real, IndexType <: Integer}
    # Struct to hold modifiable parameters for prosumer edges.
    # Immutable struct, constituent arrays can still be modified:
    #   eg. edge_parameters.massflow[i] = new_value
    # No modifiable parameters for pipe edges.

    # Prosumer edge parameters
    massflow::sp.SparseVector{ValueType, IndexType}
    deltaP::sp.SparseVector{ValueType, IndexType}
    deltaT::sp.SparseVector{ValueType, IndexType}
end


Base.@kwdef mutable struct NodeParameters{ValueType <: Real}
    # mutable struct: p_ref, T_fixed can be changed
    p_ref::ValueType
    T_fixed::ValueType
end


Base.@kwdef mutable struct GlobalParameters{ValueType <: Real}
    # mutable struct: T_ambient can vary over time
    density::ValueType
    T_ambient::ValueType
end


Base.@kwdef struct Parameters{ValueType <: Float64, IndexType <: Integer}
    # Structure of structures of arrays, for modification via DifferentialEquations callbacks:
    #   integrator.p.edge_parameters.whatever[index] = new_value
    # ValueType must be Float64

    global_parameters::GlobalParameters{ValueType}
    node_parameters::NodeParameters{ValueType}
    edge_parameters::EdgeParameters{ValueType, IndexType}
end


# Custom constructors

function NodeParameters(node_dict::gp.ComponentDict{IdxType, nc.Node}) where {IdxType <: Integer}
    # parse_gml() -> ComponentDict |> this constructor -> NodeParameters struct
    if (length(node_dict.indices[:fixed]) != 1) # Should be false, checked already in parse_gml()
        throw(ArgumentError("multiple fixed nodes defined"))
    end
    fixed_idx = node_dict.indices[:fixed][1]::IdxType
    return NodeParameters(p_ref=node_dict.components[fixed_idx].pressure,
                          T_fixed=node_dict.components[fixed_idx].temperature
                         )
end


function EdgeParameters(edge_dict::gp.ComponentDict{IdxType, nc.Edge}, ::Type{ValueType}=Float64) where {IdxType <: Integer, ValueType <: Real}
    # -> EdgeParameters{ValueType, IdxType}
    # parse_gml() -> ComponentDict |> this constructor -> EdgeParameters struct
    # !! edge_dict must be sorted in the same order as the edges in Graphs.jl object !!

    # Since all edge types have parameters (unlike nodes, where only the fixed node has parameters),
    # storing edge parameters in sparse vectors, structure of arrays (SoA) layout.

    sparsevec_length = Base.length(edge_dict.components) # SoA => all sparse vectors have the same length; TODO: UndefVarError unless Base.length()?
    # Vector of indices, for construction of sparse vectors
    prosumer_idxs = sort([edge_dict.indices[:massflow]; edge_dict.indices[:deltaP]]) # Concatenate and sort

    # Convenience functions
    @inline function get_edge_values(field::Symbol, idxs::Vector{IdxType}) # -> Vector{ValueType}
        # Handle no deltaP, no massflow or no pipe edges while maintaining consistent SparseVector{IdxType, ValueType} to EdgeParameters inner ctor
        result = [getfield(edge_dict.components[i], field) for i in idxs]
        return (isempty(result) ? zeros(ValueType, 0) : result)::Vector{ValueType}
    end

    @inline function fwd_to_ctor(idxs::Vector{IdxType}, value_symbol::Symbol) # -> Vector{SparseVector{ValueType, IdxType}
        return sp.SparseVector(sparsevec_length, idxs,
                               get_edge_values(value_symbol, idxs))::sp.SparseVector{ValueType, IdxType}
    end

    # Sparse vectors
    diameter = fwd_to_ctor(edge_dict.indices[:pipe], :diameter)
    length = fwd_to_ctor(edge_dict.indices[:pipe], :length)
    dx = fwd_to_ctor(edge_dict.indices[:pipe], :dx)
    massflow = fwd_to_ctor(edge_dict.indices[:massflow], :massflow)
    deltaP = fwd_to_ctor(edge_dict.indices[:deltaP], :deltaP)
    deltaT = fwd_to_ctor(prosumer_idxs, :deltaT)

    return EdgeParameters(diameter=diameter, length=length, dx=dx,
                          massflow=massflow, deltaP=deltaP, deltaT=deltaT)
end


end # (sub)module
