using MKL, LinearAlgebra

include("utils.jl")

abstract type AbstractOpSum end

mutable struct SpinOpSum <: AbstractOpSum
    type::DataType
    opvec::AbstractVector
end

function apply!(hmat::Matrix{T}, ops::SpinOpSum, num::Int) where T <: Number
    for op in ops.opvec
        len = length(op)
    end
end


function makeHamiltonian(ops::AbstractOpSum, basis::Vector{<:Int})
    dim = length(basis)
    hmat = Matrix{ops.type}(undef, dim, dim)
    for j in 1:dim
        apply!(hmat, ops, basis[j])
    end
    return hmat
end