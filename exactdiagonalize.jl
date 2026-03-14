using MKL, LinearAlgebra
using SparseArrays

include("utils.jl")

abstract type AbstractOpSum end

mutable struct SpinOpSum <: AbstractOpSum
    type::DataType
    opvec::AbstractVector
end

function act(opstr::String, loc::Int, bits::Int, T::DataType)
    if opstr == "Z"
        return bits, T(2 * readbit(bits, loc) - 1)
    elseif opstr == "X"
        return flip(bits, loc), one(T)
    elseif opstr == "σ+"
        return flip(bits, loc), T(! readbit(bits, loc))
    elseif opstr == "σ-"
        return flip(bits, loc), T(readbit(bits, loc))
    # We do not specify the Y operator to keep type stability
    else
        error("Operator not specified yet!")
    end
end

function act(opstr::String, loc::Tuple{Int, Int}, bits::Int, T::DataType)
    if opstr == "CZ"

    end
end

function trans(op::Tuple, bits::Int)
    element = op[1]
    newbits = bits
    oplen = length(op)

    for s in 2:2:oplen
        tmp = act(op[s], op[s+1], newbits, typeof(element))
        newbits = tmp[1]
        element *= tmp[2]
    end
    return newbits, element
end

function apply!(hmat::Matrix{T}, ops::SpinOpSum, basis::Vector{<:Int}, idx::Int, num::Int) where T <: Number
    bits = basis[idx]
    for op in ops.opvec
        newbits, element = trans(op, bits)
        if count_ones(newbits) == num && ! iszero(element)
            newidx = searchsortedfirst(basis, newbits)
            hmat[newidx, idx] += element
        end
    end
end


function makeHamiltonian(ops::AbstractOpSum, basis::Vector{<:Int}, num::Int)
    dim = length(basis)
    hmat = zeros(ops.type, dim, dim)
    for j in 1:dim
        apply!(hmat, ops, basis, j, num)
    end
    return hmat
end
