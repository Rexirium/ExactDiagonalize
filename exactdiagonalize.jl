using MKL, LinearAlgebra
using SparseArrays

include("utils.jl")

abstract type AbstractOpSum end
abstract type AbstractBasis end

mutable struct SpinOpSum <: AbstractOpSum
    type::DataType
    opvec::AbstractVector
end

mutable struct NumBasis <: AbstractBasis
    num::Int
    bitsvec::Vector{<:Int}

    NumBasis(lsize::Int, num::Int) = new(num, numbitbasis(lsize, num))
end

function act(opstr::String, loc::Int, bits::Int, T::DataType)
    """
    act a single qubit operator on the state `bits`=|1001011⟩ for bits=(1001011)₂
    |1⟩ = (1, 0)ᵀ = |↑⟩, |0⟩ = (0, 1)ᵀ = |↓⟩
    I do not specify the Y operator (has complex element) to keep type stability.
    """
    if opstr == "Z"
        return bits, T(2 * readbit(bits, loc) - 1)
    elseif opstr == "X"
        return flip(bits, loc), one(T)
    elseif opstr == "Sp" # means simplectic matrix [0 1 ; -1 0], Sp = iY
        return flip(bits, loc), T(1 - 2 * readbit(bits, loc))
    elseif opstr == "σ+"
        return flip(bits, loc), T(! readbit(bits, loc))
    elseif opstr == "σ-"
        return flip(bits, loc), T(readbit(bits, loc))
    else
        error("Operator not specified yet!")
    end
end

function act(opstr::String, loc::Tuple{Int, Int}, bits::Int, T::DataType)
    if opstr == "CX"
        c, t = loc
        bitc = readbit(bits, c)
        return flip(bits, t, bitc), one(T)
    elseif opstr == "CZ"
        c, t = loc
        i1, i2 = minmax(c, t)
        b1, b2 = readbit(bits, i1, i2)
        return bits, T(2 * (b1 ^ b2) - 1)
    else
        error("Operator not specified yet!")
    end
end

function apply(op::Tuple, bits::Int, T::DataType)
    element = op[1]
    newbits = bits
    oplen = length(op)

    for s in 2:2:oplen
        tmp = act(op[s], op[s+1], newbits, T)
        newbits = tmp[1]
        element *= tmp[2]
    end
    return newbits, element
end

function makeHamiltonian(ops::AbstractOpSum, basis::NumBasis)
    dim = length(basis.bitsvec)
    hmat = zeros(ops.type, dim, dim)
    for (j, bits) in enumerate(basis.bitsvec)
        for op in ops.opvec
            newbits, element = apply(op, bits, ops.type)
            if count_ones(newbits) == basis.num && !iszero(element)
                i = searchsortedfirst(basis.bitsvec, newbits)
                hmat[i, j] += element
            end
        end
    end
    return Hermitian(hmat)
end

let 
    L, N = 10, 5

    os = Tuple[]
    for j in 1:L
        nj = mod1(j+1, L)
        push!(os, (1.0, "Z", j, "Z", nj))
        push!(os, (1.0, "X", j, "X", nj))
        push!(os, (-1.0, "Sp", j, "Sp", nj))
    end
    # push!(os, (1.0, "X", L))
    ops = SpinOpSum(Float64, os)

    basis = NumBasis(L, N)
    @time Hnum = makeHamiltonian(ops, basis)
    nothing
end
