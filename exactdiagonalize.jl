include("operators.jl")

function timeEvolve(ops::AbstractOpSum, init::AbstractState, tf::Real)
    hmat = makeHamiltonian(ops, init.basis)
    eigenergy, U = eigen(hmat)
    phases = cos.(tf * eigenergy) .- im * sin.(tf * energy)
    expEt = Diagonal(phases)
    final = U * expEt * U' * (init.vector)
    return State(init.basis, final)
end

function timeEvolve(ops::AbstractOpSum, init::AbstractState, ts::Vector{<:Real}, watcher::Tuple)
    hmat = makeHamiltonian(ops, init.basis)
    eigenergy, U = eigen(hmat)
    expEt = Diagonal{ComplexF64}(similar(eigenergy))
    psi = Vector{ComplexF64}(similar(init.vector))

    for t in ts
        phases .= cos.(t * eigenergy) .- im * sin.(t * eigenergy)
        expEt[diagind(expEt)] .= phases
        psi .= U * expEt * U' * (init.vector)
    end
end

let 
    L, N = 6, 3

    os = Tuple[]
    for j in 1:L
        nj = mod1(j+1, L)
        push!(os, (1.0, :Z, j, :Z, nj))
        push!(os, (1.0, :X, j, :X, nj))
        push!(os, (-1.0, :iY, j, :iY, nj))
    end
    # push!(os, (1.0, :X, L))
    ops = SpinOpSum(Float64, os)

    op = (1.0, :Z, 2)
    
    basis = NumBasis(L, N)
    vs = randn(length(basis.bitsvec))

    psi = State(basis, vs)
    normalize!(psi)
    norm(psi.vector)
end