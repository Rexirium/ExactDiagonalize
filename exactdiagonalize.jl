
function spectrum(ops::OpSum, basis::AbstractBasis)
    hmat = makeHamiltonian(ops, basis)
    return eigvals(Hermitian(hmat))
end

function spectrum_numconserved(ops::OpSum, lsize::Int)
    energies = Float64[]
    sizehint!(energies, 1 << lsize)
    for num in 0:lsize
        basis = NumBasis(lsize, num)
        hmat = makeHamiltonian(ops, basis)
        eigs = eigvals!(Hermitian(hmat))
        append!(energies, eigs)
    end
    return energies
end

function timeEvolve_exact(ops::OpSum, init::AbstractState, tf::Real)
    hmat = makeHamiltonian(ops, init.basis)
    eigenergy, U = eigen(Hermitian(hmat))
    phases = complex.(cos.(tf * eigenergy), - sin.(tf * eigenergy))
    expEt = Diagonal(phases)
    final = U * expEt * U' * (init.vector)
    return State(init.basis, final)
end

function timeEvolve_exact(ops::OpSum, init::AbstractState, ts::AbstractVector, obs::AbstractObserver)
    hmat = makeHamiltonian(ops, init.basis)
    eigenergy, U = eigen(Hermitian(hmat))
    dim = length(eigenergy)

    phases = Vector{ComplexF64}(undef, dim)
    expEt = Diagonal{ComplexF64}(undef, dim)
    psi = ComplexF64.(init.vector)

    record!(obs, psi)
    for t in ts[2:end]
        phases .= complex.(cos.(t * eigenergy), - sin.(t * eigenergy))
        expEt.diag .= phases
        mul!(psi, U * expEt * U', init.vector)
        record!(obs, psi)
    end
    return State(init.basis, psi)
end