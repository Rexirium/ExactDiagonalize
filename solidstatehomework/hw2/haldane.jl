using MKL, LinearAlgebra

const Am = [ 1.0 0.0 ;
    -1/2 √3/2 ;
    -1/2 -√3/2]
const Bm = [0.0 √3 ;
    -3/2 -√3/2 ;
    3/2 -√3/2]

function makeHaldaneHamiltonian(kv::Vector{Float64}, t1::Real, t2::Number)
    hmat = Matrix{ComplexF64}(undef, 2, 2)
    hmat[1, 1] = 2 * real(t2) * sum(cos.(Bm * kv))
    hmat[2, 2] = - 2 * imag(t2) * sum(sin.(Bm * kv))
    hmat[2, 1] = t1 * sum(cis.( Am * kv))
    hmat[1, 2] = conj(hmat[2, 1])
    return hmat
end

function updateHaldaneHamiltonian!(hmat::Matrix{ComplexF64}, kv::Vector{Float64}, t1::Real, t2::Number)
    hmat[1, 1] = 2 * real(t2) * sum(cos.(Bm * kv))
    hmat[2, 2] = - 2 * imag(t2) * sum(sin.(Bm * kv))
    hmat[2, 1] = t1 * sum(cis.( Am * kv))
    hmat[1, 2] = conj(hmat[2, 1])
end