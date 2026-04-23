using MKL, LinearAlgebra
using HDF5

const Am = [ 1.0 0.0 ;
    -1/2 √3/2 ;
    -1/2 -√3/2]
const Bm = [0.0 √3 ;
    -3/2 -√3/2 ;
    3/2 -√3/2]

function makeHaldaneHamiltonian(kv::Vector{Float64}, t1::Real, t2::Number; m2::Real=0.0)
    hmat = Matrix{ComplexF64}(undef, 2, 2)
    hmat[1, 1] = 2 * real(t2) * sum(cos.(Bm * kv)) + m2
    hmat[2, 2] = - 2 * imag(t2) * sum(sin.(Bm * kv)) - m2
    hmat[2, 1] = t1 * sum(cis.( Am * kv))
    hmat[1, 2] = conj(hmat[2, 1])
    return hmat
end

function updateHaldaneHamiltonian!(hmat::Matrix{ComplexF64}, kv::Vector{Float64}, t1::Real, t2::Number; m2::Real=0.0)
    hmat[1, 1] = 2 * real(t2) * sum(cos.(Bm * kv)) + m2
    hmat[2, 2] = - 2 * imag(t2) * sum(sin.(Bm * kv)) - m2
    hmat[2, 1] = t1 * sum(cis.( Am * kv))
    hmat[1, 2] = conj(hmat[2, 1])
end

function eigenHaldane(kx::Matrix, ky::Matrix, t1::Real, t2::Number; m2::Real=0.0)
    hmat = Matrix{ComplexF64}(undef, 2, 2)
    lowerband = similar(kx)
    upperband = similar(kx)
    lowervecs = Array{ComplexF64}(undef, 2, size(kx)...)

    for (idx, kv) in enumerate(zip(kx, ky))
        updateHaldaneHamiltonian!(hmat, collect(kv), t1, t2; m2=m2)
        eigs, eigvs = eigen(hmat)
        lowerband[idx] = real(eigs[1])
        upperband[idx] = real(eigs[2])
        lowervecs[2idx-1 : 2idx] = eigvs[:, 1]
    end
    return lowerband, upperband, lowervecs
end

function computeBerryCurvature(lowervecs::Array{ComplexF64}, area::Real)
    nx, ny = size(lowervecs, 2) - 1, size(lowervecs, 3) - 1
    curvature = Matrix{Float64}(undef, nx, ny)
    U = Matrix{ComplexF64}(undef, 2, 2)
    for j in 1:ny
        for i in 1:nx
            ψ1 = lowervecs[:, i, j]
            ψ2 = lowervecs[:, i+1, j]
            ψ3 = lowervecs[:, i+1, j+1]
            ψ4 = lowervecs[:, i, j+1]

            U[1, 1] = dot(ψ1, ψ2) / abs(dot(ψ1, ψ2))
            U[2, 1] = dot(ψ2, ψ3) / abs(dot(ψ2, ψ3))
            U[2, 2] = dot(ψ3, ψ4) / abs(dot(ψ3, ψ4))
            U[1, 2] = dot(ψ4, ψ1) / abs(dot(ψ4, ψ1))

            curvature[i, j] = angle(prod(U)) / area
        end
    end
    return curvature
end

let 
    t1 = 1.0
    t2 = 0.0
    m2 = 3√3 * 0.2 * t1
    num = 120
    kxs = range(-4π/3, 0, num + 1)
    kys = range(0, 4π/3, num + 1)
    dkx, dky = step(kxs), step(kys)
    kx_grid = kxs .+ 1/2 * kys'
    ky_grid = zeros(num+1) .+ √3/2 * kys'

    lowerband, upperband, lowervecs = eigenHaldane(kx_grid, ky_grid, t1, t2; m2=m2)
    curvature = computeBerryCurvature(lowervecs, dkx * dky)
    chernnum = sum(curvature) * dkx * dky / (2π)
    println("Chern number: ", chernnum)

    h5open("solidstatehomework/hw2/haldane_data_2.h5", "w") do file
        write(file, "kx_grid", kx_grid)
        write(file, "ky_grid", ky_grid)
        write(file, "lowerband", lowerband)
        write(file, "upperband", upperband)
        write(file, "lowervecs", lowervecs)
        write(file, "curvature", curvature)
        write(file, "chernnum", chernnum)
    end
end