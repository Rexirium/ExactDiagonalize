
function timeEvolve_rk4(ops::AbstractOpSum, init::AbstractState, ts::AbstractVector, obs::AbstractObserver)
    hmat = makeHamiltonian(ops, init.basis; sparsed=true)
    psi = copy(init.vector)
    dim = length(psi)
    
    k1 = Vector{ComplexF64}(undef, dim)
    k2 = Vector{ComplexF64}(undef, dim)
    k3 = Vector{ComplexF64}(undef, dim)
    k4 = Vector{ComplexF64}(undef, dim)
    tmp = Vector{ComplexF64}(undef, dim)

    for (i, t) in enumerate(ts)
        record!(obs, psi)
        if i == length(ts)
            break
        end
        h = ts[i+1] - t
        h_2 = h / 2

        k1 .= -im * (hmat * psi)
        @. tmp = psi + h_2 * k1
        k2 .= -im * hmat * tmp
        @. tmp = psi + h_2 * k2
        k3 .= -im * hmat * tmp
        @. tmp = psi + h * k3
        k4 .= -im * hmat * tmp

        @. psi += (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    end
    return State(init.basis, psi)

end
