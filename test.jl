include("operators.jl")
include("exactdiagonalize.jl")
include("ode_solver.jl")

using CairoMakie

let 
    L, N = 10, 1
    Δ = 0.5
    init = NumState("1000000000")
    
    os = Tuple[]
    for j in 1:L
        nj = mod1(j + 1, L)
        push!(os, (Δ, :Z, j, :Z, nj))
        push!(os, (1.0, :X, j, :X, nj))
        push!(os, (-1.0, :iY, j, :iY, nj))
    end
    # push!(os, (1.0, :X, L))
    ops = SpinOpSum(Float64, os)

    obs = OperatorObserver((1.0, :Z, L), init.basis)
    
    #makeHamiltonian(ops, init.basis)
    
    ts = 0.0:0.02:10.0
    timeEvolve_exact(ops, init, ts, obs)
    fig, ax = lines(ts, obs.data)
    
end
