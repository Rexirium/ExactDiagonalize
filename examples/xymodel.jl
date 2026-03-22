#= Load ExactDiagonalize package if not already loaded
if !isdefined(Main, :ExactDiagonalize)
    include("../src/ExactDiagonalize.jl")
    using .ExactDiagonalize
end
=#
using Revise
using ExactDiagonalize
# For plotting
using CairoMakie

let 
    set_systype(:Spin)  # Set system type to spin
    L, N = 10, 1        # System size and particle number
    Δ = 1.0             # Interaction parameter

    # Initial state: single up spin at site 1
    init = NumState("1000000000")
    
    # Build Hamiltonian terms for XY model
    opsum = OpSum(Float64)
    for j in 1:L
        nj = mod1(j + 1, L)
        opsum += (Δ, :Z, j, :Z, nj)
        opsum += (1.0, :X, j, :X, nj)
        opsum += (-1.0, :iY, j, :iY, nj)
    end

    # opsum2: sum of Z operators (not used here)
    os2 = [(1.0, :Z, j) for j in 1:L]
    opsum2 = OpSum(os2, Float64)

    # Observable: Z at last site
    obs = OperatorObserver((1.0, :Z, L), init.basis)
    
    # Time points for evolution
    ts = 0.0:0.05:10.0
    @time timeEvolve(opsum, init, ts, obs)

    # Plot observable vs time
    fig = Figure()
    ax = Axis(fig[1,1],

    )
    lines!(ax, ts, obs.data)
    fig
end
