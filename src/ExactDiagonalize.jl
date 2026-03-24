module ExactDiagonalize

using MKL, LinearAlgebra
import SparseArrays: spzeros
import SparseArrays: SparseMatrixCSC as SpMatrix

include("operators.jl")
include("observers.jl")
include("exactdiag.jl")
include("ode_solver.jl")
include("sparsemat.jl")

# export public API
export NumBasis, FullBasis, AbstractState, NumState, FullState, State
export set_systype, get_systype, AbstractOp, SpinOp, get_optype, OpSum, makeHamiltonian
export apply, apply!, expected, inner
export AbstractObserver, OperatorObserver, OpSumObserver, ZObserver, XObserver, record!
export spectrum, exact, rk4, spmat, timeEvolve

end