import numpy as np
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian, quantum_LinearOperator
import matplotlib.pyplot as plt

def pxp_spectrum(lsize:int, g):
    basis = spin_basis_1d(lsize, a=2, kblock=0, pblock=1)
    
    

