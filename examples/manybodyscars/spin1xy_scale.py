from spin1xy import *

Ls = [10, 16, 24, 32]
Stot = 2
J, h = 1.0, 1.0
nt = 501
ts = np.geomspace(0.1, 1e8, nt)

dts = np.geomspace(1e-2, 1e4, nt)
profile = [11 if i < (nt // 5) else 11 for i in range(nt)]
profile[-1] = 51
ts_full = time_expand(ts, dts, profile)

entropies = np.zeros((nt, len(Ls)))
entropies_full = np.zeros_like(ts_full)

for n, L in enumerate(Ls):
    b = L // 2
    basis = spin_basis_1d(L=L, S="1", Nup = Stot)
    E, U = spin1xy_spectrum(basis, L, J, h)
    print("L = {} spectrum solved".format(L))
    
    psi0 = np.zeros(basis.Ns)
    psi0[basis.Ns // 2] = 1.0
    
    psi_t = ED_state_vs_time(psi0, E, U, ts_full, iterate=True)
    
    subA = tuple(range(L // 2))
    for i, psi in enumerate(psi_t):
        entr = my_ent_entropy(basis.states, 3, psi, b, density=False)
        entropies_full[i] = entr
        
    entropies[:, n] = latetime_average(entropies_full, profile)
    print("L = {} entropy obtained".format(L))

np.savez(f"examples/manybodyscars/spin1xy_Sz={Stot}_Lmax={Ls[-1]}.npz", 
        Ls = np.array(Ls), ts = ts,
        params = np.array([J, h]),  
        entropies=entropies
    )