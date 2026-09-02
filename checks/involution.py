"""Which state pairs with which, and what relation do their energies satisfy?
The verified statement is that the multiset {V_ii : i occupied} is symmetric
under negation. That asks for an involution sigma on the occupied states with
V_{sigma(i)} = -V_i. Find it, then look at the energies it relates."""
import numpy as np
from selection_rule import famH
from trace_rule import V

def data(N, lEO=0.039, lR1=0.12, lR2=0.15, phi=0.25, pei=0):
    ph = 2*np.pi*phi/N
    w, U = np.linalg.eigh(famH(N, lEO, lR1, lR2, ph, 0.0, pei))
    d = np.real(np.diag(U.conj().T @ V(N, ph, 0.0, pei) @ U))
    return w, d, U

for N in (6, 10, 14):
    w, d, U = data(N)
    occ = list(range(N))
    used, pairs = set(), []
    for i in occ:
        if i in used: continue
        cand = [j for j in occ if j not in used and j != i and abs(d[j] + d[i]) < 1e-9]
        if abs(d[i]) < 1e-12:
            used.add(i); pairs.append((i, i)); continue
        j = cand[0]; used |= {i, j}; pairs.append((i, j))
    print(f"N={N}  (espectro completo simetrico E<->-E?  "
          f"{np.abs(np.sort(w)+np.sort(w)[::-1]).max():.2e})")
    for (i, j) in pairs:
        print(f"   i={i:2d} E={w[i]:+8.5f} V={d[i]:+8.5f}   <->   "
              f"j={j:2d} E={w[j]:+8.5f} V={d[j]:+8.5f}    E_i+E_j={w[i]+w[j]:+8.5f}")
    print()
