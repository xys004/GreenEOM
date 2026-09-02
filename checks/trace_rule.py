"""Hellmann-Feynman: dE_occ/d(lambda_R2) = sum over occupied i of V_ii, with
V = dH/d(lambda_R2) and V_ii its diagonal in the eigenbasis of H. So the whole
question is the pattern of those diagonal elements."""
import numpy as np
from selection_rule import famH

def V(N, ph, phs, pei):
    """dH/d(lambda_R2): the range-two Rashba term at unit coupling."""
    d = 2 * N; M = np.zeros((d, d), complex)
    ix = lambda n, s: 2 * n + (0 if s == 1 else 1)
    th = lambda n, m: 2 * np.pi * n / N + m * np.pi / N
    for n in range(N):
        for s in (1, -1):
            p2 = np.exp(2j * (ph + s * phs)) if pei else 1.0
            v = -(0.5) * s * np.exp(-1j * s * th(n, 2)) * p2
            i, j = ix(n, s), ix((n + 2) % N, -s)
            M[i, j] += v; M[j, i] += np.conj(v)
    return M

print("V_ii en la base propia de H (lEO=0.039, lR1=0.12, lR2=0.15, phi/phi0=1/4)")
print("sin fase de Peierls:")
print()
for N in (6, 8, 10, 14):
    ph = 2 * np.pi * 0.25 / N
    H = famH(N, 0.039, 0.12, 0.15, ph, 0.0, 0)
    w, U = np.linalg.eigh(H)
    Vd = np.real(np.diag(U.conj().T @ V(N, ph, 0.0, 0) @ U))
    occ, emp = Vd[:N], Vd[N:]
    print(f"  N={N:2d} (N%4={N%4})")
    print(f"     V_ii ocupados : {np.array2string(occ, precision=4, suppress_small=True, max_line_width=200)}")
    print(f"     suma ocupada  : {occ.sum():+.3e}      suma vacia: {emp.sum():+.3e}     traza total: {Vd.sum():+.3e}")
    # is the occupied list antisymmetric about its middle?
    mir = np.abs(occ + occ[::-1]).max()
    print(f"     |V_ii + V_(N+1-i,N+1-i)| max sobre ocupados = {mir:.2e}"
          f"{'   <-- antisimetrica dentro del bloque ocupado' if mir < 1e-9 else ''}")
    print()
