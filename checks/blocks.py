"""Tr_k V = 0 in every 2x2 block, for every N. So sum_occ V_ii = 0 exactly when
the occupied set is a union of COMPLETE blocks -- i.e. no block is split by the
Fermi level. Count the split blocks."""
import numpy as np
from selection_rule import famH
from trace_rule import V

def Tgen(N, chi):
    T = np.zeros((N, N))
    for n in range(N): T[(n+1) % N, n] = 1.0
    return np.kron(T, np.diag([np.exp(1j*chi/2), np.exp(-1j*chi/2)]))

def split_blocks(N, lEO, lR1, lR2, phi, pei):
    ph = 2*np.pi*phi/N
    H = famH(N, lEO, lR1, lR2, ph, 0.0, pei); Vm = V(N, ph, 0.0, pei)
    G = Tgen(N, -2*np.pi/N)
    comm = np.linalg.norm(G@H - H@G)/np.linalg.norm(H)
    ev, W = np.linalg.eig(G); o = np.argsort(np.angle(ev)); ev, W = ev[o], W[:, o]
    ph_ = np.angle(ev); groups=[]; i=0
    while i < len(ph_):
        j=i
        while j+1<len(ph_) and abs(ph_[j+1]-ph_[i])<1e-8: j+=1
        groups.append(list(range(i,j+1))); i=j+1
    w = np.sort(np.linalg.eigvalsh(H).real); Ef = (w[N-1]+w[N])/2
    split = 0; trmax = 0.0
    for g in groups:
        B = np.linalg.qr(W[:, g])[0]
        eb = np.linalg.eigvalsh(B.conj().T @ H @ B).real
        trmax = max(trmax, abs(np.trace(B.conj().T @ Vm @ B).real))
        below = int((eb < Ef).sum())
        if 0 < below < len(g): split += 1
    return comm, trmax, split, len(groups)

print("bloques partidos por el nivel de Fermi a semi-llenado")
print()
print("  caso SIN fase de Peierls")
print("    N  N%4   [T,H]     max|Tr_k V|   bloques partidos / total")
for N in (5,6,7,8,9,10,11,12,13,14,18):
    c,t,s,g = split_blocks(N, 0.039, 0.12, 0.15, 0.25, 0)
    print(f"   {N:3d}   {N%4}   {c:.0e}   {t:.1e}        {s} / {g}"
          + ("     <-- ninguno" if s==0 else ""))
print()
print("  caso CON fase restaurada")
print("    N  N%4   [T,H]     max|Tr_k V|   bloques partidos / total")
for N in (6,10,14):
    c,t,s,g = split_blocks(N, 0.039, 0.12, 0.15, 0.25, 1)
    print(f"   {N:3d}   {N%4}   {c:.0e}   {t:.1e}        {s} / {g}")
