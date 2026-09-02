"""Every N has split blocks. What differs is whether the occupied member of a
split block carries V_ii = 0. In a 2x2 block V is traceless, so V_block = a.sigma
in some basis: its diagonal in the H eigenbasis vanishes iff the Bloch vector of
V is orthogonal to that of H. Which k are split, and where is the orthogonality?"""
import numpy as np
from selection_rule import famH
from trace_rule import V

def blocks(N, lEO, lR1, lR2, phi, pei):
    ph = 2*np.pi*phi/N
    H = famH(N, lEO, lR1, lR2, ph, 0.0, pei); Vm = V(N, ph, 0.0, pei)
    T = np.zeros((N, N))
    for n in range(N): T[(n+1) % N, n] = 1.0
    G = np.kron(T, np.diag([np.exp(-1j*np.pi/N), np.exp(1j*np.pi/N)]))
    ev, W = np.linalg.eig(G); o = np.argsort(np.angle(ev)); ev, W = ev[o], W[:, o]
    a = np.angle(ev); out=[]; i=0
    w = np.sort(np.linalg.eigvalsh(H).real); Ef = (w[N-1]+w[N])/2
    while i < len(a):
        j=i
        while j+1<len(a) and abs(a[j+1]-a[i])<1e-8: j+=1
        B = np.linalg.qr(W[:, i:j+1])[0]
        hb = B.conj().T @ H @ B; vb = B.conj().T @ Vm @ B
        eb, ub = np.linalg.eigh(hb)
        vdiag = np.real(np.diag(ub.conj().T @ vb @ ub))
        below = int((eb < Ef).sum())
        out.append(dict(k=a[i], split=(0 < below < len(eb)), E=eb,
                        vdiag=vdiag, trV=np.real(np.trace(vb)),
                        normV=np.linalg.norm(vb - np.trace(vb)/2*np.eye(2))))
        i=j+1
    return out

for N in (6, 8, 10):
    print(f"N={N} (N%4={N%4})")
    for b in blocks(N, 0.039, 0.12, 0.15, 0.25, 0):
        mark = " PARTIDO" if b['split'] else ""
        print(f"   k={b['k']:+7.4f} = {b['k']/np.pi:+6.3f} pi   "
              f"V_ii = [{b['vdiag'][0]:+8.5f} {b['vdiag'][1]:+8.5f}]   "
              f"||V_traceless||={b['normV']:.2e}{mark}")
    print()
