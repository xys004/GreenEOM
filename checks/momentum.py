"""Hypothesis: H and V both commute with a generalised translation
   Tgen = T_1 (x) exp(i chi sigma_z / 2),
the Rashba ring's rotating-frame symmetry. If so both are block diagonal in its
eigenvalue (a generalised momentum k), each block 2x2, and the +- pairing of the
V_ii is just Tr(V restricted to a k block) = 0. The occupied sum then vanishes
iff the occupied set is a union of COMPLETE k blocks -- which is where N mod 4
should enter."""
import numpy as np
from selection_rule import famH
from trace_rule import V

def Tgen(N, chi):
    T = np.zeros((N, N))
    for n in range(N): T[(n+1) % N, n] = 1.0
    R = np.diag([np.exp(1j*chi/2), np.exp(-1j*chi/2)])
    return np.kron(T, R)

print("1. buscar chi tal que [Tgen, H] = 0   (phaseless)")
for N in (6, 8, 10):
    ph = 2*np.pi*0.25/N
    H = famH(N, 0.039, 0.12, 0.15, ph, 0.0, 0)
    Vm = V(N, ph, 0.0, 0)
    best = None
    for m in (-2, -1, 0, 1, 2):
        chi = m*2*np.pi/N
        G = Tgen(N, chi)
        rh = np.linalg.norm(G@H - H@G)/np.linalg.norm(H)
        rv = np.linalg.norm(G@Vm - Vm@G)/max(np.linalg.norm(Vm), 1e-30)
        if best is None or rh < best[1]: best = (m, rh, rv)
    m, rh, rv = best
    print(f"   N={N:2d}  mejor m={m:+d} (chi={m}*2pi/N):  ||[T,H]||={rh:.2e}   ||[T,V]||={rv:.2e}")

print()
print("2. si conmutan: bloques de k, traza de V por bloque, y ocupacion")
for N in (6, 8, 10, 14):
    ph = 2*np.pi*0.25/N
    H = famH(N, 0.039, 0.12, 0.15, ph, 0.0, 0); Vm = V(N, ph, 0.0, 0)
    chi = -2*np.pi/N
    G = Tgen(N, chi)
    if np.linalg.norm(G@H - H@G)/np.linalg.norm(H) > 1e-10:
        chi = 2*np.pi/N; G = Tgen(N, chi)
    r = np.linalg.norm(G@H - H@G)/np.linalg.norm(H)
    # simultaneous block structure: diagonalise G, group by eigenvalue
    ev, W = np.linalg.eig(G)
    order = np.argsort(np.angle(ev)); ev, W = ev[order], W[:, order]
    Hk = W.conj().T @ np.linalg.inv(W.conj().T @ W) @ W.conj().T @ H @ W  # ill-conditioned; use QR per block instead
    # cleaner: group indices by eigenvalue phase, orthonormalise each group
    phases = np.angle(ev); groups = []
    i = 0
    while i < len(phases):
        j = i
        while j+1 < len(phases) and abs(phases[j+1]-phases[i]) < 1e-8: j += 1
        groups.append(list(range(i, j+1))); i = j+1
    traces, sizes = [], []
    ok = True
    for g in groups:
        B = np.linalg.qr(W[:, g])[0]
        hb = B.conj().T @ H @ B
        if np.linalg.norm(B @ hb @ B.conj().T @ B - H @ B)/np.linalg.norm(H) > 1e-8: ok = False
        traces.append(np.real(np.trace(B.conj().T @ Vm @ B))); sizes.append(len(g))
    print(f"   N={N:2d} (N%4={N%4})  ||[T,H]||={r:.1e}  bloques={len(groups)} de tamanos {sorted(set(sizes))}"
          f"  max|Tr_k V|={max(abs(np.array(traces))):.2e}")
