"""Does the projective-algebra lemma explain the selection rule?

A lemma on the SPINLESS ring (projective_algebra.py) shows: S = Pi T^{N/2} with Pi = diag((-1)^n)
obeys {Pi, T^{N/2}} = 0 exactly for N = 2 (mod 4). The mod-4 structure is right,
but our statement is about the 2N-dimensional SPINFUL ring with lambda_R1 and
lambda_R2 present, and about J_s^z = -dE/d(phi_s) being independent of
lambda_R2 -- which is stronger than E being even in it.

This script tests the connection rather than assuming it.
"""
import numpy as np
import itertools

def famH(N, lEO, lR1, lR2, ph, phs, pei):
    """Same H as scans.wl / figure4.wl, canonical bond angle."""
    d = 2 * N
    H = np.zeros((d, d), complex)
    ix = lambda n, s: 2 * n + (0 if s == 1 else 1)
    th = lambda n, m: 2 * np.pi * n / N + m * np.pi / N
    def add(i, j, v):
        H[i, j] += v; H[j, i] += np.conj(v)
    for n in range(N):
        for s in (1, -1):
            nx, n2 = (n + 1) % N, (n + 2) % N
            p1 = np.exp(1j * (ph + s * phs))
            p2 = np.exp(2j * (ph + s * phs)) if pei else 1.0
            add(ix(n, s), ix(nx, s), p1)
            add(ix(n, s), ix(nx, -s), -(lR1 / 2) * s * np.exp(-1j * s * th(n, 1)) * p1)
            add(ix(n, s), ix(n2, s), 1j * lEO * (s / 2) * p2)
            add(ix(n, s), ix(n2, -s), -(lR2 / 2) * s * np.exp(-1j * s * th(n, 2)) * p2)
    return H

def occ(N, *a, nf=None):
    e = np.sort(np.linalg.eigvalsh(famH(N, *a)).real)
    return e[: (nf if nf is not None else N)].sum()

def js(N, lEO, lR1, lR2, ph, pei, d=5e-4):
    return -(occ(N, lEO, lR1, lR2, ph, d, pei) - occ(N, lEO, lR1, lR2, ph, -d, pei)) / (2 * d)

# --- 1. is E even in lambda_R2, or is J_s genuinely independent of it? -------
print("1. E(lambda_R2) vs E(-lambda_R2), and J_s across the range")
for N in (6, 8, 10):
    ph = 2 * np.pi * 0.25 / N
    ev = max(abs(occ(N, 0.039, 0.12, x, ph, 0.0, 0) - occ(N, 0.039, 0.12, -x, ph, 0.0, 0))
             for x in (0.05, 0.15, 0.3))
    j0 = js(N, 0.039, 0.12, 0.0, ph, 0)
    jspan = max(abs(js(N, 0.039, 0.12, x, ph, 0) - j0) for x in (0.05, 0.15, 0.3))
    espan = max(abs(occ(N, 0.039, 0.12, x, ph, 0.0, 0) - occ(N, 0.039, 0.12, 0.0, ph, 0.0, 0))
                for x in (0.05, 0.15, 0.3))
    print(f"   N={N:2d} (N%4={N%4})  |E(x)-E(-x)|max={ev:.2e}   "
          f"|E(x)-E(0)|max={espan:.3e}   |J_s(x)-J_s(0)|max={jspan:.2e}")

# --- 2. does ASTRA's S, lifted to spin, act on our H? -----------------------
print()
print("2. candidate operators U with U H U^-1 = H(lambda_R2 -> -lambda_R2)?")
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1, -1]).astype(complex)
id2 = np.eye(2, dtype=complex)

def lift(N, spinop):
    T = np.zeros((N, N)); 
    for n in range(N): T[(n + 1) % N, n] = 1.0
    Tp = np.linalg.matrix_power(T, N // 2)
    Pi = np.diag([(-1.0) ** n for n in range(N)])
    return np.kron(Pi @ Tp, spinop)

for N in (6, 8, 10):
    ph = 2 * np.pi * 0.25 / N
    H_p = famH(N, 0.039, 0.12, 0.15, ph, 0.0, 0)
    H_m = famH(N, 0.039, 0.12, -0.15, ph, 0.0, 0)
    out = []
    for name, sop in (("1", id2), ("sx", sx), ("sy", sy), ("sz", sz)):
        U = lift(N, sop)
        r = np.linalg.norm(U @ H_p @ np.linalg.inv(U) - H_m) / np.linalg.norm(H_m)
        rc = np.linalg.norm(U @ H_p.conj() @ np.linalg.inv(U) - H_m) / np.linalg.norm(H_m)
        if r < 1e-10 or rc < 1e-10:
            out.append(f"{name}{'(K)' if rc < 1e-10 else ''}")
    print(f"   N={N:2d}  operators that work: {out if out else 'NONE'}")

# --- 3. the honest question: does ANY of this give J_s independence? --------
print()
print("3. what the spinless lemma does and does not settle")
print("   E even in lambda_R2 would give dJ_s/d(lambda_R2)=0 at 0, not independence")
print("   over a range. Section 1 above measures which of the two we actually have.")
