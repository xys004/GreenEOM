"""Exact symbolic confirmation of the selection-rule derivation.

Independent of the numpy checks in this directory (floating point) and of the
Mathematica package (which produces the paper's numbers). The linchpin is that

    T = T_1 (x) exp(-i sigma_z pi/N)

commutes with H and with V = dH/d(lambda_R2) as an IDENTITY in the couplings
lEO, lR1, lR2 and the flux ph -- not a sampled numerical fact -- which is what
validates the block-diagonalisation in generalised momentum. Also verified:
T^N = -1, and that pi/2 is an allowed momentum iff N = 2 (mod 4).

Run locally with `python exact_symmetry.py`. Also run on the ASTRUM cluster
(SymPy 1.14) as a second machine / second engine: 13/13 exact checks, PASS.
"""
import sympy as sp

I = sp.I
lEO, lR1, lR2, ph = sp.symbols('lEO lR1 lR2 ph', real=True)

def theta(N, n, m): return sp.Rational(2, 1)*sp.pi*n/N + m*sp.pi/N
def ixf(n, s): return 2*n + (0 if s == 1 else 1)

def build(N, lR2v):
    d = 2*N
    H = sp.zeros(d, d)
    def add(i, j, v):
        H[i, j] += v
        H[j, i] += sp.conjugate(v)
    for n in range(N):
        for s in (1, -1):
            nx, n2 = (n+1) % N, (n+2) % N
            add(ixf(n, s), ixf(nx, s), sp.exp(I*ph))
            add(ixf(n, s), ixf(nx, -s), -sp.Rational(1,2)*lR1*s*sp.exp(-I*s*theta(N, n, 1))*sp.exp(I*ph))
            add(ixf(n, s), ixf(n2, s), I*lEO*sp.Rational(s,2))
            add(ixf(n, s), ixf(n2, -s), -sp.Rational(1,2)*lR2v*s*sp.exp(-I*s*theta(N, n, 2)))
    return H

def Tgen(N):
    d = 2*N
    T = sp.zeros(d, d)
    up, dn = sp.exp(-I*sp.pi/N), sp.exp(I*sp.pi/N)
    for n in range(N):
        m = (n+1) % N
        T[ixf(m, 1),  ixf(n, 1)]  = up
        T[ixf(m, -1), ixf(n, -1)] = dn
    return T

def iszero(M):
    return sp.simplify(sp.expand_complex(M)).is_zero_matrix

results = []
def check(name, ok):
    results.append((name, bool(ok)))
    print(("  OK  " if ok else " FAIL ") + name, flush=True)

for N in (6, 8, 10):
    H = build(N, lR2)
    V = build(N, sp.Integer(1)) - build(N, sp.Integer(0))
    check(f"N={N}: V == dH/d(lR2)", iszero(V - sp.diff(H, lR2)))
    T = Tgen(N)
    check(f"N={N}: T^N = -1", iszero(T**N + sp.eye(2*N)))
    check(f"N={N}: [T,H] = 0  (identity in lEO,lR1,lR2,ph)", iszero(T*H - H*T))
    check(f"N={N}: [T,V] = 0", iszero(T*V - V*T))

def half_pi_allowed(N):
    return any(sp.simplify(sp.pi*(2*l+1)/N - sp.pi/2) == 0 for l in range(N))
check("k = pi/2 allowed  <=>  N = 2 (mod 4)",
      all(half_pi_allowed(N) == (N % 4 == 2) for N in range(3, 41)))

nfail = sum(1 for _, ok in results if not ok)
print(f"\n{len(results)-nfail}/{len(results)} exact checks passed")
print("VERDICT:", "PASS" if nfail == 0 else "FAIL")
