# ASTRA_ORACLE: local
# ASTRA_EST_RUNTIME: short
"""
Verification of the anomalous projective algebra  {Pi, T_{N/2}} = 0  on a
tight-binding ring, and its three consequences:
  (1) explicit symmetry operator S = Pi * T^{N/2} exists (S^2 = -I) ONLY for N = 2 (mod 4);
  (2) prediction for ODD N: the operator cannot be built and Pi is not chiral -> symmetry ABSENT;
  (3) restoring the Peierls phase destroys the (antiunitary) symmetry C = S*K.
Legs: (a) exact integer matrices via SymPy; (b) numeric on random states + exact diagonalization;
      (c) limit/degenerate cases N=2, N=4; plus a Z3 proof of the mod-4 arithmetic core.
"""
import numpy as np
import sympy as sp
from z3 import Int, Solver, Not, sat, unsat

TOL = 1e-9
GRID = [6, 8, 10, 12]          # the paper's lattice
TESTSET = [2, 4, 6, 8, 10, 12, 14]
I_ = lambda N: np.eye(N)

# ---------- object constructors ----------
def Tmat(N):                    # translation |n> -> |n+1 mod N>
    T = np.zeros((N, N)); [T.__setitem__(((n + 1) % N, n), 1.0) for n in range(N)]
    return T
def Pimat(N):                   # sublattice / chiral operator diag((-1)^n)
    return np.diag([(-1.0) ** n for n in range(N)])
def Tpow(N):                    # T^{N/2} (integer power only for even N)
    return np.linalg.matrix_power(Tmat(N), N // 2)
def H0(N):                      # bare real NN ring, t = 1
    H = np.zeros((N, N))
    for n in range(N):
        H[(n + 1) % N, n] -= 1.0; H[n, (n + 1) % N] -= 1.0
    return H
def Hphi(N, th):                # Peierls-threaded ring: bond n->n+1 carries e^{i th}
    H = np.zeros((N, N), complex)
    for n in range(N):
        H[(n + 1) % N, n] += -np.exp(1j * th); H[n, (n + 1) % N] += -np.exp(-1j * th)
    return H
def acomm(A, B): return A @ B + B @ A

# ---------- leg (a): EXACT integer matrices (no floating point) ----------
def sym_bracket_zero(N):        # exact SymPy: is {Pi, T^{N/2}} the literal zero matrix?
    T = sp.zeros(N)
    for n in range(N): T[(n + 1) % N, n] = 1
    Tp = T ** (N // 2)
    P = sp.diag(*[(-1) ** n for n in range(N)])
    return bool((P * Tp + Tp * P).is_zero_matrix)

check_sym = all(sym_bracket_zero(N) == (N % 4 == 2) for N in TESTSET if N % 2 == 0)

# ---------- leg (b1): numeric identity on random states (independent code path) ----------
rng = np.random.default_rng(7)
check_num = True
for N in TESTSET:
    if N % 2:      # odd N: T^{N/2} is not an integer power -> skip here, handled below
        continue
    P, Tp = Pimat(N), Tpow(N)
    v = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    is_zero = np.linalg.norm(acomm(P, Tp) @ v) < TOL
    check_num &= (is_zero == (N % 4 == 2))

# ---------- Z3 proof of the arithmetic core:  N=2m,  bracket=0 <=> m odd <=> N=2 (mod4) ----------
m = Int('m')
prop = ((m % 2 == 1) == ((2 * m) % 4 == 2))   # the equivalence used above, for all integers
s = Solver(); s.add(Not(prop)); res = s.check()
if res not in (sat, unsat):
    raise RuntimeError("Z3 returned 'unknown' -> operational, not a refutation")
check_z3 = (res == unsat)

# ---------- leg (c): limit / degenerate cases ----------
# smallest anomalous ring N=2 :  S^2 = -I ;   degenerate failing ring N=4 : S^2 = +I (bracket != 0)
S2 = Pimat(2) @ Tpow(2); S4 = Pimat(4) @ Tpow(4)
check_limit = (np.linalg.norm(S2 @ S2 + I_(2)) < TOL and
               np.linalg.norm(S4 @ S4 - I_(4)) < TOL and
               not sym_bracket_zero(4) and sym_bracket_zero(2))

# ---------- deliverable (3): Peierls phase destroys the antiunitary symmetry C = S*K ----------
# C H C^-1 = U conj(H) U^-1 with U = S (real, orthogonal); symmetry means C H C^-1 = -H.
N = 6; U = Pimat(N) @ Tpow(N); Uinv = np.linalg.inv(U)
trans = lambda H: U @ np.conjugate(H) @ Uinv
theta = 0.3
rel0 = np.linalg.norm(trans(H0(N)) + H0(N))                 # bare ring: exact symmetry
relθ = np.linalg.norm(trans(Hphi(N, theta)) + Hphi(N, theta))  # threaded: broken
matches_flip = np.linalg.norm(trans(Hphi(N, theta)) + Hphi(N, -theta)) < TOL  # C H(θ)C^-1=-H(-θ)
csq = np.linalg.norm((U @ np.conjugate(U)) + I_(N))         # C^2 = U U* = -I (Kramers) for N=6
check_peierls = (rel0 < TOL and relθ > 1e-3 and matches_flip and csq < TOL)

# ---------- deliverable (2): ODD N prediction -> symmetry ABSENT ----------
check_odd = True
for N in [5, 7, 9]:
    P = Pimat(N); T = Tmat(N)
    no_anticomm = all(np.linalg.norm(acomm(P, np.linalg.matrix_power(T, p))) > TOL
                      for p in range(1, N))           # no integer power anticommutes with Pi
    not_chiral = np.linalg.norm(acomm(P, H0(N))) > TOL  # frustrated (odd) ring: {Pi,H} != 0
    check_odd &= (no_anticomm and not_chiral)

# ---------- exact diagonalization cross-check on the paper grid ----------
check_diag = True
for N in GRID:
    E = np.sort(np.linalg.eigvalsh(H0(N)))
    chiral_pairing = np.allclose(E, -E[::-1], atol=TOL)   # E <-> -E from {Pi,H0}=0
    S = Pimat(N) @ Tpow(N)
    anomalous = np.linalg.norm(S @ S + I_(N)) < TOL       # S^2 = -I ?
    check_diag &= (chiral_pairing and (anomalous == (N % 4 == 2)))

# ---------- self-refutation report ----------
print("Symmetry operator:  S = Pi * T^{N/2},  Pi = diag((-1)^n),  T|n>=|n+1>")
print("Closed form: {Pi,T^{N/2}} = (1+(-1)^{N/2}) * (signed shift)  ->  0  iff  N = 2 (mod 4)")
print(f"CHECK sym_exact   : {'OK' if check_sym    else 'FAIL'}  (SymPy integer matrices, N in {TESTSET})")
print(f"CHECK num_random  : {'OK' if check_num    else 'FAIL'}  (random-state operator norm, seed=7)")
print(f"CHECK z3_mod4     : {'OK' if check_z3     else 'FAIL'}  (negation unsat: m odd <=> 2m=2 mod4)")
print(f"CHECK limit_2_4   : {'OK' if check_limit  else 'FAIL'}  (N=2: S^2=-I ; N=4: S^2=+I, bracket!=0)")
print(f"CHECK peierls_kill: {'OK' if check_peierls else 'FAIL'}  (rel0={rel0:.2e}, rel_theta={relθ:.3f}, C^2+I={csq:.2e})")
print(f"CHECK odd_absent  : {'OK' if check_odd    else 'FAIL'}  (N=5,7,9: no anticomm power & {{Pi,H0}}!=0)")
print(f"CHECK exact_diag  : {'OK' if check_diag   else 'FAIL'}  (grid {GRID}: E<->-E & S^2=-I iff N=2 mod4)")

all_ok = all([check_sym, check_num, check_z3, check_limit, check_peierls, check_odd, check_diag])
if all_ok:
    print("VERDICT: PASS  -- anomalous algebra {Pi,T_{N/2}}=0 exists exactly for N=2(mod4);")
    print("               absent for odd N; destroyed by the Peierls phase. (paper grid: only N=6,10 carry it)")
else:
    print("VERDICT: FAIL  -- a decisive check refuted the conjecture; offending flags:",
          {k: v for k, v in dict(sym=check_sym, num=check_num, z3=check_z3, limit=check_limit,
                                 peierls=check_peierls, odd=check_odd, diag=check_diag).items() if not v})