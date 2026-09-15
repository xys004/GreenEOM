"""Two questions a referee will ask, answered numerically.

(1) Compare two FLUX-BASED equilibrium expectation values:
        spin-flux route: -Tr(P dH/d(phi_s))
        flux anticommutator: -Tr(P {Sigma_z, dH/d(phi)}/2).
    These are not a test of the conventional all-bond spin current.
    In the phaseless model dH/d(phi) omits range-two hops, whereas the
    oriented-bond velocity includes them. The old manuscript incorrectly
    described the second quantity as the conventional bond current.
    The two flux operators differ on spin-flip hops, but that difference
    has zero equilibrium expectation in the clean block-diagonal model.
    See current_definition_audit.py for an analytic three-way comparison.

(2) The selection rule's crux: with no Peierls phase, V = dH/d(lambda_R2)
    vanishes IDENTICALLY on the k = +-pi/2 blocks of the generalised
    translation.  The paper asserts this.  Here we exhibit it AND identify the
    closed form that explains it:

        || V restricted to the block at k ||_F  =  sqrt(2) |sin(2k - 2 phi)|

    Phaseless (phi -> 0) this is sqrt(2)|sin 2k|, which vanishes at 2k = +-pi,
    i.e. exactly at k = +-pi/2 -- a momentum that the half-integer ladder
    k = pi(2l+1)/N contains only when N = 2 (mod 4).  Restoring the phase
    shifts the argument to 2k - 2phi and moves the zero off k = +-pi/2.
"""
import numpy as np
from selection_rule import famH
from trace_rule import V


def d_dphis(N, a, ph, pei, d=1e-6):
    return (famH(N, *a, ph, d, pei) - famH(N, *a, ph, -d, pei)) / (2 * d)


def d_dphi(N, a, ph, pei, d=1e-6):
    return (famH(N, *a, ph + d, 0.0, pei) - famH(N, *a, ph - d, 0.0, pei)) / (2 * d)


def occ_expect(H, Op, nf):
    w, U = np.linalg.eigh(H)
    P = U[:, :nf]
    return float(np.real(np.trace(P.conj().T @ Op @ P)))


def js_two_ways(N, a, ph, pei):
    H = famH(N, *a, ph, 0.0, pei)
    Sz = np.kron(np.eye(N), np.diag([1.0, -1.0])).astype(complex)
    Ops = d_dphis(N, a, ph, pei)
    Opc = d_dphi(N, a, ph, pei)
    Opa = (Sz @ Opc + Opc @ Sz) / 2.0
    return -occ_expect(H, Ops, N), -occ_expect(H, Opa, N)


print("=" * 76)
print(" (1) spin-flux derivative  vs  FLUX anticommutator (not all-bond current)")
print("=" * 76)
print("      N  phase  lEO    lR1    lR2    phi/phi0    -dE/dphi_s    flux-anticomm."
      "      diff")
worst = 0.0
for N in (6, 8, 10, 14):
    for pei in (0, 1):
        for (lEO, lR1, lR2, p) in ((0.039, 0.12, 0.15, 0.25),
                                   (0.2, 0.3, 0.3, 0.4),
                                   (0.0, 0.12, 0.3, 0.125)):
            ph = 2 * np.pi * p / N
            jf, ja = js_two_ways(N, (lEO, lR1, lR2), ph, pei)
            worst = max(worst, abs(jf - ja))
            print(f"     {N:2d}   {'on ' if pei else 'off'}  {lEO:5.3f} {lR1:5.3f}"
                  f" {lR2:5.3f}   {p:5.3f}    {jf:+.10f}  {ja:+.10f}  {jf-ja:+.1e}")
print(f"\n   worst disagreement over all cases: {worst:.2e}")

print("\n   why it is not trivial (N=10, lR2=0.15, phase off):")
N, a, ph = 10, (0.039, 0.12, 0.15), 2 * np.pi * 0.25 / 10
Sz = np.kron(np.eye(N), np.diag([1.0, -1.0])).astype(complex)
Ops, Opc = d_dphis(N, a, ph, 0), d_dphi(N, a, ph, 0)
Opa = (Sz @ Opc + Opc @ Sz) / 2.0
flip = np.ones_like(Ops, dtype=bool)
for i in range(2 * N):
    for j in range(2 * N):
        if (i % 2) == (j % 2):
            flip[i, j] = False
print(f"     the two OPERATORS differ by            {np.linalg.norm(Ops-Opa):.4e}")
print(f"       on spin-diagonal entries             {np.linalg.norm((Ops-Opa)*~flip):.4e}")
print(f"       on spin-flip entries                 {np.linalg.norm((Ops-Opa)*flip):.4e}")
print(f"     anticommutator on spin-flip entries    {np.linalg.norm(Opa*flip):.4e}")
print("     so the flip part is carried only by the flux route, and its")
print("     occupied expectation value vanishes.")

print()
print("=" * 76)
print(" (2) V = dH/dlambda_R2 on the blocks of the generalised translation")
print("=" * 76)


def Tgen(N):
    T = np.zeros((N, N))
    for n in range(N):
        T[(n + 1) % N, n] = 1.0
    chi = -2 * np.pi / N
    return np.kron(T, np.diag([np.exp(1j * chi / 2), np.exp(-1j * chi / 2)]))


def blocks(N, ph, pei):
    G, Vm = Tgen(N), V(N, ph, 0.0, pei)
    ev, W = np.linalg.eig(G)
    o = np.argsort(np.angle(ev))
    ev, W = ev[o], W[:, o]
    ang = np.angle(ev)
    out, i = [], 0
    while i < len(ang):
        j = i
        while j + 1 < len(ang) and abs(ang[j + 1] - ang[i]) < 1e-8:
            j += 1
        B = np.linalg.qr(W[:, list(range(i, j + 1))])[0]
        out.append((ang[i], np.linalg.norm(B.conj().T @ Vm @ B)))
        i = j + 1
    return out


for N in (10, 6, 14):
    for pei in (0, 1):
        ph = 2 * np.pi * 0.25 / N
        print(f"\n   N={N}, Peierls phase {'restored' if pei else 'omitted'}"
              f"   (phi = 2*pi*0.25/{N})")
        print("       k/pi     ||V block||    sqrt(2)|sin(2k-2phi)|   ratio")
        for k, nrm in blocks(N, ph, pei):
            pred = np.sqrt(2) * abs(np.sin(2 * k - 2 * (ph if pei else 0.0)))
            r = nrm / pred if pred > 1e-12 else float('nan')
            tag = "  <-- k=+-pi/2" if abs(abs(k) - np.pi / 2) < 1e-8 else ""
            print(f"     {k/np.pi:+7.4f}   {nrm:.6e}   {pred:.6e}   "
                  f"{r:8.5f}{tag}")

print("\n   the zero at k=+-pi/2 is destroyed by the phase as sqrt(2)|sin 2phi|:")
print("      phi/phi0   ||V block at k=+pi/2||   sqrt(2)|sin 2phi|")
for p in (0.0, 0.05, 0.125, 0.25, 0.5):
    ph = 2 * np.pi * p / 10
    got = [n for k, n in blocks(10, ph, 1) if abs(k - np.pi / 2) < 1e-8][0]
    print(f"       {p:5.3f}      {got:.6e}          "
          f"{np.sqrt(2)*abs(np.sin(2*ph)):.6e}")
