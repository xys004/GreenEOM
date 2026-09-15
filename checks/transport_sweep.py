"""Landauer transport through the open ring, in Python.

Two jobs.

  (1) VALIDATE.  The transport numbers reported so far come from the Wolfram
      implementation.  This is an independent code path in Python, checked
      against the deposited data/family_transmission.csv before anything new
      is computed with it.  If the two disagree, nothing below is worth
      reading.

  (2) SWEEP.  The claim that a conductance measurement separates the buckled
      ring from the flat one currently rests on a single parameter point
      (N = 6, Gamma = 0.3 t, omega = 0.3 t).  Here it is swept over contact
      strength, injection energy and ring size, and then over the coupling
      along the specified two-coupling path, without material calibration.

Conventions follow the manuscript and open_family.wl: wide-band reservoirs on
diametrically opposite sites, self-energy Sigma = -(i/2) (Gamma_L + Gamma_R),
retarded Green function G = (omega + i eta - H - Sigma)^-1, and

    M[out, in] = Re Tr[ (Gamma P_out) G (Gamma P_in) G^dag ],

with the charge transmission the sum over both spin labels and the transmitted
polarisation (T_up - T_down) / T_total.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io as _io, contextlib as _ctx
with _ctx.redirect_stdout(_io.StringIO()):
    from selection_rule import famH

ETA = 1e-8


def _proj(N, n, s):
    """Projector onto site n (0-indexed), spin s = +-1."""
    P = np.zeros((2 * N, 2 * N), complex)
    i = 2 * n + (0 if s == 1 else 1)
    P[i, i] = 1.0
    return P


def tmat(N, lEO, lR1, lR2, ph, pei, gam, w):
    """The 2x2 spin-resolved transmission matrix M[out, in]."""
    H = famH(N, lEO, lR1, lR2, ph, 0.0, pei)
    nr = N // 2                                   # site opposite site 0
    GL = gam * (_proj(N, 0, 1) + _proj(N, 0, -1))
    GR = gam * (_proj(N, nr, 1) + _proj(N, nr, -1))
    sig = -0.5j * (GL + GR)                       # wide band
    G = np.linalg.inv((w + 1j * ETA) * np.eye(2 * N) - H - sig)
    Gd = G.conj().T
    M = np.zeros((2, 2))
    for a, so in enumerate((1, -1)):
        for b, si in enumerate((1, -1)):
            M[a, b] = np.real(np.trace(
                (gam * _proj(N, nr, so)) @ G @ (gam * _proj(N, 0, si)) @ Gd))
    return M


def tchg(*a):
    return float(tmat(*a).sum())


def pz(*a):
    M = tmat(*a)
    return float((M[0, 0] + M[0, 1] - M[1, 0] - M[1, 1]) / M.sum())


# ---------------------------------------------------------------------------
# (1) validation against the deposited Wolfram output
# ---------------------------------------------------------------------------
def validate():
    import csv
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(here, "data", "family_transmission.csv")
    N, gam = 6, 0.3
    ph = 2 * np.pi * 0.25 / N
    gra = (0.2, 0.12, 0.0)
    sil = (0.2, 0.12, 0.15)

    worst = {"gra": 0.0, "sil": 0.0, "silNP": 0.0}
    rows = 0
    with open(path) as fh:
        for r in csv.DictReader(fh):
            w = float(r["w"])
            got = {
                "gra":   tchg(N, *gra, ph, 1, gam, w),
                "sil":   tchg(N, *sil, ph, 1, gam, w),
                "silNP": tchg(N, *sil, ph, 0, gam, w),
            }
            for k in worst:
                worst[k] = max(worst[k], abs(got[k] - float(r[k])))
            rows += 1

    print("=" * 74)
    print(" (1) Python transport vs the deposited Wolfram output")
    print("=" * 74)
    print("     N=6, Gamma=0.3t, phi/phi0=1/4, over %d injection energies"
          % rows)
    for k in ("gra", "sil", "silNP"):
        print("       %-6s worst |Python - Wolfram| = %.3e" % (k, worst[k]))
    ok = max(worst.values()) < 1e-9
    print("     %s" % ("AGREE to better than 1e-9, the sweep below is usable"
                       if ok else "DISAGREE -- stop, do not use the numbers below"))
    return ok


# ---------------------------------------------------------------------------
# (2a) does the separation survive contact strength, energy and ring size?
# ---------------------------------------------------------------------------
def sweep_geometry():
    print()
    print("=" * 74)
    print(" (2a) separation between phaseless and restored, over Gamma, omega, N")
    print("=" * 74)
    print("     lEO=0.2t, lR1=0.12t, lR2=0.15t, phi/phi0=1/4")
    print()
    print("      N   Gamma      min|dT|     max|dT|     min|dPz|    max|dPz|")
    lEO, lR1, lR2 = 0.2, 0.12, 0.15
    out = []
    for N in (6, 8, 10):
        ph = 2 * np.pi * 0.25 / N
        for gam in (0.1, 0.3, 0.6, 1.0):
            dT, dP = [], []
            for w in np.linspace(-2.0, 2.0, 81):
                a = (N, lEO, lR1, lR2, ph)
                dT.append(abs(tchg(*a, 0, gam, w) - tchg(*a, 1, gam, w)))
                dP.append(abs(pz(*a, 0, gam, w) - pz(*a, 1, gam, w)))
            print("     %2d   %4.2f    %.3e   %.3e   %.3e   %.3e"
                  % (N, gam, min(dT), max(dT), min(dP), max(dP)))
            out.append((N, gam, min(dT), max(dT), min(dP), max(dP)))
    return out


# ---------------------------------------------------------------------------
# (2b) a path varying both intrinsic couplings
#
# A guard first.  max|dPz| over all omega is NOT a usable statistic: Pz is a
# ratio, and where the transmission vanishes the ratio swings without meaning.
# Refining the grid shows it: at lR2 = 1e-4 the unrestricted maximum runs
# 0.044, 0.044, 0.144, 0.487 for 201, 801, 3201, 12801 points, still climbing.
# Restricting to energies that carry real current converges to five digits.
# ---------------------------------------------------------------------------
FRAC = 0.1          # keep omega where BOTH models transmit >= 10% of the peak


def separation(N, lEO, lR1, lR2, ph, gam, npts=12801, wmax=2.0):
    """Charge and polarisation separation between phaseless and restored."""
    a = (N, lEO, lR1, lR2, ph)
    ws = np.linspace(-wmax, wmax, npts)
    T0 = np.array([tchg(*a, 0, gam, w) for w in ws])
    T1 = np.array([tchg(*a, 1, gam, w) for w in ws])
    P0 = np.array([pz(*a, 0, gam, w) for w in ws])
    P1 = np.array([pz(*a, 1, gam, w) for w in ws])
    thr = FRAC * max(T0.max(), T1.max())
    keep = (T0 >= thr) & (T1 >= thr)
    return {
        "peakT": float(max(T0.max(), T1.max())),
        "dT": float(np.abs(T0 - T1).max()),
        "dPz": float(np.abs(P0 - P1)[keep].max()) if keep.any() else float("nan"),
    }


def sweep_path():
    print()
    print("=" * 74)
    print(" (2b) the separation along lEO = (4/3) lR2")
    print("=" * 74)
    print("     N=6, Gamma=0.3t, phi/phi0=1/4, lR1=0.12t, lEO scaled with lR2")
    print("     Illustrative 1D parameters; no mapping to a honeycomb material is assumed.")
    print("     polarisation taken where both models transmit >= %d%% of the peak"
          % (FRAC * 100))
    print()
    print("      lR2 (t)     peak T     max|dT|      max|dPz|     dPz/lR2")
    N, gam = 6, 0.3
    ph = 2 * np.pi * 0.25 / N
    out = []
    for lR2 in (0.15, 5e-2, 5e-3, 2.4e-3, 1e-3, 1e-4):
        lEO = 0.2 if lR2 > 0.1 else lR2 * (0.2 / 0.15)
        r = separation(N, lEO, 0.12, lR2, ph, gam)
        print("      %-10.2e  %.4f     %.4e   %.4e   %7.2f"
              % (lR2, r["peakT"], r["dT"], r["dPz"], r["dPz"] / lR2))
        out.append((lR2, r["peakT"], r["dT"], r["dPz"]))
    print()
    print("     The ratio tends to about 15.5 along lEO=(4/3)lR2.")
    print("     This is not the partial response to lR2 at fixed lEO.")
    print("     Independent-axis controls are in occupied_state_audit.py.")
    return out


if __name__ == "__main__":
    ok = validate()
    if not ok:
        sys.exit(1)
    sweep_geometry()
    sweep_path()
