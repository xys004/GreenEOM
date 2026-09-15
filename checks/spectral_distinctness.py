"""How distinct are the spectra whose occupied state is common?

The title of the paper says "spectrally distinct rings". The text asserts six
times that the excitation spectrum changes while the occupied Slater state does
not, and nowhere puts a number on it. This measures the three levels at the
same parameter point used in the transport comparison:

    level 1   the occupied projector P          -- invariant, by the theorem
    level 2   the single-particle spectrum      -- NOT invariant, measured here
    level 3   the finite-energy transmission    -- NOT invariant, Sec. V

Parameters follow the transport figure: N = N_e = 10, e/t = 0.039, r/t = 0.12,
Phi/Phi_0 = 1/4, rho = 0 (phaseless range-two hop), lambda from 0 to 0.3 t.
"""
import os
import sys
import io
import contextlib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):
    from selection_rule import famH

N, NE = 10, 10
E_INTR, R1 = 0.039, 0.12
PH = 2 * np.pi * 0.25 / N           # Phi/Phi_0 = 1/4
PEI = 0                             # rho = 0, phaseless range-two hop


def levels_and_projector(lam):
    H = famH(N, E_INTR, R1, lam, PH, 0.0, PEI)
    w, U = np.linalg.eigh(H)
    order = np.argsort(w.real)
    w, U = w.real[order], U[:, order]
    P = U[:, :NE] @ U[:, :NE].conj().T
    return w, P


def main():
    lams = np.linspace(0.0, 0.3, 61)
    w0, P0 = levels_and_projector(0.0)

    proj, eocc, shift_occ, shift_emp, shift_exc, gaps = [], [], [], [], [], []
    for lam in lams:
        w, P = levels_and_projector(lam)
        proj.append(np.linalg.norm(P - P0))
        eocc.append(abs(w[:NE].sum() - w0[:NE].sum()))
        shift_occ.append(np.abs(w[:NE] - w0[:NE]).max())
        shift_emp.append(np.abs(w[NE:] - w0[NE:]).max())
        # excitation energies, every empty level above every occupied one
        exc = w[NE:][None, :] - w[:NE][:, None]
        exc0 = w0[NE:][None, :] - w0[:NE][:, None]
        shift_exc.append(np.abs(exc - exc0).max())
        gaps.append(w[NE] - w[NE - 1])

    print("=" * 74)
    print(" Three levels at one parameter point, N = N_e = 10, phaseless")
    print(" e/t = 0.039, r/t = 0.12, Phi/Phi_0 = 1/4, lambda in [0, 0.3] t")
    print("=" * 74)
    print()
    print(" level 1   the occupied state is common")
    print("   max_lambda || P(lambda) - P(0) ||_F      = %.3e" % max(proj))
    print("   max_lambda | E_occ(lambda) - E_occ(0) |  = %.3e t" % max(eocc))
    print("   smallest Fermi gap on the grid           = %.7f t" % min(gaps))
    print()
    print(" level 2   the spectrum is not")
    print("   max shift of an occupied level           = %.6f t" % max(shift_occ))
    print("   max shift of an empty level              = %.6f t" % max(shift_emp))
    print("   max shift of an excitation energy        = %.6f t" % max(shift_exc))
    print()
    ratio = max(shift_exc) / max(proj) if max(proj) > 0 else float("inf")
    print("   The excitation spectrum moves by %.3f t while the occupied" % max(shift_exc))
    print("   projector moves by %.1e. The two differ by %.1e in magnitude."
          % (max(proj), ratio))
    print()

    # the individual occupied levels move, yet their sum does not
    w_end, _ = levels_and_projector(0.3)
    print(" why both can be true at once")
    print("   occupied levels at lambda = 0     :",
          " ".join("%+.4f" % x for x in w0[:NE]))
    print("   occupied levels at lambda = 0.3 t :",
          " ".join("%+.4f" % x for x in w_end[:NE]))
    print("   their sums                        : %+.10f  vs  %+.10f"
          % (w0[:NE].sum(), w_end[:NE].sum()))
    print("   Individual levels move; the sum, the projector and hence every")
    print("   equal-time observable of the occupied state do not.")


if __name__ == "__main__":
    main()
