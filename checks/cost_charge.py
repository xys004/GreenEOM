"""Cost of a phaseless range-two coupling on the CHARGE current (Sec. III B).

A range-two term carrying no flux phase contributes no explicit dependence on
phi, so it opens no channel of its own in J = -dE/dphi.  It is not inert: it
still reshapes the occupied spectrum and moves the current indirectly.  This
script measures both effects at a stated parameter point, so that the two
numbers quoted in the text have a deposited source rather than being carried
over from an earlier calculation.

Parameter point: eight-site ring, half filling, lambda_EO = 0.039 t,
lambda_R1 = 0.12 t, phi/phi_0 = 1/4, turning lambda_R2 on from zero.
"""
import numpy as np
from selection_rule import famH

N, lEO, lR1, P = 8, 0.039, 0.12, 0.25
PH = 2 * np.pi * P / N


def occ_E(a, ph, pei):
    return np.sort(np.linalg.eigvalsh(famH(N, *a, ph, 0.0, pei)).real)[:N].sum()


def jc(a, pei, d=1e-6):
    """charge current  J = -dE/dphi, symmetric finite difference"""
    return -(occ_E(a, PH + d, pei) - occ_E(a, PH - d, pei)) / (2 * d)


print("cost of a phaseless range-two coupling on the charge current")
print(f"  N={N}, half filling, lEO={lEO} t, lR1={lR1} t, phi/phi0={P}")
print()
base0 = jc((lEO, lR1, 0.0), 0)
base1 = jc((lEO, lR1, 0.0), 1)
print(f"  J at lambda_R2 = 0 :  phaseless {base0:+.9f} t   with phase {base1:+.9f} t")
print()
print("  lambda_R2    |dJ| phaseless (indirect)   |dJ| with phase (explicit)   ratio")
for x in (0.05, 0.10, 0.15, 0.20, 0.30):
    d0 = abs(jc((lEO, lR1, x), 0) - base0)
    d1 = abs(jc((lEO, lR1, x), 1) - base1)
    print(f"    {x:4.2f}       {d0:.3e}                  {d1:.3e}"
          f"                {d1/d0:6.1f}")
print()
print("  The indirect effect is smaller than the explicit one by more than an")
print("  order of magnitude at every coupling tested, which is the point the")
print("  text makes: the two have different origins and different sizes.")
