"""David's point 16: turn the numerical particle-hole result into an algebraic
statement. Gamma = diag((-1)^n) (x) 1_spin. Claim: {Gamma, H_range1} = 0 exactly
(so the range-1 model is particle-hole symmetric), while H_range2 COMMUTES with
Gamma instead of anticommuting -- so adding it breaks the anticommutation."""
import numpy as np
from selection_rule import famH

def Gamma(N):
    g = np.array([(-1.0)**n for n in range(N)])
    return np.kron(np.diag(g), np.eye(2))

def parts(N, ph):
    """H split into its range-1 and range-2 pieces (phaseless range 2)."""
    H_all = famH(N, 0.039, 0.12, 0.15, ph, 0.0, 0)      # lEO,lR1,lR2 all on
    H_r1  = famH(N, 0.0,   0.12, 0.0,  ph, 0.0, 0)       # only range-1 Rashba+hop
    H_r2  = H_all - H_r1                                  # the range-two remainder
    return H_all, H_r1, H_r2

for N in (6, 8, 10):
    ph = 2*np.pi*0.25/N
    G = Gamma(N)
    H, H1, H2 = parts(N, ph)
    a1 = np.linalg.norm(G@H1 + H1@G)     # {G,H1}
    c2 = np.linalg.norm(G@H2 - H2@G)     # [G,H2]
    a2 = np.linalg.norm(G@H2 + H2@G)     # {G,H2}
    aa = np.linalg.norm(G@H  + H@G)      # {G,H} full
    print(f"N={N:2d}:  ||{{G,H_range1}}|| = {a1:.1e}   (should be 0: range-1 anticommutes)")
    print(f"        ||[G,H_range2]|| = {c2:.1e}   ||{{G,H_range2}}|| = {a2:.1e}   "
          f"(range-2 COMMUTES, does not anticommute)")
    print(f"        ||{{G,H_full}}||   = {aa:.1e}   (nonzero: PH broken by range-2)")
    print()
