"""Do the Hellmann-Feynman diagonal elements pair as +-x inside the occupied
block? If they do, the cancellation IS level-by-level, and the earlier reading
of the finite shifts was an artefact of re-sorting eigenvalues at each coupling."""
import numpy as np
from selection_rule import famH
from trace_rule import V

def pairing_defect(N, lEO, lR1, lR2, phi):
    ph = 2 * np.pi * phi / N
    w, U = np.linalg.eigh(famH(N, lEO, lR1, lR2, ph, 0.0, 0))
    d = np.sort(np.real(np.diag(U.conj().T @ V(N, ph, 0.0, 0) @ U))[:N])
    return np.abs(d + d[::-1]).max(), d.sum()

print("defecto de emparejamiento  max|d_i + d_(N+1-i)|  sobre el bloque OCUPADO")
print("(cero = los V_ii vienen en pares +-x exactos)")
print()
print("   N   N%4      lEO   lR1   lR2   phi      defecto      suma")
for N in (6, 8, 10, 12, 14, 18):
    for (lEO, lR1, lR2, phi) in ((0.039, 0.12, 0.15, 0.25),
                                 (0.2, 0.3, 0.3, 0.4),
                                 (0.0, 0.05, 0.05, 0.125)):
        df, s = pairing_defect(N, lEO, lR1, lR2, phi)
        flag = "  <-- pares exactos" if df < 1e-9 else ""
        print(f"  {N:2d}    {N%4}     {lEO:5.3f} {lR1:5.2f} {lR2:5.2f} {phi:5.3f}"
              f"   {df:9.2e}  {s:+9.2e}{flag}")
    print()

print("control: con la fase de Peierls restaurada, mismo test")
def pairing_defect_ph(N, lEO, lR1, lR2, phi):
    ph = 2 * np.pi * phi / N
    w, U = np.linalg.eigh(famH(N, lEO, lR1, lR2, ph, 0.0, 1))
    d = np.sort(np.real(np.diag(U.conj().T @ V(N, ph, 0.0, 1) @ U))[:N])
    return np.abs(d + d[::-1]).max(), d.sum()
for N in (6, 10, 14):
    df, s = pairing_defect_ph(N, 0.039, 0.12, 0.15, 0.25)
    print(f"  N={N:2d}  defecto = {df:.2e}   suma = {s:+.2e}")
