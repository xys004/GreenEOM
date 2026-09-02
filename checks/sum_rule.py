"""Is the whole spectrum blind to the phaseless lambda_R2, or only the
occupied sum? The two have different mechanisms and different consequences."""
import numpy as np
from selection_rule import famH

print("phaseless lambda_R2: what exactly is invariant?")
print()
for N in (6, 8, 10, 14):
    ph = 2 * np.pi * 0.25 / N
    e0 = np.sort(np.linalg.eigvalsh(famH(N, 0.039, 0.12, 0.0,  ph, 0.0, 0)).real)
    ex = np.sort(np.linalg.eigvalsh(famH(N, 0.039, 0.12, 0.30, ph, 0.0, 0)).real)
    full = np.abs(ex - e0).max()
    occ  = abs(ex[:N].sum() - e0[:N].sum())
    # where does the movement sit?
    moved = np.abs(ex - e0)
    lo, hi = moved[:N].max(), moved[N:].max()
    print(f"  N={N:2d} (N%4={N%4})  max|dE_i| todo el espectro = {full:.3e}")
    print(f"          max|dE_i| ocupados = {lo:.3e}   vacios = {hi:.3e}")
    print(f"          |suma ocupada movida| = {occ:.3e}")
    if full > 1e-10 and occ < 1e-10:
        s = np.sort(ex[:N] - e0[:N])
        print(f"          -> los niveles SI se mueven y la suma NO: cancelacion por pares")
        print(f"             desplazamientos ocupados: {np.array2string(s, precision=4, max_line_width=100)}")
    elif full < 1e-10:
        print(f"          -> el espectro entero es ciego al termino")
    print()
