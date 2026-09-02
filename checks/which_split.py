import numpy as np
from selection_rule import famH
from trace_rule import V
from derivation import Tgen

def split_ks(N,lEO,lR1,lR2,phi,pei=0):
    ph=2*np.pi*phi/N
    H=famH(N,lEO,lR1,lR2,ph,0.0,pei); G=Tgen(N)
    ev,W=np.linalg.eig(G); o=np.argsort(np.angle(ev)); ev,W=ev[o],W[:,o]
    a=np.angle(ev); w=np.sort(np.linalg.eigvalsh(H).real); Ef=(w[N-1]+w[N])/2
    out=[]; i=0
    while i<len(a):
        j=i
        while j+1<len(a) and abs(a[j+1]-a[i])<1e-8: j+=1
        B=np.linalg.qr(W[:,i:j+1])[0]; eb=np.linalg.eigvalsh(B.conj().T@H@B).real
        if 0<int((eb<Ef).sum())<len(eb): out.append(a[i]/np.pi)
        i=j+1
    return out

print("que bloques parte el nivel de Fermi (k/pi)")
print()
print("  acoplo debil (lEO=0.039, lR1=0.12, lR2=0.15, phi=0.25):")
for N in (6,10,14,18):
    print(f"     N={N:2d}: {[round(x,3) for x in split_ks(N,0.039,0.12,0.15,0.25)]}")
print()
print("  acoplo fuerte (lEO=0.2, lR1=0.3, lR2=0.3, phi=0.4)  <- donde la regla falla:")
for N in (6,10,14,18):
    print(f"     N={N:2d}: {[round(x,3) for x in split_ks(N,0.2,0.3,0.3,0.4)]}")
