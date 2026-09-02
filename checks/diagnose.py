import numpy as np
from selection_rule import famH
from trace_rule import V
from derivation import Tgen

SIZES=[5,6,7,8,9,10,11,12,13,14,17,18]
PARS=[(0.039,0.12,0.15,0.25),(0.2,0.3,0.3,0.4),(0.0,0.05,0.05,0.125)]

print("A) con la fase restaurada, que pasa con Tgen y con V en k=+-pi/2?")
for N in (6,10,14):
    ph=2*np.pi*0.25/N
    Vp=V(N,ph,0.0,1); G=Tgen(N)
    r=np.linalg.norm(G@Vp-Vp@G)/np.linalg.norm(Vp)
    Hp=famH(N,0.039,0.12,0.15,ph,0.0,1)
    rh=np.linalg.norm(G@Hp-Hp@G)/np.linalg.norm(Hp)
    ev,W=np.linalg.eig(G); o=np.argsort(np.angle(ev)); ev,W=ev[o],W[:,o]
    a=np.angle(ev); nz=None; i=0
    while i<len(a):
        j=i
        while j+1<len(a) and abs(a[j+1]-a[i])<1e-8: j+=1
        if abs(abs(a[i])-np.pi/2)<1e-8:
            B=np.linalg.qr(W[:,i:j+1])[0]; nz=np.linalg.norm(B.conj().T@Vp@B)
        i=j+1
    print(f"   N={N:2d}  ||[T,V_fase]||={r:.1e}  ||[T,H_fase]||={rh:.1e}  ||V|| en k=pi/2 = {nz:.2e}")

print()
print("B) donde falla el bicondicional")
for N in SIZES:
    for p in PARS:
        lEO,lR1,lR2,phi=p; ph=2*np.pi*phi/N
        w,U=np.linalg.eigh(famH(N,lEO,lR1,lR2,ph,0.0,0))
        s=np.real(np.diag(U.conj().T@V(N,ph,0.0,0)@U))[:N].sum()
        gap=w[N]-w[N-1]
        if gap>1e-2 and (abs(s)<1e-9)!=(N%4==2):
            print(f"   N={N:2d} (N%4={N%4}) lEO={lEO} lR1={lR1} lR2={lR2} phi={phi}"
                  f"   suma={s:+.3e}  gap={gap:.3e}")
