"""The selection rule, derived.

  (1) Tgen = T_1 (x) exp(-i sigma_z pi/N) commutes with H and with
      V = dH/d(lambda_R2) when the range-two terms carry no flux phase.
  (2) Tgen^N = -1, so the allowed momenta are k = pi(2l+1)/N: HALF-INTEGER.
  (3) k = pi/2 is in that set  <=>  N = 2(2l+1)  <=>  N = 2 (mod 4).
  (4) V vanishes identically on the k = +-pi/2 blocks.
  (5) Tr_k V = 0 in every block, so complete blocks contribute +- pairs.
  (6) At half filling exactly the two k = +-pi/2 blocks are split by the Fermi
      level. When they exist (N = 2 mod 4) they contribute nothing because V is
      zero there, and everything else cancels in pairs -> dE/d(lambda_R2) = 0.
      When they do not (N = 0 mod 4, N odd) the split blocks carry V =/= 0 and
      the sum survives.
  (7) Restoring the flux phase does NOT break (1) -- Tgen still commutes with
      both. What it breaks is (4): V is no longer zero on the k = +-pi/2 blocks.
  (8) The rule therefore holds exactly when the ONLY blocks the Fermi level
      splits are k = +-pi/2. At weak coupling that is the case; at strong
      coupling the bands overlap, further blocks straddle E_F, and V does not
      vanish on those.
"""
import numpy as np
from selection_rule import famH
from trace_rule import V

def Tgen(N):
    T = np.zeros((N, N))
    for n in range(N): T[(n+1) % N, n] = 1.0
    return np.kron(T, np.diag([np.exp(-1j*np.pi/N), np.exp(1j*np.pi/N)]))

ok = lambda b: "OK " if b else "FAIL"
P = []
SIZES = [5,6,7,8,9,10,11,12,13,14,17,18]
PARS = [(0.039,0.12,0.15,0.25),(0.2,0.3,0.3,0.4),(0.0,0.05,0.05,0.125)]

# (2) Tgen^N = -1
P.append(("Tgen^N = -1", all(
    np.linalg.norm(np.linalg.matrix_power(Tgen(N), N) + np.eye(2*N)) < 1e-9
    for N in SIZES)))

# (3) arithmetic: pi/2 allowed iff N = 2 mod 4
def has_half_pi(N):
    ks = [np.pi*(2*l+1)/N for l in range(N)]
    return any(abs(k - np.pi/2) < 1e-12 for k in ks)
P.append(("k=pi/2 permitido  <=>  N = 2 (mod 4)",
          all(has_half_pi(N) == (N % 4 == 2) for N in SIZES)))

# (1),(4),(5),(7)
c1 = c4 = c5 = c7 = True
for N in SIZES:
    for (lEO,lR1,lR2,phi) in PARS:
        ph = 2*np.pi*phi/N
        H = famH(N,lEO,lR1,lR2,ph,0.0,0); Vm = V(N,ph,0.0,0); G = Tgen(N)
        c1 &= np.linalg.norm(G@H-H@G)/np.linalg.norm(H) < 1e-10
        c1 &= np.linalg.norm(G@Vm-Vm@G)/max(np.linalg.norm(Vm),1e-30) < 1e-10
        ev,W = np.linalg.eig(G); o=np.argsort(np.angle(ev)); ev,W = ev[o],W[:,o]
        a=np.angle(ev); i=0
        while i < len(a):
            j=i
            while j+1<len(a) and abs(a[j+1]-a[i])<1e-8: j+=1
            B=np.linalg.qr(W[:,i:j+1])[0]; vb=B.conj().T@Vm@B
            c5 &= abs(np.trace(vb).real) < 1e-9
            if abs(abs(a[i])-np.pi/2) < 1e-8:
                c4 &= np.linalg.norm(vb) < 1e-9
            i=j+1
        Vp = V(N,ph,0.0,1)
        if np.linalg.norm(Vp) > 1e-12 and N % 4 == 2:
            ev2,W2 = np.linalg.eig(G); o2=np.argsort(np.angle(ev2))
            a2=np.angle(ev2[o2]); W2=W2[:,o2]; i2=0
            while i2 < len(a2):
                j2=i2
                while j2+1<len(a2) and abs(a2[j2+1]-a2[i2])<1e-8: j2+=1
                if abs(abs(a2[i2])-np.pi/2) < 1e-8:
                    B2=np.linalg.qr(W2[:,i2:j2+1])[0]
                    c7 &= np.linalg.norm(B2.conj().T@Vp@B2) > 1e-3
                i2=j2+1
P += [("[Tgen,H]=[Tgen,V]=0 (con y sin la fase)", c1),
      ("V = 0 identicamente en los bloques k=+-pi/2", c4),
      ("Tr_k V = 0 en TODOS los bloques", c5),
      ("con la fase restaurada V ya NO se anula en k=+-pi/2", c7)]

# (6) the payoff
c6 = True
for N in SIZES:
    for (lEO,lR1,lR2,phi) in PARS:
        ph = 2*np.pi*phi/N
        G = Tgen(N)
        w,U = np.linalg.eigh(famH(N,lEO,lR1,lR2,ph,0.0,0))
        s = np.real(np.diag(U.conj().T@V(N,ph,0.0,0)@U))[:N].sum()
        gap = w[N]-w[N-1]
        if gap > 1e-2:                       # Fermi level well defined
            ev3,W3 = np.linalg.eig(G); o3=np.argsort(np.angle(ev3))
            a3=np.angle(ev3[o3]); W3=W3[:,o3]
            wf=np.sort(np.linalg.eigvalsh(famH(N,lEO,lR1,lR2,ph,0.0,0)).real)
            Ef=(wf[N-1]+wf[N])/2; only_hp=True; i3=0
            while i3 < len(a3):
                j3=i3
                while j3+1<len(a3) and abs(a3[j3+1]-a3[i3])<1e-8: j3+=1
                B3=np.linalg.qr(W3[:,i3:j3+1])[0]
                eb3=np.linalg.eigvalsh(B3.conj().T@famH(N,lEO,lR1,lR2,ph,0.0,0)@B3).real
                if 0 < int((eb3<Ef).sum()) < len(eb3):
                    if abs(abs(a3[i3])-np.pi/2) > 1e-8: only_hp=False
                i3=j3+1
            c6 &= (abs(s) < 1e-9) == (N % 4 == 2 and only_hp)
P.append(("dE/d(lambda_R2)=0  <=>  N=2(mod4) Y solo k=+-pi/2 partidos", c6))

for name, v in P: print(f"  {ok(v)}  {name}")
print()
print("VERDICT:", "PASS" if all(v for _,v in P) else "FAIL")
