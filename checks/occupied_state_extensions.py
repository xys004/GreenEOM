"""Independent numerical illustrations of the common occupied state in manuscript Sec. II."""
from occupied_state_audit import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

N=10; e=.039; r1=.12; phi=2*np.pi*.25/N
H0=famH(N,e,r1,0.,phi,0.,0); P,E=projector(H0,N)
V=famH(N,e,r1,1.,phi,0.,0)-H0

def range_term(m):
    H=np.zeros((2*N,2*N),complex)
    for n in range(N):
        for j,s in enumerate((1,-1)):
            i=2*n+j; k=2*((n+m)%N)+1-j
            value=-s/2*np.exp(-1j*s*(2*np.pi*n/N+m*np.pi/N))
            H[i,k]+=value; H[k,i]+=value.conjugate()
    return H

out={'annihilator_certificate':{'trace_PV':float(abs(np.trace(P@V))),
      'empty_occupied_block_norm':float(np.linalg.norm((np.eye(2*N)-P)@V@P))}}
out['ranges']=[]
for m in (2,3,4,6):
    W=range_term(m); Q,EE=projector(H0+.03*W,N)
    residual=max(np.linalg.norm(W@basis(N,k)-basis(N,k)@(-np.sin(m*k)*sy)) for k in np.pi*(2*np.arange(N)+1)/N)
    out['ranges'].append(dict(range=m,block_residual=float(residual),
        projector_change=float(np.linalg.norm(Q-P)),trace_PV=float(abs(np.trace(P@W))),
        empty_occupied_block_norm=float(np.linalg.norm((np.eye(2*N)-P)@W@P))))
    assert residual<1e-12
    if m%2==0: assert np.linalg.norm(Q-P)<1e-12

# A sufficient analytic interval obtained by retaining one fixed chemical potential.
mu=-e*np.sin(2*np.pi/N); lower=-np.inf; upper=np.inf
for k in np.pi*(2*np.arange(N)+1)/N:
    hb=block(N,k,e,r1,0.,phi,0.,0)
    ev=np.linalg.eigvalsh(hb)
    if ev[0]<mu<ev[1]: continue
    d0=float(np.trace(hb).real/2); dz=float((hb[0,0]-hb[1,1]).real/2)
    b=r1*np.sin(k-phi); v=np.sin(2*k)
    bound=np.sqrt((mu-d0)**2-dz**2)
    if abs(v)>1e-10:
        limits=sorted(((-bound-b)/v,(bound-b)/v))
        lower=max(lower,limits[0]); upper=min(upper,limits[1])
out['sufficient_fixed_mu_interval']={'mu':float(mu),'lower':float(lower),'upper':float(upper),
    'status':'Analytic inequality evaluated in float64, sufficient, not claimed to be the maximal canonical interval.'}

# No dynamics from changing only lambda_R2, even beyond the ground-state plateau.
# For a Slater determinant it suffices to propagate its occupied projector.
quench=[]
for r2 in (.3,2.):
    H=H0+r2*V; ev,U=np.linalg.eigh(H)
    Ut=(U*np.exp(-1j*ev*3.7))@U.conj().T
    Q,eg=projector(H,N)
    quench.append(dict(final_r2=r2,time=3.7,projector_evolution=float(np.linalg.norm(Ut@P@Ut.conj().T-P)),
        final_ground_projector_difference=float(np.linalg.norm(Q-P)),
        energy_above_final_ground=float(np.trace(P@H).real-eg[:N].sum())))
out['quenches']=quench

xs=np.linspace(0,.3,101)
pd=[]; ed=[]; ps=[]; es=[]
for x in xs:
    Q,ee=projector(H0+x*V,N)
    pd.append(np.linalg.norm(Q-P)); ed.append(np.max(np.abs(ee-E)))
    Hb=famH(N,e,r1,0.,phi,0.,1); Pb,Eb=projector(Hb,N)
    Qb,eeb=projector(famH(N,e,r1,x,phi,0.,1),N)
    ps.append(np.linalg.norm(Qb-Pb)); es.append(np.max(np.abs(eeb-Eb)))
fig,axs=plt.subplots(1,2,figsize=(9.2,3.6),layout='constrained')
axs[0].semilogy(xs,np.maximum(pd,1e-16),label='Phaseless',color='#B3341F')
axs[0].semilogy(xs,np.maximum(ps,1e-16),label='Phase restored',color='#2E7D57')
axs[0].set(ylabel=r'$\|P_{\rm occ}(\lambda)-P_{\rm occ}(0)\|_F$',xlabel=r'$\lambda_{R2}/t$',title='Occupied state')
axs[0].legend(loc='center right',frameon=False)
axs[1].plot(xs,ed,label='Phaseless',color='#B3341F')
axs[1].plot(xs,es,label='Phase restored',color='#2E7D57')
axs[1].set(ylabel=r'$\max_j|E_j(\lambda)-E_j(0)|/t$',xlabel=r'$\lambda_{R2}/t$',title='Single-particle spectrum')
axs[1].legend(loc='upper left',frameon=False)
for ax in axs: ax.grid(color='#dddddd',linewidth=.5)
for ext in ('pdf','eps'):
    fig.savefig(ROOT / 'figures' / ('Fig2.'+ext))
plt.close(fig)
(ROOT / 'data' / 'occupied_state_extensions.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
