"""Compare the same isolated occupied state with finite-energy transport.

Run from the extracted supplement: python -B checks/common_state_transport.py
The contacts define a separate quadratic scattering problem; no invariance
of an open-system many-body state is asserted. NumPy and matplotlib required.
"""
from pathlib import Path
import contextlib, io, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()):
    from selection_rule import famH
    from transport_sweep import tmat

ROOT=Path(__file__).resolve().parents[1]
N=10; e=.039; r=.12; phi=2*np.pi*.25/N; gamma=.3; omega=.3
couplings=np.linspace(0,.3,61)
def occupied(H):
    E,U=np.linalg.eigh(H)
    return U[:,:N]@U[:,:N].conj().T,E
H0=famH(N,e,r,0.,phi,0.,0)
P0,E0=occupied(H0)
rows=[]
for lam in couplings:
    H=famH(N,e,r,lam,phi,0.,0)
    P,E=occupied(H)
    M=tmat(N,e,r,lam,phi,0,gamma,omega)
    # Independent two-by-two contact-block Caroli trace.
    GL=np.zeros_like(H); GR=np.zeros_like(H)
    GL[:2,:2]=gamma*np.eye(2)
    GR[N:N+2,N:N+2]=gamma*np.eye(2)
    G=np.linalg.inv((omega+1j*1e-8)*np.eye(2*N)-H+.5j*(GL+GR))
    T=float(np.trace(GR@G@GL@G.conj().T).real)
    assert abs(T-M.sum())<1e-12
    rows.append(dict(lam=float(lam),projector_residual=float(np.linalg.norm(P-P0)),
        energy_residual=float(abs(E[:N].sum()-E0[:N].sum())),gap=float(E[N]-E[N-1]),
        transmission=T,polarization=float((M[0].sum()-M[1].sum())/T)))
assert max(q['projector_residual'] for q in rows)<2e-14
assert min(q['gap'] for q in rows)>1e-8
assert abs(rows[-1]['transmission']-rows[0]['transmission'])>1e-5
out=dict(parameters=dict(N=N,Ne=N,e=e,r=r,flux_fraction=.25,rho=0,gamma=gamma,omega=omega,eta=1e-8),
    status='Numerical illustration of the analytically proved common state, not its proof.',
    max_projector_residual=max(q['projector_residual'] for q in rows),
    max_energy_residual=max(q['energy_residual'] for q in rows),
    min_fermi_gap=min(q['gap'] for q in rows),rows=rows)
(ROOT/'data/common_state_transport.json').write_text(json.dumps(out,indent=2)+'\n')
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(7.1,2.55),layout='constrained')
axs[0].plot(couplings,[max(q['projector_residual'],1e-16) for q in rows],color='#25624f')
axs[0].set_yscale('log'); axs[0].set_ylim(1e-16,1e-13)
axs[0].set_ylabel(r'$\|P(\lambda)-P(0)\|_{\rm F}$')
axs[0].set_title('(a) Isolated occupied state',fontsize=10)
axs[1].plot(couplings,[q['transmission'] for q in rows],color='#b54732')
axs[1].set_ylabel(r'$T(\omega=0.3t)$')
axs[1].set_title('(b) With two contacts',fontsize=10)
for ax in axs:
    ax.set_xlabel(r'$\lambda/t$');ax.grid(alpha=.15)
for ext in ('pdf','eps'):
    fig.savefig(ROOT/'figures'/('Fig_common_transport.'+ext))
plt.close(fig)
print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
print('Endpoints:',json.dumps([rows[0],rows[-1]]))
print('PASS: common isolated state and changing transport at identical parameters.')
