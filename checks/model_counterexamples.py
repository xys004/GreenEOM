"""Focused counterexamples to overbroad statements in the manuscript."""
from occupied_state_audit import *


def customH(N,e,r1,r2,phi,conv='canonical', geometric_probe=0.):
    H=np.zeros((2*N,2*N),complex)
    def add(i,j,v): H[i,j]+=v; H[j,i]+=v.conjugate()
    def th(n,m):
        if conv=='site': return 2*np.pi*n/N
        if conv=='average': return np.pi*(n+(n+m)%N)/N
        return 2*np.pi*n/N+m*np.pi/N
    for n in range(N):
        for j,s in enumerate((1,-1)):
            ix=2*n+j; jx=2*((n+1)%N)+j; kx=2*((n+2)%N)+j
            p1=np.exp(1j*(phi+geometric_probe)); p2=np.exp(2j*geometric_probe)
            add(ix,jx,p1)
            add(ix,jx+1-2*j,-r1/2*s*np.exp(-1j*s*th(n,1))*p1)
            add(ix,kx,1j*e*s/2*p2)
            add(ix,kx+1-2*j,-r2/2*s*np.exp(-1j*s*th(n,2))*p2)
    return H


out={}
N=10; phi=2*np.pi*.25/N
spectra={c:np.linalg.eigvalsh(customH(N,0.,.3,0.,phi,c)) for c in ('canonical','site','average')}
out['range_one_angle_spectral_difference']={c:float(np.max(np.abs(spectra[c]-spectra['canonical']))) for c in ('site','average')}

H=famH(N,.039,.12,.15,phi,0.,0); P,E=projector(H,N)
Sz=np.kron(np.eye(N),sz)
d=1e-6
Oflux=(famH(N,.039,.12,.15,phi+d,0.,0)-famH(N,.039,.12,.15,phi-d,0.,0))/(2*d)
Ogeom=(customH(N,.039,.12,.15,phi,geometric_probe=d)-customH(N,.039,.12,.15,phi,geometric_probe=-d))/(2*d)
out['spin_current_raw_convention']={
    'flux_anticommutator':float(-np.trace(P@(Sz@Oflux+Oflux@Sz)/2).real),
    'all_bonds_anticommutator':float(-np.trace(P@(Sz@Ogeom+Ogeom@Sz)/2).real),
    'operator_difference_norm':float(np.linalg.norm(Ogeom-Oflux)),
    'note':'Both raw quantities use sigma_z, with no physical 1/(2N) factor. The extra probe weights EVERY directed bond by its range, even when the AB phase is omitted.'}

out['spin_current_other_points']=[]
for N,e,r1,r2,p,nf in ((8,.039,.12,.15,.25,8),(10,.2,.3,.3,.4,10),(10,.039,.12,.15,.25,8)):
    phi=2*np.pi*p/N; Sz=np.kron(np.eye(N),sz)
    H=famH(N,e,r1,r2,phi,0.,0); P,E=projector(H,nf)
    Of=(famH(N,e,r1,r2,phi+d,0.,0)-famH(N,e,r1,r2,phi-d,0.,0))/(2*d)
    Og=(customH(N,e,r1,r2,phi,geometric_probe=d)-customH(N,e,r1,r2,phi,geometric_probe=-d))/(2*d)
    out['spin_current_other_points'].append(dict(N=N,nf=nf,lEO=e,lR1=r1,lR2=r2,flux_fraction=p,
        flux=float(-np.trace(P@(Sz@Of+Of@Sz)/2).real),
        all_bonds=float(-np.trace(P@(Sz@Og+Og@Sz)/2).real)))

# A large gauge transformation closes on the ring only if N*chi=2*pi*integer.
def gauge_residual(chi,pei):
    U=np.diag(np.repeat(np.exp(1j*chi*np.arange(N)),2))
    A=famH(N,.039,.12,.15,phi,0.,pei)
    B=famH(N,.039,.12,.15,phi+chi,0.,pei)
    return float(np.linalg.norm(U.conj().T@A@U-B))
N=10; phi=2*np.pi*.25/N
out['ring_gauge_transform']={
    'arbitrary_chi_correct_model':gauge_residual(.071,1),
    'quantized_chi_correct_model':gauge_residual(2*np.pi/N,1),
    'quantized_chi_phaseless_model':gauge_residual(2*np.pi/N,0)}

# PT counterexample verified as a matrix, not merely a characteristic polynomial.
Hpt=sp.Matrix([[2*sp.I,1],[1,-2*sp.I]]); parity=sp.Matrix([[0,1],[1,0]])
out['PT_counterexample']={'PT_residual':str(parity*sp.conjugate(Hpt)*parity-Hpt),
                           'eigenvalues':[str(v) for v in Hpt.eigenvals()]}
out['quadratic_closure_counterexample']='H=epsilon*c1dag*c1 on an N-site spinful ring: seed G(c1,c1dag) closes on 1 unknown; H is quadratic and periodic.'
out['support_counterexample']='Changing H_12=t to exp(i*phi)*t with t nonzero changes a coefficient but not the support {G_2j}.'
(ROOT / 'data' / 'model_counterexamples.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
