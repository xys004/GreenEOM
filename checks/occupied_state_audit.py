"""Reproducible audit. All assertions are about the explicitly encoded 1D model.

Run: python -B checks/occupied_state_audit.py [path-to-GreenEOM]
Results are written to data/occupied_state_audit.json.
"""
from pathlib import Path
import contextlib
import io
import json
import sys
import numpy as np
import sympy as sp
import z3

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'checks'))
with contextlib.redirect_stdout(io.StringIO()):
    from selection_rule import famH
    from transport_sweep import tmat, validate

sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]], complex)
sz = np.diag([1., -1.])
eye = np.eye(2)
results = {'versions': {'numpy': np.__version__, 'sympy': sp.__version__, 'z3': z3.get_version_string()}}


def block(N, k, e, r1, r2, phi, psi=0., pei=0):
    a = np.pi / N
    q = k - phi
    q2 = k - pei * phi
    b2 = a - pei * psi
    d0 = 2*np.cos(q)*np.cos(a-psi) + e*np.cos(2*q2)*np.sin(2*b2)
    dz = -2*np.sin(q)*np.sin(a-psi) + e*np.sin(2*q2)*np.cos(2*b2)
    off = 1j*(r1*np.sin(q)*np.exp(1j*psi) + r2*np.sin(2*q2)*np.exp(2j*pei*psi))
    return np.array([[d0+dz, off], [off.conjugate(), d0-dz]])


def basis(N, k):
    B = np.zeros((2*N, 2), complex)
    for n in range(N):
        for j, s in enumerate((1, -1)):
            B[2*n+j, j] = np.exp(-1j*(k+s*np.pi/N)*n)/np.sqrt(N)
    return B


def projector(H, nf):
    E, U = np.linalg.eigh(H)
    return U[:, :nf] @ U[:, :nf].conj().T, E


def symbolic():
    k, a, phi, psi, e, r1, r2 = sp.symbols('k a phi psi e r1 r2', real=True)
    q = k-phi
    checks = {}
    # Derive forward plus backward bond contributions independently of block().
    for m in (1, 2, 4):
        lam, chi, eta = sp.symbols('lam chi eta', real=True)
        forward = -lam/2*sp.exp(sp.I*m*(chi-k))*sp.exp(sp.I*m*eta)
        backward = lam/2*sp.exp(-sp.I*m*(chi-k))*sp.exp(sp.I*m*eta)
        checks[f'rashba_bond_range_{m}'] = sp.simplify(sp.expand_complex(
            forward+backward-sp.I*lam*sp.sin(m*(k-chi))*sp.exp(sp.I*m*eta))) == 0
    for pei in (0, 1):
        f = r1*sp.sin(q)*sp.exp(sp.I*psi) + r2*sp.sin(2*(k-pei*phi))*sp.exp(2*sp.I*pei*psi)
        h = sp.Matrix([[0, sp.I*f], [-sp.I*sp.conjugate(f), 0]])
        dh = sp.simplify(sp.diff(h, r2))
        norm2 = sp.simplify(sp.trace(dh.conjugate().T*dh))
        checks[f'frobenius_norm_squared_pei_{pei}'] = sp.simplify(norm2-2*sp.sin(2*(k-pei*phi))**2) == 0
        # At psi=0 the spin-flip flux derivative is along sigma_x;
        # the Hamiltonian has only sigma_y and sigma_z components.
        checks[f'spin_flip_orthogonality_pei_{pei}'] = sp.simplify(
            sp.trace(h.subs(psi, 0)*sp.diff(h, psi).subs(psi, 0))) == 0
        if pei == 0:
            checks['zero_V_at_pi2_for_all_spin_flux'] = sp.simplify(dh.subs(k, sp.pi/2)) == sp.zeros(2)
    # Algebraic projector for nondegenerate H=dI+hx*sx+hy*sy+hz*sz.
    x, y, z, rad = sp.symbols('x y z rad', real=True, nonzero=True)
    M = sp.Matrix([[z, x-sp.I*y], [x+sp.I*y, -z]])
    P = (sp.eye(2)-M/rad)/2
    residual = sp.simplify(P*P-P)
    checks['spectral_projector_idempotency'] = all(
        sp.cancel(v).subs(rad**2, x*x+y*y+z*z).simplify() == 0 for v in residual)
    results['symbolic'] = checks
    assert all(checks.values()), checks


def verify_blocks():
    worst = 0.
    for N in (5, 6, 8, 10, 14):
        for pei in (0, 1):
            for psi in (0., .013):
                phi = 2*np.pi*.25/N
                H = famH(N, .039, .12, .15, phi, psi, pei)
                for l in range(N):
                    k = np.pi*(2*l+1)/N
                    B = basis(N, k)
                    hb = block(N, k, .039, .12, .15, phi, psi, pei)
                    worst = max(worst, float(np.linalg.norm(H@B-B@hb)))
    results['block_reduction_max_residual'] = worst
    assert worst < 1e-12


def projectors_and_counterexamples():
    rows = []
    for N, p, e, r1 in ((6,.25,.039,.12), (10,.25,.039,.12), (14,.25,.039,.12),
                         (10,.4,.2,.3), (8,.25,.039,.12), (8,0.,.039,.12)):
        for pei in (0, 1):
            phi = 2*np.pi*p/N
            P0, E0 = projector(famH(N,e,r1,0.,phi,0.,pei), N)
            maxP, maxE, gapmin, spectral = 0., 0., 1e9, 0.
            split = set()
            for r2 in np.linspace(0,.3,31):
                H = famH(N,e,r1,r2,phi,0.,pei)
                P, E = projector(H,N)
                mu = (E[N-1]+E[N])/2
                gapmin = min(gapmin,float(E[N]-E[N-1]))
                maxP = max(maxP,float(np.linalg.norm(P-P0)))
                maxE = max(maxE,abs(float(E[:N].sum()-E0[:N].sum())))
                spectral = max(spectral,float(np.max(np.abs(E-E0))))
                for l in range(N):
                    k = np.pi*(2*l+1)/N
                    eb = np.linalg.eigvalsh(block(N,k,e,r1,r2,phi,0.,pei))
                    if eb[0] < mu < eb[1]: split.add(round(k/np.pi,6))
            rows.append(dict(N=N,flux=p,lEO=e,lR1=r1,pei=pei,max_projector_change=maxP,
                             max_ground_energy_change=maxE,min_gap=gapmin,max_spectral_change=spectral,
                             split_k_over_pi=sorted(split)))
            if N in (6,10,14) and p==.25 and pei==0:
                assert maxP < 1e-12 and maxE < 1e-12 and spectral > 1e-3
    results['projectors'] = rows
    # The gap is genuinely finite: an interval, not an ambiguous occupation.
    # Finite T makes fully occupied blocks imperfectly occupied and loses the plateau.
    thermal = []
    for beta in (100.,20.,5.):
        mats=[]
        for r2 in (0.,.3):
            E,U=np.linalg.eigh(famH(10,.039,.12,r2,2*np.pi*.25/10,0.,0))
            lo,hi=E.min()-10,E.max()+10
            for _ in range(100):
                mu=(lo+hi)/2
                f=1/(1+np.exp(np.clip(beta*(E-mu),-700,700)))
                if f.sum()>10: hi=mu
                else: lo=mu
            mats.append((U*f)@U.conj().T)
        thermal.append({'beta_t':beta,'density_matrix_change':float(np.linalg.norm(mats[1]-mats[0]))})
    results['finite_temperature'] = thermal


def z3_checks():
    z = {}
    # Fixed occupations: trace-zero in full blocks, zero perturbation in split ones.
    for n in (0,1,2):
        v0,v1=z3.Reals('v0 v1')
        S=z3.Solver()
        S.add(v0+v1==0)
        if n==1: S.add(v0==0,v1==0)
        contribution = 0 if n==0 else (v0 if n==1 else v0+v1)
        S.add(contribution!=0)
        z[f'occupation_{n}']=str(S.check())
        assert S.check()==z3.unsat
    # H_PT = [[i*g,1],[1,-i*g]] has E^2=1-g^2.
    g,x=z3.Reals('g x')
    S=z3.Solver(); S.add(g==2,x*x==1-g*g)
    z['PT_does_not_imply_real_spectrum']=str(S.check())
    assert S.check()==z3.unsat
    # Invalid inference E(lambda,0) constant => d_psi E independent of lambda.
    lam,psi=sp.symbols('lam psi',real=True)
    E=lam*psi
    z['mixed_derivative_counterexample']=str(sp.diff(E,psi))
    results['z3'] = z


def transport_audit():
    capture=io.StringIO()
    with contextlib.redirect_stdout(capture): ok=validate()
    assert ok
    results['transport_validation']=capture.getvalue()
    N=6; phi=2*np.pi*.25/N
    # Vectorized energy sweep using exactly the same Hamiltonian/contact convention.
    ws=np.linspace(-2,2,12801)
    def metrics(e,r2):
        TT=[]; PP=[]
        for pei in (0,1):
            H=famH(N,e,.12,r2,phi,0.,pei)
            sig=np.zeros(2*N,complex); sig[[0,1,N,N+1]]=-.15j
            A=(ws[:,None,None]+1e-8j)*np.eye(2*N)-H-np.diag(sig)
            G=np.linalg.inv(A)
            M=.3**2*np.abs(G[:,N:N+2,:2])**2
            T=M.sum(axis=(1,2)); P=(M[:,0,:].sum(axis=1)-M[:,1,:].sum(axis=1))/T
            TT.append(T); PP.append(P)
        threshold=.1*max(T.max() for T in TT)
        mask=(TT[0]>=threshold)&(TT[1]>=threshold)
        iT=np.argmax(np.abs(TT[0]-TT[1])); ids=np.flatnonzero(mask)
        iP=ids[np.argmax(np.abs(PP[0]-PP[1])[mask])]
        return dict(lEO=e,lR2=r2,max_dT=float(abs(TT[0][iT]-TT[1][iT])),
                    max_dP=float(abs(PP[0][iP]-PP[1][iP])),omega_dT=float(ws[iT]),omega_dP=float(ws[iP]),
                    threshold=float(threshold))
    results['transport_paths']=[metrics(.2/.15*r,r) for r in (.0024,.001,.0001)]
    # Vary intrinsic SO and Rashba independently. This is NOT a material calibration.
    results['transport_axis_controls']=[metrics(e,r) for e,r in
        ((.0032,0.),(0.,.0024),(.0032,.0024),(.0039/1.6,.0007/1.6))]


if __name__=='__main__':
    symbolic(); verify_blocks(); projectors_and_counterexamples(); z3_checks(); transport_audit()
    out=ROOT / 'data' / 'occupied_state_audit.json'
    out.write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps(results,indent=2))
