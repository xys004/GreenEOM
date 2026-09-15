"""Distinguish spin-flux, flux-anticommutator and all-bond currents.

Analytic directed-hop derivatives, with t=hbar=1. The raw all-bond
quantity is -2*N times the conventional normalized spin current.
Run: python -B checks/current_definition_audit.py
"""
from pathlib import Path
import contextlib
import io
import json
import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    from selection_rule import famH


def matrices(N, e, r, lam, phi, rho):
    H, Ds, Df, Vb, V2intr, V2rashba = [np.zeros((2*N, 2*N), complex) for _ in range(6)]
    def add(A, i, j, amplitude):
        A[i, j] += amplitude
        A[j, i] += amplitude.conjugate()
    for n in range(N):
        for j, s in enumerate((1, -1)):
            for m in (1, 2):
                phase_weight = 1 if m == 1 else 2*rho
                phase = np.exp(1j*phase_weight*phi)
                diagonal = 1 if m == 1 else 1j*e*s/2
                flip = -(r if m == 1 else lam)*s/2*np.exp(-1j*s*(2*np.pi*n/N+m*np.pi/N))
                for spin_flip, amplitude in ((False, diagonal), (True, flip)):
                    i = 2*n+j
                    k = 2*((n+m)%N)+(1-j if spin_flip else j)
                    h = amplitude*phase
                    for A, factor in ((H, 1), (Ds, 1j*s*phase_weight), (Df, 1j*phase_weight), (Vb, 1j*m)):
                        add(A, i, k, factor*h)
                    if m == 2:
                        add(V2rashba if spin_flip else V2intr, i, k, 1j*m*h)
    return H, Ds, Df, Vb, V2intr, V2rashba


results = []
for N, rho in ((8, 0), (8, 1), (10, 0)):
    e, r, lam, fraction = .039, .12, .15, .25
    phi = 2*np.pi*fraction/N
    H, Ds, Df, Vb, Vi, Vr = matrices(N, e, r, lam, phi, rho)
    reference = famH(N, e, r, lam, phi, 0., rho)
    h_error = float(np.linalg.norm(H-reference))
    assert h_error < 1e-13
    energies, U = np.linalg.eigh(reference)
    P = U[:, :N]@U[:, :N].conj().T
    Sigma = np.kron(np.eye(N), np.diag([1., -1.]))
    def ac(A):
        return (Sigma@A+A@Sigma)/2
    def mean(A):
        return float(np.trace(P@A).real)
    A, B, C = -mean(Ds), -mean(ac(Df)), -mean(ac(Vb))
    missing_intrinsic = -mean(ac(Vi)) if rho == 0 else 0.
    assert abs(A-B) < 1e-12
    assert abs(C-A-missing_intrinsic) < 1e-12
    assert np.linalg.norm(ac(Vr)) < 1e-14
    if N == 8 and rho == 0:
        assert abs(C-A) > .0017
    else:
        assert abs(C-A) < 1e-12
    # Independent finite-difference comparison to the existing Hamiltonian.
    derivative_errors = []
    for step in (1e-4, 1e-5, 1e-6):
        fd = (famH(N,e,r,lam,phi,step,rho)-famH(N,e,r,lam,phi,-step,rho))/(2*step)
        derivative_errors.append(float(np.linalg.norm(fd-Ds)))
    assert max(derivative_errors) < 1e-7
    results.append(dict(N=N, Ne=N, rho=rho, e=e, r=r, lam=lam,
        flux_fraction=fraction, spin_flux=A, flux_anticommutator=B,
        all_bond_raw=C, missing_intrinsic=missing_intrinsic,
        relative_difference=(C-A)/A, hamiltonian_residual=h_error,
        fermi_gap=float(energies[N]-energies[N-1]),
        finite_difference_derivative_errors=derivative_errors))

output = dict(method='Analytic derivatives of every directed hopping term; occupied projector from famH.',
              convention='t=hbar=1; all_bond_raw=-2*N*<conventional spin current>',
              cases=results)
path = Path(__file__).resolve().parents[1]/'data/current_definition_audit.json'
path.write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
print(json.dumps(output, indent=2))
print('PASS: three current definitions, missing intrinsic term, restored-phase control.')
