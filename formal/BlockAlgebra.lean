import Mathlib

/- Algebraic certificates for the audit. These do not formalize the spectral
theorem, the real-space Fourier transform, or the many-body ground-state claim.
Those links are supplied explicitly in the accompanying derivation. -/
namespace GreenEOMAudit

abbrev M2 := Matrix (Fin 2) (Fin 2) ℂ

def sigmaX : M2 := !![0, 1; 1, 0]
def yz (y z : ℂ) : M2 := !![z, -Complex.I*y; Complex.I*y, -z]
def vBlock (s : ℂ) : M2 := !![0, Complex.I*s; -Complex.I*s, 0]

theorem trace_v (s : ℂ) : Matrix.trace (vBlock s) = 0 := by
  simp [Matrix.trace, vBlock, Fin.sum_univ_two]

theorem v_zero : vBlock 0 = 0 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [vBlock]

theorem v_square (s : ℂ) : vBlock s * vBlock s = s^2 • (1 : M2) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [vBlock, Matrix.mul_apply, Fin.sum_univ_two] <;> ring_nf <;>
    simp [Complex.I_sq]

theorem yz_anticommutes_x (y z : ℂ) :
    yz y z * sigmaX + sigmaX * yz y z = 0 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [yz, sigmaX]

/-- The sigma_x contribution has zero trace against a normalized two-level
spectral projector whose polarization lies in the yz plane. The physical
projector interpretation additionally needs real y,z, r^2=y^2+z^2, r>0. -/
theorem spin_flip_expectation_zero (y z r : ℂ) :
    Matrix.trace (((1 : M2) - r⁻¹ • yz y z) * sigmaX) = 0 := by
  simp [Matrix.trace, Matrix.mul_apply, Fin.sum_univ_two, yz, sigmaX]

/-- Every even-range phaseless Rashba term vanishes at k=pi/2. -/
theorem even_range_zero (m : ℕ) :
    Real.sin ((2 * (m : ℝ)) * (Real.pi / 2)) = 0 := by
  rw [show (2 * (m : ℝ)) * (Real.pi / 2) = (m : ℝ) * Real.pi by ring]
  exact Real.sin_nat_mul_pi m

/-- A perturbation preserving an idempotent occupied projector has no
occupied-to-empty matrix elements. Together with trace(P V)=0 this is the
one-particle certificate that its second quantization annihilates the Slater
determinant. The latter Fock-space implication is not formalized here. -/
theorem no_particle_hole_mixing {R : Type*} [Ring R] (P V : R)
    (hid : P * P = P) (hcomm : V * P = P * V) :
    (1-P)*V*P = 0 := by
  calc
    (1-P)*V*P = V*P - P*(V*P) := by noncomm_ring
    _ = V*P - P*(P*V) := by simp only [hcomm]
    _ = V*P - (P*P)*V := by rw [mul_assoc]
    _ = V*P - P*V := by rw [hid]
    _ = 0 := by rw [hcomm, sub_self]

end GreenEOMAudit

#print axioms GreenEOMAudit.trace_v
#print axioms GreenEOMAudit.v_zero
#print axioms GreenEOMAudit.v_square
#print axioms GreenEOMAudit.yz_anticommutes_x
#print axioms GreenEOMAudit.spin_flip_expectation_zero
#print axioms GreenEOMAudit.even_range_zero
#print axioms GreenEOMAudit.no_particle_hole_mixing
