/-
Lean 4 + mathlib formalisation of the arithmetic and trigonometric core of the
selection rule in "Hermiticity and covariance in equations of motion for
spin-orbit-coupled rings".

The physical chain is:

  T = T_1 (x) exp(-i sigma_z pi/N)  obeys  T^N = -1,
  so its eigenvalue phases form the half-integer ladder  k = pi(2l+1)/N;

  with the Peierls phase OMITTED on the range-two hops,
      ||V(k)|| = sqrt 2 * |sin 2k|,
  which vanishes exactly at k = +-pi/2;

  hence the persistent spin current loses its lambda_R2 dependence exactly
  when the ladder CONTAINS k = +-pi/2  --  the paper's condition N = 2 (mod 4);

  restoring the phase shifts the argument to 2k - 2phi, and at k = pi/2 the
  norm becomes sqrt 2 * |sin 2phi|, which is non-zero for 0 < phi < pi/2, so
  the corrected model responds.

Four statements below are the exact content of that chain that is not linear
algebra.  `ladder_pi_two` and `ladder_neg_pi_two` are the arithmetic; the
`phaseless_*` and `restored_*` lemmas are the trigonometry.
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic

namespace SelectionRule

open Real

/-! ### The ladder

The rungs are `k_l = pi (2l+1) / N` for `l = 0, ..., N-1`.  Asking whether a
given rung equals a given angle clears the denominator into a statement about
integers, which is what we state here. -/

/-- The half-integer ladder contains `k = pi/2` -- that is, `pi(2l+1)/N = pi/2`
for some rung `l` -- if and only if `N ≡ 2 (mod 4)`.  This is the paper's
condition on the ring size. -/
theorem ladder_pi_two (N : ℕ) (hN : 0 < N) :
    (∃ l : ℕ, l < N ∧ 2 * (2 * l + 1) = N) ↔ N % 4 = 2 := by
  constructor
  · rintro ⟨l, -, hl⟩
    omega
  · intro h
    exact ⟨(N - 2) / 4, by omega, by omega⟩

/-- The same condition governs `k = -pi/2`, which appears on the ladder
`l = 0, ..., N-1` as the rung `3pi/2`.  The two zeros of `|sin 2k|` therefore
enter and leave together, so the mechanism is not an artefact of picking one
sign of `k`. -/
theorem ladder_neg_pi_two (N : ℕ) (hN : 0 < N) :
    (∃ l : ℕ, l < N ∧ 2 * (2 * l + 1) = 3 * N) ↔ N % 4 = 2 := by
  constructor
  · rintro ⟨l, -, hl⟩
    omega
  · intro h
    exact ⟨(3 * N - 2) / 4, by omega, by omega⟩

/-- The condition is a genuine restriction: it holds for `N = 10` (the ring the
paper plots) and fails for `N = 8` (the ring where the paper reports a
response). -/
theorem ladder_ten_yes : (10 : ℕ) % 4 = 2 := by norm_num

theorem ladder_eight_no : (8 : ℕ) % 4 ≠ 2 := by norm_num

/-! ### The zero of the phaseless coupling

With the phase omitted the block norm is `sqrt 2 * |sin 2k|`. -/

/-- Phaseless, the block norm vanishes at `k = pi/2`. -/
theorem phaseless_zero_at_pi_two :
    Real.sin (2 * (π / 2)) = 0 := by
  rw [show (2 : ℝ) * (π / 2) = π by ring]
  exact Real.sin_pi

/-- Phaseless, it vanishes at `k = -pi/2` as well. -/
theorem phaseless_zero_at_neg_pi_two :
    Real.sin (2 * (-(π / 2))) = 0 := by
  rw [show (2 : ℝ) * (-(π / 2)) = -π by ring, Real.sin_neg, Real.sin_pi, neg_zero]

/-! ### Restoring the phase moves the zero

With the phase restored the argument becomes `2k - 2phi`, and at `k = pi/2` the
sine collapses to `sin 2phi`.  This is the closed form the paper reports for the
residual. -/

/-- At `k = pi/2` the restored coupling has norm governed by `sin 2phi`, not by
zero: `sin(pi - 2phi) = sin 2phi`. -/
theorem restored_at_pi_two (φ : ℝ) :
    Real.sin (2 * (π / 2) - 2 * φ) = Real.sin (2 * φ) := by
  rw [show (2 : ℝ) * (π / 2) - 2 * φ = π - 2 * φ by ring]
  exact Real.sin_pi_sub _

/-- Consequently the zero is destroyed for every flux strictly between the two
values that would restore it: for `0 < phi < pi/2` the block does not vanish at
`k = pi/2`, so the corrected model responds to the coupling. -/
theorem restored_ne_zero {φ : ℝ} (h0 : 0 < φ) (h1 : φ < π / 2) :
    Real.sin (2 * (π / 2) - 2 * φ) ≠ 0 := by
  rw [restored_at_pi_two]
  exact ne_of_gt (Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith))

/-- The two halves stated together at the ring the paper plots: `N = 10` puts
`k = pi/2` on the ladder, the phaseless coupling vanishes there, and any flux in
`(0, pi/2)` makes the restored coupling non-zero at the same momentum. -/
theorem selection_rule_at_ten {φ : ℝ} (h0 : 0 < φ) (h1 : φ < π / 2) :
    (∃ l : ℕ, l < 10 ∧ 2 * (2 * l + 1) = 10)
    ∧ Real.sin (2 * (π / 2)) = 0
    ∧ Real.sin (2 * (π / 2) - 2 * φ) ≠ 0 :=
  ⟨(ladder_pi_two 10 (by norm_num)).2 (by norm_num),
   phaseless_zero_at_pi_two,
   restored_ne_zero h0 h1⟩

end SelectionRule

-- No `sorry` anywhere; confirm the proofs rest only on the standard axioms.
#print axioms SelectionRule.ladder_pi_two
#print axioms SelectionRule.ladder_neg_pi_two
#print axioms SelectionRule.phaseless_zero_at_pi_two
#print axioms SelectionRule.restored_at_pi_two
#print axioms SelectionRule.restored_ne_zero
#print axioms SelectionRule.selection_rule_at_ten
