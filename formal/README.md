# Formal certificates and their scope

The manuscript proves a common occupied Slater state under explicit block-occupation hypotheses. This directory certifies components of that argument, not the entire physical theorem.

- SelectionRule.lean and selection_rule_z3.py prove that the half-integer momentum ladder contains both pi/2 and 3pi/2 precisely when N = 2 mod 4. This is an arithmetic equivalence, not an if-and-only-if classification of vanishing currents. Lean also checks the trigonometric phase shift.
- BlockAlgebra.lean proves seven algebraic statements: trace zero, the zero block, its square, anticommutation with sigma_x, the corresponding zero expectation, even-range trigonometric zeros, and absence of occupied-to-empty mixing for a commuting idempotent projector.
- The operator norm is |sin(2q)|. The Frobenius norm is sqrt(2)|sin(2q)|. The manuscript derives the real-space block reduction, and checks/occupied_state_audit.py verifies its symbolic bond identities and numerical reconstruction independently.

Neither Lean file formalizes the real-space Fourier transform, the spectral theorem, the second-quantized Slater-state implication, or the occupation inequalities as a certified numerical enclosure. Those links are supplied as mathematical arguments in Sec. II. Restoring a phase changes the block zero at generic flux; it does not force a nonzero observable response at every parameter point.

Run the arithmetic checks with Python and z3-solver:

    python formal/selection_rule_z3.py

With Lean 4.32.2 and the pinned mathlib available, from formal/:

    lake env lean SelectionRule.lean
    lake env lean BlockAlgebra.lean

Both files print their axiom dependencies. No sorryAx is used. The integration was checked against an existing mathlib v4.32.2 cache; a clean dependency download was not repeated. The local lakefile pins this version and includes both libraries. A fresh setup needs the dependencies and their compiled cache before lake build.

Z3's pi/2 and 3pi/2 encodings use linear integer arithmetic. An alternative congruence with a product of two symbolic integers need not be decidable by that solver; an unknown result would not refute the theorem.
