# Formal verification of the selection rule

The paper's central claim is a selection rule: at half filling of a ring with
`N ≡ 2 (mod 4)`, and at couplings weak enough that the Fermi level splits only
the `k = ±π/2` blocks, the persistent spin current of the uncorrected
buckled-lattice model is exactly independent of `λ_R2`.

The argument has four steps:

1. the generalised translation `T = T₁ ⊗ exp(−i σ_z π/N)` obeys `T^N = −1`, so
   its eigenvalue phases form the half-integer ladder `k = π(2ℓ+1)/N`;
2. with the Peierls phase **omitted** on the range-two hops, the coupling
   restricted to the block at `k` has norm `√2 |sin 2k|`, which vanishes exactly
   at `k = ±π/2`;
3. therefore the current loses its `λ_R2` dependence exactly when the ladder
   *contains* `k = ±π/2`, and that happens if and only if `N ≡ 2 (mod 4)`;
4. restoring the phase shifts the argument to `2k − 2φ`, so at `k = π/2` the
   norm becomes `√2 |sin 2φ|` and the zero is destroyed.

Steps 1, 2 and 4 are linear algebra on explicit matrices, and they are checked
numerically in `../checks/js_definition_and_pi2.py`, which also exhibits the
closed form of step 2 against the computed block norms.

**Step 3 is arithmetic**, and it is the step the condition on `N` rests on.
Arithmetic over the integers is decidable, so this directory proves it rather
than sampling it. The trigonometric collapse in step 4 is proved here too.

## What is in here

| File | Tool | Content |
|---|---|---|
| `selection_rule_z3.py` | Z3 | Six claims, each discharged by asserting its negation and requiring `unsat` |
| `SelectionRule.lean` | Lean 4 + mathlib | The same arithmetic as theorems, plus the trigonometry of step 4 |
| `lakefile.toml`, `lean-toolchain` | Lake | The mathlib dependency, pinned |

The two tools overlap deliberately. They are independent implementations of the
same statement, and the point of running both is that they agree.

## What is proved

Both tools establish the two arithmetic facts:

- the ladder contains `k = +π/2` — that is, `π(2ℓ+1)/N = π/2` for some rung
  `ℓ < N` — **if and only if** `N ≡ 2 (mod 4)`;
- the same condition governs `k = −π/2`, which appears on the ladder
  `ℓ = 0, …, N−1` as the rung `3π/2`. The two zeros therefore enter and leave
  together, so the mechanism is not an artefact of choosing one sign of `k`.

Lean additionally proves the trigonometry of step 4: `sin(2 · π/2) = 0`
phaseless, `sin(π − 2φ) = sin 2φ` with the phase restored, and hence that the
block does **not** vanish at `k = π/2` for any `0 < φ < π/2`.

Z3 additionally checks that the condition is a genuine restriction — even rings
exist on both sides of it — that no rung of the ladder is an integer multiple of
`π`, and that the four ring sizes quoted in the paper (`N = 6, 8, 10, 14`) fall
on the side the paper says they do.

## What is *not* proved here

The block-norm identity `‖V(k)‖ = √2 |sin(2k − 2φ)|` itself. That is a statement
about explicit `2N × 2N` matrices, and it is verified numerically in
`../checks/js_definition_and_pi2.py`, where the computed block norms are printed
against the closed form. What this directory adds is that *given* that identity,
the condition on `N` follows exactly, with no appeal to the sizes swept.

## Running them

Z3 (`pip install z3-solver`):

```
python selection_rule_z3.py
```

Every claim prints `PASS`; a failure would print a counterexample. Run with
z3 4.16.0.

Lean, against an existing mathlib v4.32.2 build. This is the route that was
actually used: the file was placed in a Lean package whose `lean-toolchain` is
`v4.32.2` and whose mathlib was already built, and compiled with

```
lake env lean SelectionRule.lean
```

It prints the axiom dependencies of each theorem and nothing else. The proofs
rest only on `propext`, `Classical.choice` and `Quot.sound`; `sorryAx` does not
appear, so nothing is admitted without proof.

From scratch, `lakefile.toml` pins the same mathlib revision:

```
lake exe cache get
lake build
```

That route fetches mathlib's prebuilt object files rather than compiling the
library, and was not re-run end to end for this deposit.

## One encoding note

Writing the `k = −π/2` case as the general congruence `2(2ℓ+1) = N(4m−1)`
multiplies two variables and leaves the decidable fragment of integer
arithmetic; z3 then answers `unknown`. That is a limit of the encoding and says
nothing about the ring. Stating the same fact as the rung `3π/2`, i.e.
`2(2ℓ+1) = 3N`, keeps it linear and it is discharged immediately. The Lean proof
takes the same route, by `omega`.
