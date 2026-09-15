"""Arithmetic certificates for the half-integer momentum ladder.
The two special rungs occur iff N = 2 mod 4. This is not a classification
of all current zeros. The physical block and occupation hypotheses are
proved separately in manuscript Sec. III; see formal/README.md.
Universal implications are checked by refutation; existence claims are
checked directly for satisfiability, and the four stated sizes are evaluated.
"""
from z3 import (Int, Solver, ForAll, Exists, Implies, And, Or, Not, sat, unsat,
                get_version_string)

results = []


def prove(name, claim, vars_to_bind):
    """Assert not(claim) and require unsat."""
    s = Solver()
    s.add(Not(ForAll(vars_to_bind, claim)) if vars_to_bind else Not(claim))
    r = s.check()
    ok = (r == unsat)
    results.append((name, ok, r))
    print(f"  [{'PASS' if ok else 'FAIL'}]  {name}")
    if not ok and r == sat:
        print(f"          counterexample: {s.model()}")
    return ok


print("=" * 74)
print(f" Z3 {get_version_string()} -- selection rule arithmetic")
print("=" * 74)
print()

N, l, m = Int('N'), Int('l'), Int('m')

# ---------------------------------------------------------------------------
# Claim 1.  k = pi(2l+1)/N equals pi/2 for some rung l in range  <=>  N = 2 mod 4.
#
#   pi(2l+1)/N = pi/2   <=>   2(2l+1) = N.
# ---------------------------------------------------------------------------
print("Claim 1  the half-integer ladder contains k = +pi/2  iff  N = 2 (mod 4)")
prove(
    "  N>0 :  (exists l, 0<=l<N and 2(2l+1)=N)  <->  N mod 4 = 2",
    Implies(N > 0,
            Exists([l], And(l >= 0, l < N, 2 * (2 * l + 1) == N)) == (N % 4 == 2)),
    [N])

# ---------------------------------------------------------------------------
# Claim 2.  Same for k = -pi/2.  The rungs run over l = 0..N-1, so the angle
#           -pi/2 appears as the rung 3pi/2:
#               pi(2l+1)/N = 3pi/2   <=>   2(2l+1) = 3N,
#           which stays inside linear integer arithmetic.  (Writing the general
#           congruence as 2(2l+1) = N(4m-1) multiplies two variables and leaves
#           the decidable fragment: z3 then answers "unknown", which is a limit
#           of the encoding and not a statement about the ring.)
# ---------------------------------------------------------------------------
print()
print("Claim 2  the ladder contains k = -pi/2 (rung 3pi/2) under the same condition")
prove(
    "  N>0 :  (exists l, 0<=l<N and 2(2l+1)=3N)  <->  N mod 4 = 2",
    Implies(N > 0,
            Exists([l], And(l >= 0, l < N, 2 * (2 * l + 1) == 3 * N))
            == (N % 4 == 2)),
    [N])

# ---------------------------------------------------------------------------
# Claim 3.  The condition is not vacuous and not universal: it must admit rings
#           and exclude rings.  (A condition true for every N would make the
#           paper's "N = 2 (mod 4)" decorative; one true for no N would make
#           the result empty.)
# ---------------------------------------------------------------------------
print()
print("Claim 3  the condition is a genuine restriction (some N in, some N out)")
s = Solver()
s.add(N > 2, N % 2 == 0, N % 4 == 2)
r_in = s.check()
s2 = Solver()
s2.add(N > 2, N % 2 == 0, N % 4 != 2)
r_out = s2.check()
ok = (r_in == sat and r_out == sat)
results.append(("  even rings both satisfying and violating N = 2 (mod 4) exist",
                ok, (r_in, r_out)))
print(f"  [{'PASS' if ok else 'FAIL'}]  even rings both satisfying and violating "
      f"N = 2 (mod 4) exist")

# ---------------------------------------------------------------------------
# Claim 4.  For an EVEN ring the ladder never contains k = 0 or k = pi, so the
#           half-integer ladder is genuinely half-integer: no rung coincides
#           with an integer-momentum value.  (This is what T^N = -1 buys, and
#           the paper leans on it when it says the blocks come in +-k pairs.)
# ---------------------------------------------------------------------------
print()
print("Claim 4  for even N, no rung is an integer multiple of pi")
prove(
    "  N>0 :  never  pi(2l+1)/N = pi*m,  i.e. 2l+1 = N*m has no solution "
    "for even N",
    Implies(And(N > 0, N % 2 == 0),
            Not(Exists([l, m], And(l >= 0, l < N, 2 * l + 1 == N * m)))),
    [N])

# ---------------------------------------------------------------------------
# Claim 5.  Half filling occupies N of the 2N states.  The paper's statement is
#           about half filling specifically; check that N states is an integer
#           count for every ring size, i.e. the condition is well posed for all
#           N (no parity obstruction hiding in "half filling").
# ---------------------------------------------------------------------------
print()
print("Claim 5  half filling of 2N states is an integer occupancy for every N")
prove("  N>0 :  2N is even and N is its half",
      Implies(N > 0, And((2 * N) % 2 == 0, 2 * N - N == N)), [N])

# ---------------------------------------------------------------------------
# Claim 6.  The paper's two exceptional sizes.  It reports the effect at N = 10
#           and cites N = 6 and N = 14 as the neighbouring cases; it reports
#           N = 8 as responding.  Check those four against the rule.
# ---------------------------------------------------------------------------
print()
print("Claim 6  the sizes quoted in the paper fall on the side the paper says")
expect = {6: True, 8: False, 10: True, 14: True}
allok = True
for n, should in expect.items():
    got = (n % 4 == 2)
    good = (got == should)
    allok &= good
    print(f"  [{'PASS' if good else 'FAIL'}]  N={n:2d}: rule says "
          f"{'zeros present' if got else 'zeros absent':11s} , paper says "
          f"{'zeros present' if should else 'zeros absent'}")
results.append(("  quoted ring sizes agree with the rule", allok, None))

print()
print("=" * 74)
npass = sum(1 for _, ok, _ in results if ok)
print(f" {npass}/{len(results)} claims discharged")
print("=" * 74)

assert npass == len(results), "An arithmetic certificate failed"
