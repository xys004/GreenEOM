"""Exact four-mode illustration of Supplementary Sec. S6, not a proof by sampling.

Run: python -B checks/two_block_ordering.py. Requires SymPy.
"""
from itertools import combinations
import sympy as s

x=s.symbols('lambda',real=True)
levels=[-3+x,-1-x,s.Integer(1),s.Integer(3)]
V=s.diag(1,-1,0,0)
P=s.diag(1,1,0,0)
assert s.trace(P*V)==0
assert (s.eye(4)-P)*V*P==s.zeros(4)
energies=[s.expand(sum(levels[i] for i in occupied)) for occupied in combinations(range(4),2)]
assert energies[0]==-4
canonical=s.reduce_inequalities([q>=energies[0] for q in energies],x)
fixed_mu=s.reduce_inequalities([levels[0]<0,levels[1]<0,levels[2]>0,levels[3]>0],x)
assert canonical.as_set()==s.Interval(-2,4)
assert fixed_mu.as_set()==s.Interval.open(-1,3)
assert min(q.subs(x,5) for q in energies)==-5
assert energies[0].subs(x,5)==-4
print('Exact two-particle sector energies:',energies)
print('Fixed-Ne ground-state region:',canonical)
print('Sufficient fixed-mu=0 region:',fixed_mu)
print('At lambda=5: selected energy -4, ground-state energy -5.')
print('PASS: exact occupation intervals and persistent excited eigenstate.')
