"""Referee #2 extras: (1) Theorem A to n<=24; (2) the 'Equivalently' chi-formulas in Theorem A; (3) brute-force
Lemma 2.1 (star products), Lemma 2.2(b) relations, and the parity of type (4,2) mentioned in the (I6) proof."""
import os
import sys, itertools
from fractions import Fraction
from math import comb
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from young import partitions, f_hook
from ref2_topend0_theoremA import chi, contents, omega, A_formulas, ff
from ref2_topend0_algebra import compose, transp, cycle_type, J, mult, add, ident

def ff_(n,i):
    r=1
    for j in range(i): r*=n-j
    return r

# (1)+(2)
bad=0; cnt=0
for n in range(0, 25):
    for lam in partitions(n):
        lam=tuple(lam); c=contents(lam); f=f_hook(list(lam))
        C1,C2,C4=sum(c),sum(x*x for x in c),sum(x**4 for x in c)
        A1,A2,A3,A4,_=A_formulas(n,C1,C2,C4)
        om=[omega(lam,r) for r in [(3,),(2,2),(5,),(3,3)]]
        cnt+=1
        if om!=[A1,A2,A3,A4]: bad+=1; print('MISMATCH A',n,lam)
        if n>=3:
            rho=tuple([3]+[1]*(n-3)); lhs=chi(lam,rho); rhs=Fraction(f*(3*C2-Fraction(3*n*(n-1),2)), ff_(n,3))
            if lhs!=rhs: bad+=1; print('MISMATCH chi3',n,lam)
        if n>=4:
            rho=tuple([2,2]+[1]*(n-4)); lhs=chi(lam,rho); rhs=Fraction(f*(C1*C1-3*C2+n*(n-1)), 6*comb(n,4))
            if lhs!=rhs: bad+=1; print('MISMATCH chi22',n,lam)
        if n>=6:  # |K_33| = (n)_6/18
            rho=tuple([3,3]+[1]*(n-6)); lhs=Fraction(ff_(n,6),18)*chi(lam,rho)/f
            if lhs!=A4: bad+=1; print('MISMATCH K33 size',n,lam)
print(f'Theorem A + chi-formulas, all lam |- n <= 24: {cnt} shapes,', 'OK' if bad==0 else f'{bad} BAD')

# (3) Lemma 2.1 for n=6: (i1 k)(i2 k)...(i_mu k) = (i_mu ... i_1 k)
def cyc(seq, n):
    p=list(range(n))
    for a,b in zip(seq, seq[1:]+seq[:1]): p[a]=b
    return tuple(p)
n=6; k=5; bad=0
for mu in range(1,5):
    for seq in itertools.permutations(range(k), mu):
        prod=ident(n)
        for i in seq: prod=mult(prod,{transp(i,k,n):1})
        want=cyc(list(reversed(seq))+[k], n)
        if prod!={want:1}: bad+=1
print('Lemma 2.1 (n=6,k=6, mu<=4):', 'OK' if bad==0 else f'{bad} BAD')
# Lemma 2.2(b): s_k J_k s_k = J_{k+1} - s_k ; s_k J_{k+1} s_k = J_k + s_k ; commutation J_k J_l = J_l J_k
n=6; bad=0
for k in range(1,n):
    s={transp(k-1,k,n):1}
    if mult(mult(s,J(k,n)),s)!=add((1,J(k+1,n)),(-1,s)): bad+=1
    if mult(mult(s,J(k+1,n)),s)!=add((1,J(k,n)),(1,s)): bad+=1
for k in range(1,n+1):
    for l in range(1,n+1):
        if mult(J(k,n),J(l,n))!=mult(J(l,n),J(k,n)): bad+=1
print('Lemma 2.2(a),(b) n=6:', 'OK' if bad==0 else f'{bad} BAD')
# parity of (4,2): sign = (-1)^(n - #cycles)
p=cyc([0,1,2,3],6); p=compose(p,cyc([4,5],6))
print('cycle type', cycle_type(p), 'sign =', (-1)**(6-2), '(the report says (4,2) is odd)')
