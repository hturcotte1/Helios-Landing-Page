"""ref2_skeptic_C123.py -- referee #2: independent recheck of Corollaries C1, C2, C3 and the H_k consequence (skeptic2.md Sec.3).
Own code throughout (own d-vector recursion, own GL_3 character arithmetic, own chain counts). Exact integers."""
import sys
from itertools import permutations
from functools import lru_cache
from math import comb
sys.setrecursionlimit(10000)
PERMS = list(permutations(range(3)))

# ---------- own d-vector ----------
@lru_cache(maxsize=None)
def dvec(lam):
    lam = tuple(x for x in lam if x > 0)
    n = sum(lam)
    if n == 0: return (1,)
    out = [0]*(n+1); out[0] = 1
    for i in range(len(lam)):
        if i+1 == len(lam) or lam[i] > lam[i+1]:
            sub = lam[:i] + (lam[i]-1,) + lam[i+1:]
            dv = dvec(tuple(x for x in sub if x > 0))
            for j, v in enumerate(dv): out[j+1] += v
    return tuple(out)

def conj(lam):
    return tuple(sum(1 for x in lam if x > c) for c in range(lam[0])) if lam else ()

# ---------- own GL_3 characters: dict weight->mult ----------
def sgnv(v):
    s = 1
    for i in range(3):
        for j in range(i+1,3):
            d = v[i]-v[j]
            if d == 0: return 0
            if d < 0: s = -s
    return s

def alt(gamma):
    """A_gamma = sum_w sgn(w) x^{w gamma}"""
    out = {}
    for s in PERMS:
        w = tuple(gamma[s[i]] for i in range(3))
        sg = 1
        for i in range(3):
            for j in range(i+1,3):
                if s[i] > s[j]: sg = -sg
        out[w] = out.get(w, 0) + sg
    return {k:v for k,v in out.items() if v}

def mul(P, Q):
    out = {}
    for k1,v1 in P.items():
        for k2,v2 in Q.items():
            k = (k1[0]+k2[0], k1[1]+k2[1], k1[2]+k2[2])
            out[k] = out.get(k,0) + v1*v2
    return {k:v for k,v in out.items() if v}

RHO = (2,1,0)
@lru_cache(maxsize=None)
def char(lam):
    """character of V_lambda for GL_3 (lambda dominant, possibly negative entries) by dividing alternants exactly."""
    num = alt((lam[0]+2, lam[1]+1, lam[2]))
    den = alt(RHO)
    # exact polynomial division in Laurent monomials: repeatedly strip the lexicographically largest term
    quot = {}
    rem = dict(num)
    dlead = max(den)  # (2,1,0)
    while rem:
        lead = max(rem)
        c = rem[lead] // den[dlead]
        assert c * den[dlead] == rem[lead]
        m = (lead[0]-dlead[0], lead[1]-dlead[1], lead[2]-dlead[2])
        quot[m] = quot.get(m,0) + c
        for k,v in den.items():
            kk = (k[0]+m[0], k[1]+m[1], k[2]+m[2])
            rem[kk] = rem.get(kk,0) - c*v
            if rem[kk] == 0: del rem[kk]
    return {k:v for k,v in quot.items() if v}

def decompose(ch):
    """decompose a (virtual) character into irreducibles by highest weight; returns dict lam->mult and total count."""
    ch = dict(ch); out = {}
    while ch:
        lead = max(ch)          # lexicographically max weight is dominant and a highest weight
        c = ch[lead]
        out[lead] = c
        for k,v in char(lead).items():
            ch[k] = ch.get(k,0) - c*v
            if ch[k] == 0: del ch[k]
    return out

def dual(ch):
    return {(-k[0],-k[1],-k[2]):v for k,v in ch.items()}

def racah_count(lam, W):
    a = (lam[0]+2, lam[1]+1, lam[2])
    return sum(m*sgnv((a[0]+b[0], a[1]+b[1], a[2]+b[2])) for b,m in W.items())

# C1 direct test: lam dominant with parts in a box, W = V_mu (irreducible) and W = some non-irreducible products
print("== C1: #const(V_lam x W) vs #const(V_lam x W^*), GL_3, direct decomposition ==")
doms = [(a,b,c) for a in range(0,6) for b in range(0,a+1) for c in range(0,b+1)]
bad = 0; cnt = 0
for lam in doms:
    cl = char(lam)
    for mu in doms:
        W = char(mu); Wd = dual(W)
        dec1 = decompose(mul(cl, W)); dec2 = decompose(mul(cl, Wd))
        assert all(v > 0 for v in dec1.values()) and all(v > 0 for v in dec2.values())
        c1 = sum(dec1.values()); c2 = sum(dec2.values())
        r1 = racah_count(lam, W); r2 = racah_count(lam, Wd)
        cnt += 1
        if not (c1 == r1 and c2 == r2):
            bad += 1; print("  Racah-Speiser count mismatch", lam, mu, c1, r1, c2, r2)
        if c1 != c2:
            bad += 1; print("  C1 FAILS", lam, mu, c1, c2)
print("C1 direct: %d (lam,mu) pairs with parts <= 5, failures = %d" % (cnt, bad))
# also W = V^{(x)j} (tensor powers of the defining rep), lam in box, j<=6, and W = V_mu (x) V_nu (non-irreducible)
V = char((1,0,0)); bad = 0; cnt = 0
for lam in doms:
    cl = char(lam); W = {(0,0,0):1}
    for j in range(1, 7):
        W = mul(W, V)
        c1 = sum(decompose(mul(cl, W)).values()); c2 = sum(decompose(mul(cl, dual(W))).values())
        cnt += 1
        if c1 != c2: bad += 1; print("  C1 FAILS for V^j", lam, j, c1, c2)
print("C1 with W = V^{(x)j}, j<=6, lam parts<=5: %d cases, failures = %d" % (cnt, bad))
bad = 0; cnt = 0
small = [(a,b,c) for a in range(0,4) for b in range(0,a+1) for c in range(0,b+1)]
for lam in doms[:20]:
    cl = char(lam)
    for mu in small:
        for nu in small:
            W = mul(char(mu), char(nu))
            c1 = sum(decompose(mul(cl, W)).values()); c2 = sum(decompose(mul(cl, dual(W))).values())
            cnt += 1
            if c1 != c2: bad += 1; print("  C1 FAILS for V_mu x V_nu", lam, mu, nu)
print("C1 with W = V_mu x V_nu (parts<=3), 20 lam: %d cases, failures = %d" % (cnt, bad))

# ---------- C2: chains through <=3-row partitions ----------
@lru_cache(maxsize=None)
def Nchain(nu, j):
    if j == 0: return 1
    tot = 0
    for i in range(3):
        if i == 0 or nu[i-1] > nu[i]:
            tot += Nchain(nu[:i] + (nu[i]+1,) + nu[i+1:], j-1)
    return tot
print("== C2: N_j(nu) = N_j(nu*) ==")
bad = cnt = 0
for n1 in range(0, 13):
    for n2 in range(0, n1+1):
        for n3 in range(0, n2+1):
            nu = (n1,n2,n3); ns = (n1-n3, n1-n2, 0)
            for j in range(0, 26):
                cnt += 1
                if Nchain(nu, j) != Nchain(ns, j): bad += 1; print("  C2 FAILS", nu, j)
print("C2: nu_1 <= 12, j <= 25: %d cases, failures = %d" % (cnt, bad))
# C2 also equals #const(V_nu x V^j) -- check against decomposition for small cases
bad = cnt = 0
for nu in small:
    W = {(0,0,0):1}
    for j in range(1, 6):
        W = mul(W, V); cnt += 1
        if sum(decompose(mul(char(nu), W)).values()) != Nchain(nu, j): bad += 1; print("  Pieri count mismatch", nu, j)
print("N_j(nu) = #const(V_nu x V^j): %d cases, failures = %d" % (cnt, bad))

# ---------- C3: R_K and the reduction d_j = N_j((3)) / N_j((3,3)) ----------
print("== C3 ==")
bad = 0
for K in range(1, 31):
    u = dvec((3,)*K + (2,2,2)); v = dvec((3,)*K + (1,1,1))
    ut = dvec(conj((3,)*K + (2,2,2))); vt = dvec(conj((3,)*K + (1,1,1)))
    assert u == ut and v == vt
    assert conj((3,)*K+(2,2,2)) == (K+3,K+3,K) and conj((3,)*K+(1,1,1)) == (K+3,K,K)
    okR = all(u[j] == v[j] for j in range(K+1))
    okred = all(u[j] == Nchain((3,0,0), j) and v[j] == Nchain((3,3,0), j) for j in range(K+1))
    first = next((j for j in range(len(u)) if u[j] != v[j]), None)
    # boundary: at j = K+1 the reduction should (generally) break: check whether it does
    red_break = (u[K+1] != Nchain((3,0,0), K+1)) or (v[K+1] != Nchain((3,3,0), K+1))
    if not (okR and okred): bad += 1
    print("K=%2d R_K=%s reduction(j<=K)=%s first_diff=%s reduction breaks at K+1: %s" % (K, okR, okred, first, red_break))
print("C3 failures (K<=30):", bad)

# ---------- H_k: Lemma B formula and window ----------
print("== H_k ==")
def C1C2(lam):
    c1 = c2 = 0
    for i,l in enumerate(lam):
        for c in range(l):
            c1 += c - i; c2 += (c-i)**2
    return c1, c2
def hooks(lam):
    lt = conj(lam); H = []
    for i,l in enumerate(lam):
        for c in range(l):
            H.append(l - c + lt[c] - i - 1)
    return sorted(H)
for k in range(1, 12):
    lam = (5,) + (3,)*k + (2,2,2); mu = (4,4) + (3,)*k + (1,1,1); n = 3*k+11
    u, v = dvec(lam), dvec(mu)
    a, b = dvec((3,)*k + (2,2,2)), dvec((3,)*k + (1,1,1))
    lemB = all(u[j] - v[j] == sum(comb(j,i)*(a[j-i]-b[j-i]) for i in range(0, min(j,2)+1)) for j in range(k+1))
    # Lemma B in its full form for both shapes: d_j(lam) = sum_i C(j,i) D_i(alpha') d_{j-i}(body)  with D=(1,1,1)
    lemB_full = all(u[j] == sum(comb(j,i)*a[j-i] for i in range(0, min(j,2)+1)) for j in range(k+1)) and \
                all(v[j] == sum(comb(j,i)*b[j-i] for i in range(0, min(j,2)+1)) for j in range(k+1))
    D = [j for j in range(n+1) if u[j] != v[j]]
    c1l, c2l = C1C2(lam); c1m, c2m = C1C2(mu)
    print("k=%2d n=%2d hooks_eq=%s f_eq=%s C2_eq=%s C1: %d,%d  LemmaB(diff)=%s LemmaB(full)=%s diff idx == {k+4..n-4}: %s  (mu != lam^t: %s)" % (
        k, n, hooks(lam)==hooks(mu), u[n]==v[n], c2l==c2m, c1l, c1m, lemB, lemB_full, D == list(range(k+4, n-3)), mu != conj(lam)))
