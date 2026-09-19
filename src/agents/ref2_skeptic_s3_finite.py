"""ref2_skeptic_s3_finite.py -- referee #2, independent finite checks of Theorem S3.
sgn(v) computed independently as parity of the sorting permutation (0 on ties)."""
import random, itertools, sys
def sgn(v):
    v = list(v)
    if len(set(v)) < len(v): return 0
    # parity of permutation sorting v into strictly decreasing order = parity of #inversions
    inv = sum(1 for i in range(len(v)) for j in range(i+1, len(v)) if v[i] < v[j])
    return -1 if inv % 2 else 1
PERMS = list(itertools.permutations(range(3)))
def LHS(a, b): return sum(sgn([a[i] + b[s[i]] for i in range(3)]) for s in PERMS)
def RHS(a, b): return sum(sgn([a[i] - b[s[i]] for i in range(3)]) for s in PERMS)

def box_check(BOX):
    cnt = bad = 0
    for A in range(1, BOX+1):
        for B in range(1, BOX+1):
            a = (A+B, B, 0)
            for p in range(0, BOX+1):
                for q in range(0, p+1):
                    cnt += 1
                    if LHS(a, (p, q, 0)) != RHS(a, (p, q, 0)):
                        bad += 1
                        if bad < 5: print("  FAIL", A, B, p, q)
    print("box %d: %d points, %d failures" % (BOX, cnt, bad)); return bad

bad = box_check(16)
bad += box_check(24)
# exhaustive un-normalised: a strictly decreasing in [-5,5]^3, b arbitrary in [-5,5]^3 (ties, negatives, unsorted)
cnt = 0
for a in itertools.product(range(-5, 6), repeat=3):
    if not (a[0] > a[1] > a[2]): continue
    for b in itertools.product(range(-5, 6), repeat=3):
        cnt += 1
        if LHS(a, b) != RHS(a, b):
            bad += 1; print("  FAIL general", a, b)
print("general exhaustive a dec in [-5,5]^3, b in [-5,5]^3: %d cases, cumulative failures %d" % (cnt, bad))
# random large
random.seed(2)
for t in range(200000):
    a = sorted(random.sample(range(-10**6, 10**6), 3), reverse=True)
    b = [random.randint(-10**6, 10**6) for _ in range(3)]
    if random.random() < 0.3: b[random.randrange(3)] = b[random.randrange(3)]  # ties in b
    if random.random() < 0.3: b = [a[0] - a[1] * random.choice([1, -1]), b[1], b[2]]  # near-tie constructions
    if LHS(a, b) != RHS(a, b): bad += 1; print("  FAIL random", a, b)
print("random large: cumulative failures", bad)
# also try a NOT strictly decreasing (ties in a) to see whether the hypothesis is needed
badties = 0
for a in itertools.product(range(-3, 4), repeat=3):
    if a[0] > a[1] > a[2]: continue
    for b in itertools.product(range(-3, 4), repeat=3):
        if LHS(a, b) != RHS(a, b): badties += 1
print("a with ties / unsorted: failures", badties, "(informational; hypothesis needed if > 0)")
print("TOTAL failures within hypotheses:", bad)
