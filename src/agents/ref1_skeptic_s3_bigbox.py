"""ref1_skeptic_s3_bigbox.py -- sign vectors of the arrangement H realised in box 32 and at random large integer points,
compared with those realised in box 16 (checks 'every cell has a point in box 16' without using Caratheodory).
Pure Python, exact integers."""
import os
import sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref1_skeptic_s3 import forms, eval_point
H = forms()
def sv(x):
    return tuple((1 if d > 0 else (-1 if d < 0 else 0)) for d in (f[0]*x[0]+f[1]*x[1]+f[2]*x[2]+f[3]*x[3] for f in H))
def box(BX):
    s = set()
    for A in range(1, BX+1):
        for B in range(1, BX+1):
            for p in range(0, BX+1):
                for q in range(0, p+1):
                    s.add(sv((A, B, p, q)))
    return s
s16 = box(16); print("box 16:", len(s16))
s32 = box(32); print("box 32:", len(s32), "new vs box16:", len(s32 - s16))
random.seed(0)
sr = set(); bad = 0
for i in range(300000):
    # mix of generic large points and points with forced coincidences (small differences)
    A = random.randint(1, 10**6); B = random.randint(1, 10**6); q = random.randint(0, 10**6)
    p = q + random.choice([0, random.randint(0, 5), random.randint(0, 10**6)])
    if random.random() < 0.3:
        B = A + random.randint(-3, 3); B = max(B, 1)
    if random.random() < 0.3:
        q = max(0, A + random.randint(-3, 3)); p = q + random.choice([0, A, B, A+B, abs(A-B), random.randint(0, 3)])
    x = (A, B, p, q); sr.add(sv(x))
    l, r = eval_point(x)
    if l != r: bad += 1
print("random large: sign vectors", len(sr), "new vs box16:", len(sr - s16), "; S3 failures:", bad)
