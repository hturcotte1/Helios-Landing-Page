"""Compare the digests of skeptic_dvec.py (independent code) with data/collisions_python.txt for n <= NMAX."""
import sys
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src/agents')
from skeptic_dvec import levels, transpose
P = (1 << 61) - 1
ref = {}
for line in open('/home/user/Helios-Landing-Page/syt_tree/data/collisions_python.txt'):
    if line.startswith('#') or not line.strip(): continue
    parts = line.split()
    ref[int(parts[0])] = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5]), int(parts[7]))
nmax = int(sys.argv[1])
for n, D in levels(nmax):
    groups = {}
    for lam, v in D.items():
        groups.setdefault(tuple(v), []).append(lam)
    nsym = sum(1 for lam in D if transpose(lam) == lam)
    npairs = (len(D) - nsym) // 2
    bad = [g for g in groups.values() if len(g) > 1 and not (len(g) == 2 and transpose(g[0]) == g[1])]
    digest = 0
    for lam, v in D.items():
        pw = 1
        for x in v:
            digest = (digest + x * pw) % P
            pw = pw * 7 % P
    mine = (len(D), len(groups), nsym, npairs, len(bad), digest)
    status = 'MATCH' if ref.get(n) == mine else ('NO-REF' if n not in ref else 'MISMATCH')
    print(n, status, mine, flush=True)
    if bad: print('   COLLISION', bad, flush=True)
