"""topend: rank of the evaluation map V_m -> Q^{T_N}, T_N = all partitions of n <= N, V_m = span{n^j C_mu: j+|mu| <= m}."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from topend_omega import monomials, evaluate, rank_and_pivot_rows_modp
from topend_mn import C, partitions
for m in (2, 4, 6, 8):
    monos = monomials(m)
    for N in range(8, 21):
        T = [lam for n in range(N + 1) for lam in partitions(n)]
        rows = [[evaluate(lam, mo, {k: C(lam, k) for k in range(11)}) for mo in monos] for lam in T]
        rank, _ = rank_and_pivot_rows_modp(rows)
        print(f"m={m} dim={len(monos)} N={N} |T|={len(T)} rank={rank}", flush=True)
        if rank == len(monos): break
