"""Phase 1: verify E_lambda(t) = exp(t + t^2/2) P_lambda(t) coefficientwise,
i.e. s_k(lam) = sum_j C(k,j) I_{k-j} d_j(lam), using the brute-force s_k."""
import sys, time
sys.path.insert(0, __import__('os').path.dirname(__file__))
from census import check_identity, s_up, s_up_by_skew, s_from_d
from young import partitions

if __name__ == "__main__":
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    max_k = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    t0 = time.time()
    ok = check_identity(max_n, max_k, verbose=True)
    print("identity holds for all |lam| <=", max_n, "and k <=", max_k, ":", ok, f"({time.time()-t0:.1f}s)")
    # second, fully independent route for s_k (enumerate mu and count skew chains)
    bad = 0
    for n in range(0, 7):
        for lam in partitions(n):
            for k in range(0, 6):
                if s_up_by_skew(lam, k) != s_from_d(lam, k):
                    bad += 1
    print("independent skew-chain route agrees for |lam|<=6, k<=5:", bad == 0)
