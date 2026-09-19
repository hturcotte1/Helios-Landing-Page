"""Symbolic check of Schur's Pfaffian identity in the form used in Theorem 4.7:
Pf( (x_i - x_j)/(1 - x_i x_j) )_{i,j<=l} = prod_{i<j} (x_i - x_j)/(1 - x_i x_j), l = 2, 4, 6."""
import sympy as sp
def pf(M, idx):
    if not idx: return sp.Integer(1)
    i = idx[0]; s = 0
    for k, j in enumerate(idx[1:]):
        rest = idx[1:k+1] + idx[k+2:]
        s += (-1) ** k * M[i][j] * pf(M, rest)
    return s
for l in (2, 4, 6):
    x = sp.symbols(f'x1:{l+1}')
    M = [[(x[i] - x[j]) / (1 - x[i] * x[j]) if i != j else 0 for j in range(l)] for i in range(l)]
    lhs = pf(M, list(range(l)))
    rhs = sp.prod([(x[i] - x[j]) / (1 - x[i] * x[j]) for i in range(l) for j in range(i + 1, l)])
    print(l, sp.simplify(sp.together(lhs - rhs)) == 0)
