import sys; sys.path.insert(0, '.')
from ref1_topend2_check import partitions, X_contents, Y_contents, addable_rows, add_row, f_hook
cnt = 0; bad = 0
for m in range(23, 29):
    for nu in partitions(m):
        X = X_contents(nu); Y = Y_contents(nu); fnu = f_hook(nu)
        for i in addable_rows(nu):
            x = (nu[i] if i < len(nu) else 0) - i
            num = m + 1
            for y in Y: num *= abs(x - y)
            den = 1
            for xp in X:
                if xp != x: den *= abs(x - xp)
            if f_hook(add_row(nu, i)) * den != fnu * num: bad += 1; print("FAIL", nu, x)
            cnt += 1
print("Lemma K m=23..28:", cnt, "checks, failures:", bad)
