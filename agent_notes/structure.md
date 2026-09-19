# Structure angle: handles on $\lambda\mapsto P_\lambda(t)$, and Conjecture B for shapes with at most 4 rows

Agent: "structure" angle. Scripts (all in `src/agents/`, all exact integer / `Fraction` arithmetic):
`structure_explore.py` (first exploration: Hermite/character expansion, $P_\lambda(-1)$ table, three-row Pfaffian
identity, separation of $3$-row shapes by content moments), `structure_verify.py` (parts (1) Theorem S3 brute force,
(2) content formulas, (5) Pfaffian $\ell=6$, (6) bilinear test, (7) the failed $\Psi$ test; its parts (3)-(4) use
`sympy` root finding and were superseded by the faster exact scripts), `structure_verify_rest.py` (driver for (5)-(7)),
`structure_recover_fast.py` and `structure_recover_small.py` (recovery algorithms of Theorems S1/S2 on all
$\le3$-row shapes with $n\le300$ and all $\le4$-row shapes with $n\le150$), `structure_extra.py` (derangement
transform, $P_\lambda(-1)$ tables, the $\Psi$ special-value test). Every VERIFIED statement below names its script. Notation as in `BRIEFING.md`
and `results.md`: $P_\lambda(t)=\sum_j d_j(\lambda)t^j/j!$, $u_i=d_{n-i}$, $C_k=\sum_{\square}c(\square)^k$.

Labels: **PROVED** = complete proof written here (and machine-checked on the stated range);
**VERIFIED** = computational statement with the stated range; **CONJECTURAL** = evidence only.

## 0. Summary

1. **Theorem S1 (PROVED).** If $\lambda,\mu\vdash n$ both have at most $3$ rows and $d(\lambda)=d(\mu)$, then
   $\lambda=\mu$ when $n\ge10$, and $\mu\in\{\lambda,\lambda^t\}$ for $n\le9$ (exhaustive check of the same
   $(n,C_1^2,C_2)$ argument; also covered by Theorem 4.1). The proof uses only the invariants $n$, $C_1^2$, $C_2$
   (Cor. 4.10.2) and a sign lemma: every shape with $\le3$ rows and $n\ge10$ has $C_1>0$.
2. **Theorem S2 (PROVED).** Same for at most $4$ rows: if $\lambda,\mu\vdash n\ge21$ both have $\le4$ rows and
   $d(\lambda)=d(\mu)$ then $\lambda=\mu$. Uses $n,C_1^2,C_2,C_4$; the extra ingredient is that Newton's
   identities make $6p_5+15p_4$ an affine function of $e_4$ with slope $-30(n-8)$.
   (By transposition the same holds for $\le3$ / $\le4$ columns.) Beyond $4$ rows the even content moments
   $C_2,C_4$ provably determined by $d$ are not enough; one more even-degree invariant is needed per row.
3. **Theorem S3 (PROVED).** $n!\,P_\lambda(t)=\sum_{\tau\in S_n}\chi^\lambda(\tau^2)\,(1+t)^{\mathrm{fix}(\tau)}$.
   So $P_\lambda$ is the fixed-point generating function of the square map, weighted by $\chi^\lambda$; the
   natural variable is $s=1+t$, and $t=-1$ gives the **derangement value**
   $n!P_\lambda(-1)=\sum_{\tau\ \text{derangement}}\chi^\lambda(\tau^2)\in\mathbb Z$
   ($=D_n$ for $\lambda=(n)$, $=(-1)^n(n-1)$ for $(n-1,1)$). Equivalent integer form:
   the **derangement transform** $\delta_m(\lambda)=\sum_{\tau\in\mathrm{Der}_m}\chi^\lambda(\tau^2\cup1^{n-m})$
   satisfies $i!\,u_i=\sum_m\binom im\delta_m$; so the $d$-vector is the binomial transform of $(\delta_m)$,
   $\delta_0=\delta_2=f^\lambda$, $\delta_1=0$, $\delta_3=2\chi^\lambda(3,1^{n-3})$.
4. **Theorem S4 (PROVED).** Pfaffian expansion along the last index: for even $\ell>\ell(\lambda)$,
   $P_\lambda=\sum_{j<\ell}(-1)^{j+1}E_{\alpha_j-1}(t)\,\mathrm{Pf}(B^{(\hat j\hat\ell)})$ with $E_k(t)=\sum_{i\le k}t^i/i!$;
   for three rows
   $$P_{(a,b,c)}=E_c\,P_{(a+2,b+2)}-E_{b+1}\,P_{(a+2,c+1)}+E_{a+2}\,P_{(b+1,c+1)} .$$
5. **Question (c) (VERIFIED negative, $n=7..10$):** there is no identity
   $\sum_cP_{\lambda-c}(t)P_{\lambda-c}(u)=\sum_{a,b}\beta_{ab}P_\lambda^{(a)}(t)P_\lambda^{(b)}(u)$ with constants
   $\beta_{ab}$ (the linear system is inconsistent), so the multiset of children polynomials is not obtained
   from $P_\lambda$ by a constant-coefficient bilinear differential rule. In the verified range ($n\le75$)
   the multiset *is* of course determined, because $\lambda$ is.
6. A false lead, recorded so nobody repeats it: the polynomial $\Psi_\lambda(s)=\sum_m s^mI_mw_m(\lambda)$
   (see §3.4) is **not** a Kronecker multiplicity ($s^{\mathrm{fix}\,\pi}$ is not the character of
   $(\mathbb C^s)^{\otimes n}$, $s^{\#\mathrm{cycles}}$ is); `structure_extra.py` (9) shows $\Psi_\lambda(-1)\ne(-1)^n$
   and `structure_verify.py` (7) finds $633$ non-integral or negative values of $\Psi_\lambda(s)$, $s\le5$, $n\le10$.

## 1. Conjecture B for shapes with at most 3 or 4 rows

### 1.1 Content power sums in row coordinates

Fix $L\ge\ell(\lambda)$, pad $\lambda$ with zeros to $L$ parts and put $x_i=\lambda_i-i$ ($1\le i\le L$), so
$x_1>x_2>\dots>x_L\ge-L$ (strict because $\lambda_i\ge\lambda_{i+1}$). Let $p_k=\sum_{i=1}^Lx_i^k$ and let
$S_k(m)=\sum_{u=1}^mu^k$ be the Faulhaber polynomial, regarded as a polynomial in $m$ (so it is defined for all
integers $m$). Facts about $S_k$ ($k\ge1$): (i) $S_k(m)-S_k(m-1)=m^k$ as polynomials, hence
$S_k(m)-S_k(m')=\sum_{u=m'+1}^mu^k$ for all integers $m\ge m'$; (ii) $S_k(-m)=(-1)^{k+1}S_k(m-1)$ for all integers
$m$. Proof of (ii): $g(m):=S_k(-m)+(-1)^kS_k(m-1)$ satisfies
$g(m)-g(m-1)=-(1-m)^k+(-1)^k(m-1)^k=0$ and $g(1)=S_k(-1)+(-1)^kS_k(0)=S_k(-1)=S_k(0)-0^k=0$, so $g\equiv0$.

**Lemma 1.1 (PROVED).** $C_k(\lambda)=\sum_{i=1}^L\big(S_k(x_i)-S_k(-i)\big)=\sum_{i=1}^LS_k(x_i)+(-1)^k\sum_{i=1}^LS_k(i-1)$.

*Proof.* The boxes of row $i$ have the consecutive contents $1-i,2-i,\dots,\lambda_i-i=x_i$, so by (i) their $k$-th
power sum is $S_k(x_i)-S_k(-i)$ (also for $\lambda_i=0$, where it is $0$). Then use (ii). $\square$

With $S_1(m)=\frac{m(m+1)}2$, $S_2(m)=\frac{m(m+1)(2m+1)}6$, $S_4(m)=\frac{6m^5+15m^4+10m^3-m}{30}$ this gives
$$n=p_1+\tfrac{L(L+1)}2,\qquad C_1=\tfrac{p_2+p_1}2-\sum_{i\le L}\tbinom i2,\qquad
C_2=\tfrac{2p_3+3p_2+p_1}6+\sum_{i\le L}S_2(i-1),\qquad C_4=\tfrac{6p_5+15p_4+10p_3-p_1}{30}+\sum_{i\le L}S_4(i-1).\tag{1.2}$$
Constants: $L=3$: $4,\ 6,\ 18$; $L=4$: $10,\ 20,\ 116$.
**VERIFIED** (`structure_verify.py` (2)): (1.2) holds for every $\lambda$ with $\le L$ rows, $n\le30$, $L=3,4,5$.

**Lemma 1.3 (sign of $C_1$; PROVED).** Row $i$ of length $\ell$ contributes $\ell(\ell+1-2i)/2$ to $C_1$. Hence for
$\lambda=(a,b,c)$ with $\le3$ rows, $2C_1=a(a-1)+b(b-3)+c(c-5)\ge a(a-1)-8$; for $\lambda=(a,b,c,d)$ with $\le4$ rows,
$2C_1\ge a(a-1)-20$ (minima of $\ell(\ell-3),\ell(\ell-5),\ell(\ell-7)$ over integers are $-2,-6,-12$). Consequently
$C_1(\lambda)\le0$ forces $a\le3$, $n\le9$ ($\le3$ rows), resp. $a\le5$, $n\le20$ ($\le4$ rows).
**VERIFIED** exhaustively (`structure_recover_fast.py`, $n\le300$ resp. $n\le150$): the shapes with $\le3$ rows and
$C_1\le0$ all have $n\le9$, those with $\le4$ rows and $C_1\le0$ all have $n\le16$ (lists in §1.2, §1.3).

### 1.2 Three rows

**Theorem S1 (PROVED).** Let $\lambda,\mu\vdash n$ have at most $3$ rows each, and $d(\lambda)=d(\mu)$. If $n\ge10$
then $\lambda=\mu$. (For $n\le9$, $\mu\in\{\lambda,\lambda^t\}$ by Theorem 4.1 of `results.md`.)

*Proof.* By Corollary 4.10.2 of `results.md` (PROVED there), the $d$-vector determines $C_1^2$ and $C_2$; so
$C_1(\lambda)^2=C_1(\mu)^2$ and $C_2(\lambda)=C_2(\mu)$. As $n\ge10$, Lemma 1.3 gives $C_1(\lambda),C_1(\mu)>0$, hence
$C_1(\lambda)=C_1(\mu)$. Take $L=3$ in (1.2): $p_1=n-6$, $p_2=2C_1+8-p_1$, $p_3=3C_2-18-\tfrac{3p_2+p_1}2$ are
therefore the same for $\lambda$ and $\mu$. Newton's identities $e_1=p_1$, $2e_2=e_1p_1-p_2$,
$3e_3=e_2p_1-e_1p_2+p_3$ then give the same monic cubic $\prod_i(T-x_i)$, hence the same multiset $\{x_1,x_2,x_3\}$,
hence (the $x_i$ being strictly decreasing) the same triple, hence $\lambda_i=x_i+i$ agree. $\square$

**VERIFIED** (`structure_recover_fast.py`, `structure_recover_small.py`): the recovery algorithm of the proof — given
$(n,C_1^2,C_2)$, try both signs of $C_1$, compute $p_1,p_2,p_3\to e_1,e_2,e_3$, find the integer roots of the cubic in
$[-3,n]$, keep strictly decreasing roots $\ge-3$ — was run on all $776{,}525$ shapes with $\le3$ rows and $n\le300$.
It returns exactly $\{\lambda\}$ for every such shape with $n\ge8$; the only shapes where it returns more are
$(2),(1,1),(3),(1^3),(3,1),(2,1,1),(3,2),(2,2,1),(3,3),(2,2,2),(3,3,1),(3,2,2)$ ($n\le7$), and there it returns
exactly $\{\lambda,\lambda^t\}$ (checked programmatically). Hence the conclusion $\mu\in\{\lambda,\lambda^t\}$ of Theorem
S1 holds for **all** $n$ with $\lambda,\mu$ of $\le3$ rows, the range $n\le9$ being covered by this exhaustive
$(n,C_1^2,C_2)$ check (independently of Theorem 4.1). The shapes with $\le3$ rows and $C_1\le0$ are exactly
$(1),(1,1),(2,1),(1^3),(2,2),(2,1,1),(3,1,1),(2,2,1),(3,2,1),(2,2,2),(3,2,2),(3,3,2),(3,3,3)$ (max $n=9$, as
Lemma 1.3 predicts).

### 1.3 Four rows

**Theorem S2 (PROVED).** Let $\lambda,\mu\vdash n$ have at most $4$ rows each, $n\ge21$, and $d(\lambda)=d(\mu)$.
Then $\lambda=\mu$.

*Proof.* By Corollary 4.10.2, $C_1^2,C_2,C_4$ agree; by Lemma 1.3 ($n\ge21$) $C_1>0$ for both, so $C_1$ agrees.
Take $L=4$ in (1.2): $p_1=n-10$, $p_2=2C_1+20-p_1$, $p_3=3C_2-60-\tfrac{3p_2+p_1}2$ agree, hence so do
$e_1,e_2,e_3$ (Newton, as above). Newton's identities for four variables,
$$p_4=e_1p_3-e_2p_2+e_3p_1-4e_4,\qquad p_5=e_1p_4-e_2p_3+e_3p_2-e_4p_1,$$
express $p_4$ and $p_5$ as affine functions of $e_4$ with slopes $-4$ and $-4e_1-p_1=-5e_1$. Hence
$30(C_4-116)=6p_5+15p_4+10p_3-p_1$ is affine in $e_4$ with slope $6(-5e_1)+15(-4)=-30(e_1+2)=-30(n-8)\ne0$.
So $e_4$ agrees, the monic quartics $\prod_i(T-x_i)$ agree, the strictly decreasing integer $4$-tuples agree, and
$\lambda=\mu$. $\square$

**VERIFIED** (`structure_recover_fast.py`, `structure_recover_small.py`): the recovery algorithm from
$(n,C_1^2,C_2,C_4)$ (both signs of $C_1$, solve for $e_4$, integer roots of the quartic in $[-4,n]$, keep strictly
decreasing roots $\ge-4$) was run on all $1{,}014{,}428$ shapes with $\le4$ rows and $n\le150$. It returns exactly
$\{\lambda\}$ for every such shape with $n\ge15$, $n\ne8$; at $n=8$ the slope $-30(n-8)$ vanishes and the algorithm
(deliberately) returns nothing; for $n\le14$ the ambiguous cases (e.g. $(4,4,4,2)$ vs $(4,4,3,3)$ at $n=14$) are
exactly the pairs $\{\lambda,\lambda^t\}$ (checked programmatically for all $n\le20$, $n\ne8$). So the conclusion
$\mu\in\{\lambda,\lambda^t\}$ holds for all $n\ne8$ within $\le4$-row shapes, and $n=8$ is covered by Theorem 4.1.
There are $42$ shapes with $\le4$ rows and $C_1\le0$, the largest being $(4,4,4,4)$ ($n=16\le20$, consistent with
Lemma 1.3, whose bound is not sharp).

*Remarks.* (a) The statement is about two shapes *both* of which have $\le4$ rows (or, transposing both, $\le4$
columns). It does not exclude a collision between a $4$-row shape and a shape with more rows; that would need
Conjecture B in general. (b) Why it stops at $4$ rows: for $L$ rows one needs $p_1,\dots,p_L$, i.e. (besides $n$)
$L-1$ independent content moments. $C_1$ is available through $C_1^2$ plus the sign lemma; $C_2,C_4$ are available
(Cor. 4.10.2). $C_4$ involves $p_5$ and $p_4$, but for $L=4$ Newton's identities tie $p_5$ to $p_4$, which is
what rescues the count. For $L=5$ one would need $C_6$ (or another independent even invariant); Theorem 4.10.4 only
gives $C_6-16C_1C_3$ at level $7$, which involves the odd moment $C_3$ and is not enough by itself. (c) The sign
lemma is the only place where the *discreteness* of $\lambda$ enters; for $n\le9$ (resp. $n\le20$) there are
$\le3$-row (resp. $\le4$-row) shapes with $C_1\le0$, but the exhaustive scan shows the wrong sign never produces a
valid competing shape in these ranges either.

## 2. $P_\lambda$ as a character sum over square roots (Theorem S3)

Classical facts used (both standard, both already used in `results.md` Lemma 4.5): (F1) Pieri:
$s_\nu p_1^{j}=\sum_{\mu\supseteq\nu,|\mu/\nu|=j}f^{\mu/\nu}s_\mu$, so $f^{\lambda/\nu}=\langle s_\lambda,s_\nu p_1^j\rangle$
for $|\lambda/\nu|=j$. (F2) Frobenius: $s_\lambda=\sum_{\rho\vdash n}z_\rho^{-1}\chi^\lambda(\rho)p_\rho$, and
$\langle s_\lambda,p_\rho\rangle=\chi^\lambda(\rho)$. (F3) Frobenius–Schur for $S_N$ (all irreducible characters real,
indicator $1$): for $\pi\in S_N$, $\sum_{\nu\vdash N}\chi^\nu(\pi)=\#\{\tau\in S_N:\tau^2=\pi\}$.

**Theorem S3 (PROVED).** For $\lambda\vdash n$,
$$n!\,P_\lambda(t)=\sum_{\tau\in S_n}\chi^\lambda(\tau^2)\,(1+t)^{\mathrm{fix}(\tau)},\qquad\text{equivalently}\qquad
d_j(\lambda)=\frac{j!}{n!}\sum_{\tau\in S_n}\chi^\lambda(\tau^2)\binom{\mathrm{fix}(\tau)}j .$$

*Proof.* Let $F_N=\sum_{\nu\vdash N}s_\nu$. By (F1), $d_j(\lambda)=\sum_{\nu\vdash n-j}f^{\lambda/\nu}=\langle s_\lambda,p_1^jF_{n-j}\rangle$
(the terms with $\nu\not\subseteq\lambda$ vanish). By (F2) applied to each $s_\nu$ and (F3),
$F_N=\sum_{\sigma\vdash N}z_\sigma^{-1}\big(\sum_\nu\chi^\nu(\sigma)\big)p_\sigma=\frac1{N!}\sum_{\pi\in S_N}\#\{\tau:\tau^2=\pi\}\,p_{\mathrm{type}(\pi)}
=\frac1{N!}\sum_{\tau\in S_N}p_{\mathrm{type}(\tau^2)}$ (using $\#\{\pi\text{ of type }\sigma\}=N!/z_\sigma$).
Hence $p_1^jF_{n-j}=\frac1{(n-j)!}\sum_{\tau\in S_{n-j}}p_{\mathrm{type}(\tau^2\cup1^j)}$ and, by (F2),
$$d_j(\lambda)=\frac1{(n-j)!}\sum_{\tau\in S_{n-j}}\chi^\lambda\big(\tau^2\cup1^{j}\big).$$
Now count pairs $(A,\tau)$ with $A\subseteq[n]$, $|A|=j$, $\tau\in S_n$, $A\subseteq\mathrm{Fix}(\tau)$: for fixed
$A$, $\tau$ is an arbitrary permutation of $[n]\setminus A$, so
$\sum_{\tau\in S_n}\chi^\lambda(\tau^2)\binom{\mathrm{fix}\tau}j=\binom nj\sum_{\tau'\in S_{n-j}}\chi^\lambda(\tau'^2\cup1^j)=\binom nj(n-j)!\,d_j=\frac{n!}{j!}d_j$.
Multiply by $t^j$ and sum over $j$: $\sum_j\binom{\mathrm{fix}\tau}jt^j=(1+t)^{\mathrm{fix}\tau}$. $\square$

**VERIFIED** (`structure_verify.py` (1)): both displayed identities, by brute force over all of $S_n$ with
Murnaghan–Nakayama characters (`check_characters.chi`), for all $\lambda\vdash n\le7$.

**Corollary S3.1 (derangement value; PROVED).** $n!\,P_\lambda(-1)=\sum_{\tau\in\mathrm{Der}_n}\chi^\lambda(\tau^2)\in\mathbb Z$.
For $\lambda=(n)$ this is the derangement number $D_n$; for $\lambda=(n-1,1)$ it is $(-1)^n(n-1)$.

*Proof.* Put $t=-1$ in S3. For $(n)$, $\chi\equiv1$. For $(n-1,1)$, $\chi(\pi)=\mathrm{fix}(\pi)-1$ and
$\mathrm{fix}(\tau^2)=2c_2(\tau)$ ($c_2$ = number of $2$-cycles; a derangement has no fixed points and its square
fixes exactly the points in its $2$-cycles). $\sum_{\tau\in\mathrm{Der}_n}c_2(\tau)=\binom n2D_{n-2}$ (choose the
$2$-cycle, derange the rest). So the value is $n(n-1)D_{n-2}-D_n$. From $D_m=mD_{m-1}+(-1)^m$:
$(n-1)D_{n-2}=D_{n-1}+(-1)^n$, so $n(n-1)D_{n-2}=nD_{n-1}+n(-1)^n=D_n+(n-1)(-1)^n$. $\square$
(Table in `structure_explore.py`: $1,0,1,2,9,44,265,1854,14833,133496$ for $(n)$, $n=0..9$, and
$1,-2,3,-4,5,-6,7,-8$ for $(n-1,1)$, $n=2..9$; also `structure_extra.py` (10) for hooks and two-row shapes to $n=12$.)

**Corollary S3.2 (derangement transform; PROVED).** Define, for $0\le m\le n$,
$\delta_m(\lambda)=\sum_{\tau\in\mathrm{Der}_m}\chi^\lambda(\tau^2\cup1^{n-m})\in\mathbb Z$ ($\delta_0=f^\lambda$). Then
$$ i!\,u_i(\lambda)=\sum_{m=0}^i\binom im\delta_m(\lambda)\qquad(u_i=d_{n-i}),\qquad\text{i.e.}\qquad
n!\,P_\lambda(s-1)=\sum_{m=0}^n\binom nm\,\delta_{m}(\lambda)\,s^{\,n-m}.$$
So the $d$-vector and the integer vector $(\delta_m)_{m\le n}$ determine each other (binomial transform), and
$\delta_0=\delta_2=f^\lambda$, $\delta_1=0$, $\delta_3=2\chi^\lambda(3,1^{n-3})$, $\delta_4=6\chi^\lambda(2,2,1^{n-4})+3f^\lambda$
(the $9$ derangements of $4$ points are $6$ four-cycles, whose squares have type $(2,2)$, and $3$ double
transpositions, whose squares are trivial).

*Proof.* In S3, group $\tau\in S_n$ by its fixed-point set $B$ ($|B|=k$): $\tau$ restricted to $[n]\setminus B$ is a
derangement of $n-k$ points, so $\sum_{\mathrm{fix}\tau=k}\chi^\lambda(\tau^2)=\binom nk\delta_{n-k}$ and
$n!P_\lambda(s-1)=\sum_k\binom nk\delta_{n-k}s^k$. Comparing with $n!P_\lambda(s-1)=\sum_j\frac{n!}{j!}d_j(s-1)^j$
and reading off, or directly: $n!d_j/j!=\sum_k\binom kj\binom nk\delta_{n-k}$, and $\frac{j!}{n!}\binom kj\binom nk=\frac1{(k-j)!(n-k)!}$,
so $d_j=\sum_m\frac{\delta_m}{m!\,(n-j-m)!}$, i.e. $(n-j)!\,d_j=\sum_m\binom{n-j}m\delta_m$. $\square$

**VERIFIED** (`structure_extra.py` (8)): $i!u_i=\sum_m\binom im\delta_m$ with $\delta_m$ computed by brute force
over $\mathrm{Der}_m$ and Murnaghan–Nakayama characters, all $\lambda\vdash n\le7$.

*Why this is a handle.* (i) Transpose invariance is visible: $\tau^2$ is always even, so
$\chi^{\lambda^t}(\tau^2)=\chi^\lambda(\tau^2)$. (ii) The level structure of `results.md` §4.10 is visible:
$\delta_m/f^\lambda=\sum_{\rho\vdash m,\ \text{no part }1}\frac{m!}{z_\rho}\,\frac{\chi^\lambda(\mathrm{sq}(\rho)\cup1^{n-m})}{f^\lambda}$
where $\mathrm{sq}(\rho)$ replaces each even part $2k$ by $k,k$; only the classes that are squares of derangements
of $m$ points enter level $m$. (iii) In the variable $s=1+t$ the polynomial $n!P_\lambda(s-1)$ has integer
coefficients $\binom nm\delta_m$; the coefficient of $s^0$ is the derangement value. (iv) Grouping instead by
$\pi=\tau^2$: square roots of $\pi$ preserve $\mathrm{Fix}(\pi)$ and act there as an involution, so
$P_\lambda(t)=\sum_{m=0}^nw_m(\lambda)\,\hat H_m(1+t)$ with $\hat H_m(x)=\sum_{\iota\in\mathrm{Inv}_m}x^{\mathrm{fix}\,\iota}=\sum_i\binom m{2i}(2i-1)!!\,x^{m-2i}$
and $w_m(\lambda)=\frac1{m!(n-m)!}\sum_{\tau\in S_{n-m}\text{ without 1- or 2-cycles}}\chi^\lambda(\tau^2\cup1^m)$
(**VERIFIED** for all $n\le9$, `structure_explore.py` (1)). This is the same identity as
$E_\lambda=e^{t+t^2/2}P_\lambda$ read backwards, and I found no new invariant in it; the attempted
representation-theoretic reading of $\sum_ms^mI_mw_m$ was wrong (summary item 6).

## 3. Pfaffian expansion and the three-row formula (Theorem S4)

Recall Theorem 4.7 of `results.md`: for even $\ell\ge\ell(\lambda)$, $\alpha_i=\lambda_i+\ell-i$,
$P_\lambda=\mathrm{Pf}(B)$, $B_{ij}=\sum_{a\le\alpha_i,b\le\alpha_j}\mathrm{sgn}(\alpha_i-a-\alpha_j+b)\frac{t^{a+b}}{a!b!}$.
Two consequences of the formula for $B_{ij}$: if $\alpha_j=0$ then $B_{ij}=\sum_{a<\alpha_i}t^a/a!=E_{\alpha_i-1}(t)$
with $E_k(t)=\sum_{i=0}^kt^i/i!$; and for $\alpha_i>\alpha_j$, $B_{ij}=P_{(\alpha_i-1,\alpha_j)}$ (Theorem 4.7 with
$\ell=2$, where $\alpha=(\lambda_1+1,\lambda_2)$).

**Lemma 3.1 (Pfaffian expansion along the last index; PROVED).** For an antisymmetric $\ell\times\ell$ matrix
$A$ with $\ell$ even, $\mathrm{Pf}(A)=\sum_{j=1}^{\ell-1}(-1)^{j+1}a_{j\ell}\,\mathrm{Pf}(A^{(\hat j\hat\ell)})$.

*Proof.* $\mathrm{Pf}(A)=\sum_M\mathrm{sgn}(M)\prod_{\{i<k\}\in M}a_{ik}$ over perfect matchings $M$ of $[\ell]$,
$\mathrm{sgn}(M)=(-1)^{\mathrm{cr}(M)}$ with $\mathrm{cr}(M)$ the number of crossing pairs $\{i<k\},\{i'<k'\}$, $i<i'<k<k'$.
Split $M=\{\{j,\ell\}\}\cup M'$. A pair of $M'$ crosses $\{j,\ell\}$ iff exactly one of its endpoints lies in the
open interval $(j,\ell)$. The $\ell-j-1$ points of $(j,\ell)$ are covered by $M'$; those matched inside come in
pairs, so the number of pairs of $M'$ with exactly one endpoint in $(j,\ell)$ has the parity of $\ell-j-1\equiv j+1$
($\ell$ even). Hence $\mathrm{sgn}(M)=(-1)^{j+1}\mathrm{sgn}(M')$ and the sum factors as claimed. $\square$

**Theorem S4 (PROVED).** Let $\ell$ be even with $\ell>\ell(\lambda)$ (so $\alpha_\ell=0$). Then
$P_\lambda(t)=\sum_{j=1}^{\ell-1}(-1)^{j+1}E_{\alpha_j-1}(t)\,\mathrm{Pf}\big(B^{(\hat j\hat\ell)}\big)$. In particular
for $\lambda=(a,b,c)$ ($\ell=4$, $\alpha=(a+3,b+2,c+1,0)$):
$$P_{(a,b,c)}(t)=E_c(t)\,P_{(a+2,b+2)}(t)-E_{b+1}(t)\,P_{(a+2,c+1)}(t)+E_{a+2}(t)\,P_{(b+1,c+1)}(t),$$
where the two-row polynomials are given explicitly by Corollary 4.8 of `results.md`.

*Proof.* Lemma 3.1 with $A=B$, $B_{j\ell}=E_{\alpha_j-1}$; for $\ell=4$, $\mathrm{Pf}(B^{(\hat j\hat4)})$ is the single
entry $B_{ik}$, $\{i<k\}=\{1,2,3\}\setminus\{j\}$, and $B_{12}=P_{(a+2,b+2)}$, $B_{13}=P_{(a+2,c+1)}$, $B_{23}=P_{(b+1,c+1)}$. $\square$

**VERIFIED**: the three-row formula against the $d$-vector recursion for all $3$-row shapes with $n\le16$
(123 shapes, `structure_explore.py` (3)); the general expansion with $\ell=6$ against both the full Pfaffian and the
recursion for all $5$-row shapes with $n\le12$ (`structure_verify.py` (5)).

*Remark.* Iterating, $P_\lambda$ for $\ell(\lambda)\le2k-1$ is an alternating sum of products of $k$ factors
(truncated exponentials and two-row polynomials), i.e. the Pfaffian is a "sum over matchings of the rows"
with the two-row polynomial as edge weight and $E_{\lambda_i+\ell-i-1}$ as the weight of the row matched to the
phantom row $\ell$. I did not find a way to read off the row lengths from this expression directly (the terms
interfere); Theorem S1 instead goes through the top of the vector.

## 4. Question (c): children polynomials

**VERIFIED (negative)** (`structure_verify.py` (6)). For $n=7,8,9,10$ the linear system in the unknowns
$\beta_{ab}$ ($0\le a,b\le n$)
$$\sum_{c}P_{\lambda-c}(t)\,P_{\lambda-c}(u)=\sum_{a,b}\beta_{ab}\,P_\lambda^{(a)}(t)\,P_\lambda^{(b)}(u)\quad\text{for all }\lambda\vdash n$$
(coefficientwise in $t,u$) is inconsistent: $\mathrm{rank}A=(n+1)^2=64,81,100,121$ but $\mathrm{rank}[A|b]=65,82,101,122$
for $n=7,8,9,10$ (exact rational ranks via `sympy.Matrix.rank`). So no fixed constant-coefficient bilinear
differential operator produces the "second moment" of the children from $P_\lambda$; only the first moment
$\sum_cP_{\lambda-c}=P_\lambda'$ is available this way. Since the $P_\lambda$, $\lambda\vdash n$, span a space of
dimension $\le n+1$ while there are $p(n)$ of them, *any* rule of the form "children data $=\Phi(P_\lambda)$" must be
genuinely nonlinear in $P_\lambda$ or use integrality; I have no candidate. (In the verified range $n\le75$ the
multiset $\{P_{\lambda-c}\}$ is trivially determined, because $\lambda$ is, up to transpose, and the children of
$\lambda^t$ are the transposes of the children of $\lambda$.)

## 5. Assessment and what to try next

* Theorems S1–S2 are, as far as I can see, the first proofs of Conjecture B on a class of shapes with arbitrarily
  many corners (all $\le4$-row shapes), but only *within* the class. The mechanism is: even content moments
  determined by $d$ (Cor. 4.10.2) $+$ a sign lemma for $C_1$ $+$ Newton's identities with the discreteness of the
  $x_i$. The obvious continuation is to find one more $d$-determined even-degree invariant that is a polynomial in
  $n,C_1,\dots,C_6$ whose $p_6$/$p_7$ content is independent of $C_1^2,C_2,C_4$ (then $\le5$ rows would follow);
  Theorem 4.10.4's level-7 invariant $C_6-16C_1C_3$ would suffice *if* $C_1C_3$ could be handled, which for
  $\le5$ rows it can: $C_3$ is a polynomial in $p_1,\dots,p_4$, all of which are already known at that point, so
  $C_6-16C_1C_3$ gives $p_7+\tfrac72p_6+\dots$, and Newton ties $p_6,p_7$ to $e_5$ — **CONJECTURAL until the level-7
  statement of Theorem 4.10.4 is proved in full (it is currently supported by rank certificates only)**.
* Theorem S3 gives a clean way to think about the whole vector: $n!P_\lambda(s-1)=\sum_m\binom nm\delta_m(\lambda)s^{n-m}$
  with $\delta_m=\sum_{\tau\in\mathrm{Der}_m}\chi^\lambda(\tau^2\cup1^{n-m})$. Conjecture B says these $n+1$ integers
  (of which $\delta_1=0$, $\delta_2=\delta_0$) determine $\lambda$ up to transpose. A possible use: $\delta_m/f^\lambda$
  is a sum of normalized characters $\hat\chi^\lambda$ at classes with support $\le m$, and the Kerov–Olshanski
  theory makes $\hat\chi^\lambda(\rho\cup1^{n-|\rho|})$ a polynomial in $n$ and contents with leading term controlled
  by the free cumulants; the leading terms might be exploited asymptotically. I did not pursue this.
* Special values: $t=-1$ (derangements) and $t=0$ are the only points where S3 collapses to something classical;
  I found no other special $t$ with combinatorial meaning.
