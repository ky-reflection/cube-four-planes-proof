# Candidate-threshold exclusions using the six center pyramids

Work in $[-1,1]^3$ with normalized volume. Use the strict alternating vertex signs and the scores $C_{ij}=-n_i\cdot v_j$ from the [affine volume formulas](affine-vertex-volumes-proof.md). Rescale each row so that its largest score is one. Write $x,y$ for its other two scores and sort them as $0<x\le y\le1$ within the row. The offset $H_i$ may range over $(-x,x)$. Set
$$
m=(1+2\sqrt2)/48,\qquad M=0.079759>m.
$$

This proof excludes every arrangement with maximum volume at most $M$ in each of the following domains.

- All six maximal-column patterns when all eight nonpivot scores belong to $[7/10,1]$.
- The reciprocal matching pattern $1032$ when the eight nonpivot scores belong to one of $[.18,.50]$, $[.25,.65]$, $[.30,.74]$, $[.40,.89]$, and $[.45,.93]$.
- The pattern $1032$ when both nonpivot scores of any one row are at most $1/60$, with the other three rows arbitrary. There are also the mixed row criteria below.

The direction coordinates and all four offsets are independent within these conditions. In particular the exclusions include the nonconcurrent arrangements with 15 blocks.

## Exact loss on a center pyramid

Let $F_B$ be the cube face whose default positive pair is $B$ in the alternating chart, and $F_{B^c}$ its opposite. The six pyramids from the cube center to these faces partition the cube up to boundaries. For each cut $i$, orient the score on $F_B$ so that the default sign holds at three corners. Its four corner scores in parallelogram order are $-a,b,d,R$, where $a,b,d$ are its three off-diagonal scores, $R=a+b+d$, and $-a$ and $R$ are opposite. The corresponding signed offset is $h=H_i$ for $i\in B$ and $h=-H_i$ otherwise. On $F_{B^c}$ the score magnitudes are the same and the signed offset is reversed.

Let $P_a(h)$ be the volume in the pyramid to $F_B$ where cut $i$ disagrees with the default sign. At radial scale $\rho$ the bad face area is the probability that
$$
(a+b)U+(a+d)V\le a+h/\rho
$$
for independent uniform $U,V$ on $[0,1]$. Use the triangle area formula with inclusion-exclusion at the edges of the square, and integrate with normalized volume weight $\rho^2/2$ on $[0,1]$. Put $D_a=(a+b)(a+d)$. The result is
$$
P_a(h)=\frac{1}{4D_a}\int_0^1
\big[(a\rho+h)_+^2-(h-b\rho)_+^2-(h-d\rho)_+^2+(h-R\rho)_+^2\big]\,d\rho .
$$
Since $|h|<\min(a,b,d)$, explicit integration gives
$$
P_a(h)=\frac{1}{4D_a}
\left[\frac{a^2}3+ah+h^2-
\begin{cases}
|h|^3/(3a),&h\le0,\\
(1/b+1/d-1/R)|h|^3/3,&h>0.
\end{cases}\right].
\tag{1}
$$
The union bound on the four bad sets in one pyramid yields
$$
V_B\ \ge\ \frac16-\sum_i P_{a_{B,i}}(\sigma_{B,i}H_i),\qquad
\sigma_{B,i}=\begin{cases}1&i\in B,\\-1&i\notin B.\end{cases}
\tag{2}
$$
This uses only the part of the sign cell inside its own pyramid. In particular the losses in (2) are necessary constraints under $V_B\le M$.

## Cap deficits and tetrahedron bounds

Reflect the cube so that the weighted offset sum is nonnegative and the all-positive cell is absent. Let $e_i=1/2-\operatorname{vol}\{n_i\cdot X>H_i\}$. Write $t$ for the volume of the all-negative tetrahedron, which is zero at concurrence, and set $G=7M-1/2$. The fourteen nonempty proper sign cells have deficits $d_J=M-V_J$, and hence
$$
0\le d_J\le M,\quad
\sum_J d_J=2G+t,\quad
e_i=\sum_{J\ni i}d_J-G,\quad
-G\le e_i\le 8M-1/2.
\tag{3}
$$
Also $|\sum_i e_i|\le 2G+3t$ by summing over the cardinalities of the proper subsets.

For a normalized row $(1,x,y)$, the exact cap formula gives
$$
\frac{e_i}{H_i}
=\frac{1+xy/(x+y)-H_i^2/[3(x+y)]}{(1+x)(1+y)},
\tag{4}
$$
with the continuous interpretation at zero. For fixed $H_i$, (4) is decreasing in each of $x,y$ on $|H_i|\le x,y\le1$. For example, the sign of its $x$-derivative is the sign of
$-x(x+2y+xy)+(H_i^2/3)(1+2x+y)$, which is nonpositive since $H_i^2\le x^2$ and $x\le1$. If $x,y\in[\alpha,\beta]$, the lower slope bound starts at $(1+\beta/3)/(1+\beta)^2$ by taking $x=y=\beta$ and $|H_i|\le\beta$. Each lower bound $k$ improves the offset bound to $h_*=\min(\beta,(8M-1/2)/k)$, and hence improves $k$ to $(1+\beta/2-h_*^2/(6\beta))/(1+\beta)^2$. The upper slope bound is $K=(1+\alpha/2)/(1+\alpha)^2$. In particular $e_i/H_i\ge1/3$ throughout the normalized chart, and $|H_i|\le 3(8M-1/2)<1/2$ under the proposed volume cap.

Let $D_i$ be the positive rooted tree weights of the evaluation matrix as in the [direction-box formula](affine-vertex-volumes-proof.md#a-box-with-eight-independently-varying-direction-coordinates). The absolute determinant of the three normal rows other than $i$ equals $D_i/4$, since the determinant of the corresponding three parity vertices has magnitude four. With $B_*=\sum_i D_iH_i\ge0$, the variables $z_j=H_j-n_j\cdot X$ for $j\ne i$ describe the tetrahedron as $z_j\ge0$, $\sum_{j\ne i}D_j z_j\le B_*$. It lies in the cube by the triple-intersection theorem. Accounting for the determinant and normalization gives
$$
t=\frac{B_*^3}{12\prod_i D_i}.
\tag{5}
$$

The rational checker improves an upper bound $T$ on $t$, starting from $T=M$. For fixed $D_i$, replace $H_i$ by the upper bound $e_i/k$ for $e_i\ge0$ and $e_i/K$ otherwise. This is convex and increasing in $e_i$. Maximize the resulting numerator using the capped simplex (3). For an upper bound it suffices to set $t=T$, since adding deficits to reach that total increases all $e_i$. The maximum over the simplex occurs at a vertex, with all but one deficit either zero or $M$.

For a fixed vector of these height upper bounds, maximizing (5) over independent bounds $[D_-,D_+]$ on the four weights needs only endpoints. After taking a cube root, dependence on one weight $z$ has the form $(a z+b)_+z^{-1/3}$ up to a positive factor. Where its numerator is positive the derivative has the sign of $2az-b$. For $a\ge0$ a critical point can only be a minimum. For $a<0$ positivity implies $b>0$, so the expression decreases. Among endpoint choices with a fixed number of $D_+$'s, the largest numerator assigns those to the largest heights. After improving $T$, maximize $\sum_i e_i^2$ on the same vertices. Its bound is convex in $t$ between successive changes of the number of full deficits. Thus one needs $0,T$ and any intermediate breakpoints where $2G+t$ is a multiple of $M$.

## Reciprocal pivots

The pivot pattern $1032$ consists of the reciprocal pairs $01$ and $23$. Keep the four faces for the other two splits $02\mid13$ and $03\mid12$ in (2). In each row the scores appearing as $a_{B,i}$ on these faces are exactly the two nonpivot scores $x,y$, each on two opposite faces. Summing (1) on these faces, with $h=|H_i|$, gives the row loss
$$
E(x,y,h)=
\frac{x^2/(1+x)+y^2/(1+y)}{6(x+y)}
+\frac{(2+x+y)h^2}{2(x+y)(1+x)(1+y)}\left(1-\frac{\gamma h}{6}\right),\qquad
\gamma=1+1/x+1/y-1/(1+x+y).
\tag{6}
$$
Consequently a counterexample must have $\sum_i E_i\ge 2/3-4M$.

The [row checker](affine_pyramid_row_certificate.py) bounds (6) directly with tensor Bernstein coefficients. Put $x=u y$, $h=s x$ and cancel the common factor $y$ in the denominator. All remaining denominators are positive on the closed parameter box. The loss increases in $h\in[0,x]$, since the derivative of the offset term is proportional to $h(1-\gamma h/4)$ and $\gamma h<3$. Thus a bound on $s$ can be used at its largest permitted value. The checker also uses (4) with the cap constraint when $s$ varies freely. The following are upper bounds on each row loss, obtained with zero unresolved boxes.

| Row condition | Upper bound |
| --- | ---: |
| $h\le (307/2500)x$, or $x\le79/500$ with any $h\le x$ | $0.0869$ |
| $y\le1/60$ with any $h\le x$ | $0.0085$ |
| $y\le27/200$ with any $h\le x$ | $0.0600$ |
| Any row under the necessary cap condition $|e_i|\le 8M-1/2$ | $0.1130$ |

For instance, one row in the second line can be combined with three unrestricted rows. Two rows in the third line can be combined with two unrestricted rows. One row in the third line, two in the first line and one unrestricted row also suffice. Each sum is less than $2/3-4M$.

To obtain a larger region away from zero, use a common interval $[\alpha,\beta]$ for the eight nonpivot scores. Each tree weight is between $4\alpha(1+\alpha)^2$ and $4\beta(1+\beta)^2$. Indeed the tree formula is coordinatewise increasing, and at equal nonpivot scores $w$ the cofactor is $4w(1+w)^2$, by the eigenvalues $4w,2+2w,2+2w$ of the matching Laplacian on the zero-sum subspace.

The centered term in (6) is at most $\beta/[6(1+\beta)]$. Dividing its offset term by $e_i^2$ gives a positive factor independent of $h$ times $(1-\gamma h/6)/(x+y+xy-h^2/3)^2$. This is nonincreasing up to the allowed $h<1/2$, since its derivative is nonpositive when $8h-\gamma h^2\le\gamma(x+y+xy)$, and the right side is at least four. At $h=0$ the ratio is
$$
c(x,y)=\frac{(2+x+y)(x+y)(1+x)(1+y)}{2(x+y+xy)^2}.
$$
It is nonincreasing in each argument. For the $x$-derivative, put $A=1+x,\ B=1+y$. Nonpositivity is equivalent to nonnegativity of
$(AB+1)(A+B)(A+B-2)-2A(AB-1)(A+B-1)$, whose degree $(3,3)$ Bernstein coefficients on $0\le x,y\le1$ are all nonnegative. Thus the total row loss is at most
$$
\frac{4\beta}{6(1+\beta)}+c(\alpha,\alpha)\max_{t\le T}\sum_i e_i^2.
$$
The [main checker](affine_pyramid_matching_certificate.py) compares this rational upper bound to $2/3-4M$ for the five intervals stated above.

## All pivot patterns near the regular directions

Now keep all six inequalities (2). Suppose every nonpivot score is at least $\alpha$. The tree bounds use $D_+=16$. For $D_-$, enumerate the three pivot choices in each of the three nonroot rows, substitute $\alpha$ for all nonpivot scores and take the minimum of the resulting tree sums. This works for every pivot pattern. Apply (3)–(5) to obtain $t\le T$ and $\sum H_i^2\le Q$.

Consider three representatives $B$ of the opposite face pairs. Define
$$
L_B=\sum_i \sigma_{B,i}\frac{a_{B,i}}{4D_{a_{B,i}}}H_i,\qquad
S_B=\sum_i\sigma_{B,i}H_i,\qquad S_0=\sum_i H_i.
$$
Both signs of $L_B$ are bounded above by the centered and even terms of (1) less $1/6-M$. The sum of the three centered row losses, before division by twelve, is
$$
\sum_a \frac{a^2}{(a+b)(a+d)}=1-\frac{2abd}{(a+b)(a+d)(b+d)}
\ \le\ 1-\frac{\alpha}{(1+\alpha)^2}.
$$
For the inequality, expand $(a+b)(a+d)(b+d)/(abd)$ as two plus the six ordered ratios. This is convex in each score with the others fixed, so its maximum on $[\alpha,1]^3$ is at endpoints. The maximum is $2(1+\alpha)^2/\alpha$.

In each row one representative face uses the pivot and two use nonpivot scores. Let $h=|H_i|$. On the pivot face the quadratic and cubic terms together are bounded above by $q_p(h^2-h^3/3)$, where $q_p=1/[4(1+\alpha)^2]$, since the cubic coefficient in brackets is at least $1/3$ for either sign. On a nonpivot face the corresponding bound is $q_n(h^2-h^3/(3\alpha))$, where $q_n=1/[8\alpha(1+\alpha)]$. For the positive sign in (1), one of $b,d$ equals one, so $1/b+1/d-1/R\ge2-1/(2+\alpha)\ge1/\alpha$ for $\alpha=7/10$. For the negative sign, replace the other nonpivot score in the denominator by $\alpha$. The factor $(1-h/(3a))/[(1+a)(a+\alpha)]$ is then decreasing in $a\in[\alpha,1]$. The logarithmic derivative of the numerator is at most $h/[\alpha(3\alpha-h)]<1$ by the cap bound on $h$, and the denominator's logarithmic derivative is at least one. Its maximum is thus at $a=\alpha$.

Put $q_*=q_p+2q_n$ and $r_*=q_p+2q_n/\alpha$. Since $\sum |H_i|^3\ge(\sum H_i^2)^{3/2}/2$, a necessary inequality is
$$
\sum_B |L_B|\le A(Q):=
\frac{1-\alpha/(1+\alpha)^2}{3}
+q_*Q-r_*Q^{3/2}/6-3(1/6-M).
\tag{7}
$$
The combined quadratic and cubic term is increasing on $0\le Q\le1$, since $r_*<4q_*$. This justifies substituting an upper bound on the actual squared norm.

To tighten $Q$, compare the coefficients $a_{B,i}/(4D_{a_{B,i}})$ with a common reference $k_*>0$. Their pivot interval is $[1/16,1/(4(1+\alpha)^2)]$. Their nonpivot interval is $[\alpha/(4(1+\alpha)^2),1/(4(1+\sqrt{\alpha})^2)]$. The upper nonpivot bound follows by minimizing $(a+1)(a+b)/a=1+b+a+b/a$ with $b\ge\alpha$. If the two intervals have maximum deviations $\Delta_p,\Delta_n$ from $k_*$, the norm of the three coefficient errors in any row is at most $\nu=\sqrt{\Delta_p^2+2\Delta_n^2}$. The triangle inequality on the three contrast coordinates and (7) give
$$
\sqrt{\sum_B S_B^2}\le \frac{A(Q)+2\nu\sqrt Q}{k_*}.
$$
For the mean coordinate, use $|\sum e_i|\le2G+3T$. The slopes in (4) are between $k_Q=\max(k,3/8-Q/24)$ and $K$, so
$$
|S_0|\le \frac{2G+3T+(K-k_Q)\sqrt Q}{(K+k_Q)/2}.
$$
The four sign rows defining $S_0$ and the $S_B$'s are orthogonal, with squared length four. This gives a new upper bound on the squared offset norm from the two preceding inequalities. At $\alpha=7/10$, the rational checker gets $T\le0.0008565925$, then ten refinements give $Q\le0.0284310223$. The right side of (7) is strictly negative at its certified bound, a contradiction.

## Reproduction

From this directory run

~~~sh
python3 -B affine_pyramid_row_certificate.py
python3 -B affine_pyramid_matching_certificate.py
~~~

All comparisons use rational arithmetic. Square roots and repeated upper bounds are rounded outward on a rational grid. The pyramid loss formula was also checked by exact integration on 160 rational arrangements, together with the six pyramid lower bounds and the central volume identity against the exact cell volumes. The exclusions establish a strict margin over the candidate threshold in the stated domains. The [reciprocal-pivot proof](reciprocal-pivot-proof.md) also covers pattern $1032$ when all nonpivot scores are at most $347/500$, with no positive lower bound on them, and gives exclusions that only restrict two of the directions. The [six affine pivot covers](reciprocal-pivot-proof.md#affine-candidate-charts) complete the alternating chart.
