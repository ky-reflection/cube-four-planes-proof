# Exact cell volumes in the candidate's centered vertex chart

Use the cube $[-1,1]^3$ with volumes normalized by $8$. Let $v_0,\dots,v_3$ be its four vertices with positive coordinate product. Label and orient the four centered cuts so that the only positive plane at $v_i$ is $i$. The opposite vertices then have complementary signs. Put

$$
C_{ij}=-n_i\cdot v_j>0\quad(i\ne j).
$$

These twelve entries determine the cuts. The omitted diagonal of the matrix of evaluations is $n_i\cdot v_i=\sum_{j\ne i}C_{ij}$. Positive rescaling of any row does not change the arrangement. This note gives the exact seven antipodal cell volumes throughout the chart and a local bound involving just the four singleton volumes.

## Corner triangles and edge overlaps

Each face has opposite corners $v_i,v_j$, with the other two corners $-v_k,-v_l$. The sign pattern away from four corner triangles is positive on $i,j$. The cut $i$ removes a negative triangle at $v_j$, and the cut $j$ removes one at $v_i$. The cut $k$ adds a positive triangle at $-v_l$, and the cut $l$ adds one at $-v_k$. A triangle is contained in the half-square on its side of the diagonal joining the two adjacent corners, so triangles at opposite corners have disjoint interiors. Only adjacent triangles can overlap, and there are no triple overlaps.

Write $T_{ij}$ for the normalized volume of the pyramid from the origin to the triangle cut by $i$ at $v_j$. Its legs as fractions of a cube edge are $C_{ij}/(C_{ij}+C_{ik})$ and $C_{ij}/(C_{ij}+C_{il})$. Therefore

$$
t_{ij}:=12 T_{ij}=\frac{C_{ij}^2}{(C_{ij}+C_{ik})(C_{ij}+C_{il})}.
\tag{1}
$$

An overlap at the edge $v_i,-v_j$ can occur on at most one of its two incident faces. On the face whose two parity-positive corners are $v_i,v_k$, the two triangles at that edge are cut by $k,l$, respectively. Their legs along the edge, from opposite endpoints, are

$$
\ell=\frac{C_{ki}}{C_{ki}+C_{kj}},
\qquad r=\frac{C_{lj}}{C_{li}+C_{lj}}.
$$

Their perpendicular legs are $a=C_{ki}/(C_{ki}+C_{kl})$ and $b=C_{lj}/(C_{lj}+C_{lk})$. In unit square coordinates the overlap is a triangle of base $(\ell+r-1)_+$. The intersection height is $(\ell+r-1)_+/(\ell/a+r/b)$, obtained by adding the two horizontal distances at that height. Thus its normalized pyramid volume is $(\ell+r-1)_+^2/[12(\ell/a+r/b)]$. The corresponding test on the other face swaps $k,l$ and replaces $\ell+r-1$ by its negative.

For each unordered $\{i,j\}$, let $O_{ij}$ be the volume of this unique overlap if positive and zero at a tie. The overlap on the antipodal edge has the same volume. Write $o_{ij}=12O_{ij}$. The exact volumes, multiplied by $12$, are

$$
\begin{aligned}
s_i&=\sum_{j\ne i}(t_{ji}-o_{ij}),\\
p_{ab}&=2-t_{ab}-t_{ba}-t_{cd}-t_{dc}
       +\sum_{i\in\{a,b\},\,j\in\{c,d\}}o_{ij},
       \qquad \{a,b\}\mid\{c,d\}.
\end{aligned}
\tag{2}
$$

Here $s_i$ denotes the singleton cell with positive index $i$, and $p_{ab}$ denotes one member of the complementary two-positive pair. The volumes of their opposite cells agree. For the first identity, each negative triangle at $v_i$ contributes its area minus its incident overlaps. For each $j\ne i$, subtract the overlap on the edge from $v_i$ to $-v_j$. Its antipode gives the corresponding correction to the complementary triple. For the second identity, the pyramid over the face with parity-positive corners $v_a,v_b$ has normalized volume $1/6$. The assigned pair occupies that face outside the union of the four triangles. Each overlap on this face adds back its area by inclusion-exclusion. If an edge incident to this face has its overlap on the other face instead, that overlap has positive set $\{a,b\}$ and contributes the same amount to the assigned pair through the neighboring pyramid. These are all face contributions since any point belonging to two triangles is an adjacent overlap. In particular $\sum_i s_i+\sum p_{ab}=6$, where the second sum runs over the three complementary splits.

The same face decomposition can be integrated on an outer part of each affine face pyramid. With signed offsets $h_i$, the scaled face at radial parameter $\rho\in[0,1]$ has the required vertex signs whenever $\rho C_{ij}>|h_i|$ for all off-diagonal entries. The legs in the scaled face are of the form $(\rho C_{ij}\pm h_i)/[\rho(C_{ij}+C_{ik})]$. Ratios of two legs of one triangle have no offset term. Multiplying the normalized face formulas by $3\rho^2\,d\rho$ then gives quadratic polynomials for the triangle terms and squares of positive parts of affine functions of $\rho$ for the overlap terms, with constant denominators. The inner pyramid, where a corner sign can change, needs a separate decomposition.

## Polynomial form of the centered inequalities

For one edge with endpoint labels $i,j$ and complementary labels $k,l$, put
$$
(A,B,E)=(C_{ki},C_{kj},C_{kl}),\qquad (C,D,F)=(C_{li},C_{lj},C_{lk}),\qquad \Delta=AD-BC.
$$
The overlap formula is equivalently
$$
o_{ij}=\frac{\Delta^2}{(A+B)(C+D)\,Q},\qquad
Q=\begin{cases}
(A+E)(C+D)+(D+F)(A+B), & \Delta\ge0,\\
(C+F)(A+B)+(B+E)(C+D), & \Delta\le0.
\end{cases}
$$
Indeed, $\ell+r-1=\Delta/[(A+B)(C+D)]$, and the other factor in the denominator comes from $\ell/a+r/b$. At $\Delta=0$ both branches give zero. Thus the entire strict vertex chart can be described by four positive rows of sum one, the signs of six quadratic minors, and polynomial equations of degree at most five for the six auxiliary overlap volumes. Equation (1) gives the twelve triangle volumes by equations of degree at most three. With these auxiliary variables, the seven proposed volume caps in (2) are linear inequalities. This supplies the system for exact box subdivision without multiplying all the denominators together. The overlap and triangle denominators are positive throughout the strict chart.

An example of a subdivision workflow for a different four-hyperplane problem is Soberón's polynomial nonvanishing certificate, which separates the Bernstein coefficient vectors from zero on each retained box. [Four hyperplanes do not always equipartition a mass in $\mathbb R^4$](https://arxiv.org/abs/2608.23312)

## The logarithmic split directions

The edge-separated family has $C_{ij}=r_i w_j$ off the diagonal, with $w_j>0$, $\sum_jw_j=1$. Write a perturbation of the logarithms as $\eta_{ij}$, and remove row constants so that each off-diagonal row sums to zero. There is a unique orthogonal decomposition

$$
\eta_{ij}=X_{ij}+y_j+y_i/3,\qquad \sum_jy_j=0,
$$

where $X$ has zero off-diagonal row and column sums. The second term changes only the edge-separated weights, and its squared norm is $(8/3)\sum_j y_j^2$. For an endpoint pair $i<j$ with complementary indices $k<l$, define the logarithmic edge difference

$$
d_{ij}=X_{ki}+X_{lj}-X_{kj}-X_{li}.
$$

The six row vectors of the map $X\mapsto d$ have Gram matrix, after changing some row signs, equal to the graph Laplacian of $K_{2,2,2}$, whose parts are pairs of complementary edges. The diagonal entries are $4$. A pair of noncomplementary rows has dot product $-1$ after changing signs, and complementary rows have dot product zero. The eigenvalues are $0,4,4,4,6,6$, seen by considering vectors constant on each part and vectors antisymmetric within a part. The nullspace of the map on all twelve entries consists of row and column constants. Consequently

$$
\sum_{i<j}d_{ij}^2\ge4\|X\|^2.
\tag{3}
$$

There is also a sup norm bound $\|X\|_\infty\le\|d\|_\infty/2$. For example, the row and column sums give $24X_{01}=-d_{02}-d_{03}+5d_{12}+5d_{13}$, and simultaneous permutations of the indices give the other entries.

At an edge-separated arrangement the denominator $\ell/a+r/b$ equals $1/(w_i+w_j)$ on either side of the edge. The first-order change in $\ell+r-1$ is $[w_iw_j/(w_i+w_j)^2]d_{ij}$. Therefore the two overlap branches have the same quadratic term

$$
o_{ij}=\frac{w_i^2w_j^2}{(w_i+w_j)^3}d_{ij}^2
       +O(\|\eta\|^3).
\tag{4}
$$

This also shows directly that the overlap functions are continuously differentiable throughout the chart and have a common second-order expansion at any edge-separated arrangement.

## A local bound from the singleton caps

At the candidate put $x=\sqrt2-1$, $b=1/(3+x)$, and $w=(xb,b,b,b)$, with index $0$ the small singleton. Set $W=2x/(1+x)$, and consider

$$
\mathcal F=s_0+W(s_1+s_2+s_3).
$$

Write $W_0=1$, $W_1=W_2=W_3=W$. In the row triangle contribution $\sum_{j\ne i}W_j t_{ij}$, the derivative with respect to the log entry corresponding to $a=C_{ij}$, with other entries $b',c$, includes the edge term

$$
\frac{a b'}{(a+b')^2(a+c)(b'+c)}
\bigl[W_j a(b'+c)-W_k b'(a+c)\bigr],
$$

where $b'=C_{ik}$, and a second term with the other pair. At the candidate the three quantities $W_j C_{ij}\sum_{k\ne i,j}C_{ik}$ in each row agree. Thus the triangle contribution has zero first variation, as do the overlaps by (4). Differentiating the displayed edge terms gives the following quadratic forms for the triangle Hessian.

- For row $0$, the form is $W/16$ times the sum of the three squared pairwise differences of its log entries.
- For a high-index row $i$, let $A_i$ be its small-column log entry minus the mean of its two remaining entries, and $B_i$ half the difference of those two. The form is $2x^2[A_i^2+(1+2x)B_i^2]/(1+x)^4$.

For the second formula, the contribution from the two pairs involving the small column is $2x^2(A_i^2+xB_i^2)/(1+x)^4$. The pair of the two high columns contributes $2x^2B_i^2/(1+x)^3$. A row with zero sum has squared norm $(2/3)A_i^2+2B_i^2$. Hence the triangle Hessian is at most $(129/1000)\|\eta\|^2$.

By (2) and (4), the overlap Hessian is subtracted with coefficient $2(W_i+W_j)w_i^2w_j^2/(w_i+w_j)^3$ for each $d_{ij}^2$. The coefficient is greater than $56/1000$ for each edge type. For these two rational comparisons one may use $4142/10000<x<4143/10000$. In the upper triangle bound the two row norms are at most $3W/16$ and $3x^2/(1+x)^4$. In the lower overlap bound the two coefficients are $2(1+W)x^2 b/(1+x)^3$ and $Wb/2$. Equation (3) then bounds the full Hessian by

$$
\frac{129}{1000}\|\eta\|^2-\frac{224}{1000}\|X\|^2.
\tag{5}
$$

To apply this to the three high singleton constraints, put

$$
A=\frac{1+2x}{2(1+x)^2},\qquad B=\frac{1-x}{4(1+x)^2},
\qquad G=\frac38+\frac{5A-B}{3}.
$$

Differentiating (1) in the row-sum-zero gauge gives, for each high $j$,

$$
ds_j=\frac38\eta_{0j}
     +A\sum_{i\ne0,j}\eta_{ij}
     +B\sum_{i\ne0,j}\eta_{i0}.
$$

On the residual subspace this is $(3/8-A)X_{0j}-B X_{j0}$ and the sum over the high indices is zero. On the edge-weight subspace the sum is $(-8A+16B)y_0/3$, so the tangent equations $ds_1=ds_2=ds_3=0$ force $y_0=0$. For $y_0=0$, the edge-weight contribution to each $ds_j$ is $G y_j$. These facts show that the three edge-weight variations have an invertible derivative onto the three high singleton values, and on their full tangent kernel

$$
\sum_j y_j^2\le
\frac{(3/8-A)^2+B^2}{G^2}\|X\|^2<\frac{11}{1000}\|X\|^2.
$$

The last rational comparison follows from the same interval for $x$, since $A>3/8$, $A>2B$, and the ratio is bounded above by using upper bounds for $A,B$ in its numerator and a lower bound for $5A-B$ in its denominator. Thus $\|\eta\|^2<1.03\|X\|^2$ on the tangent kernel, and (5) is strictly negative there.

It follows that the candidate strictly maximizes the sum of the four singleton volumes locally subject to the three high singleton volumes being at most $m=(1+2\sqrt2)/48$. Indeed, use $X$ and the vector $z_j=s_j-12m$ for the three high indices as local coordinates by the inverse function theorem. The first variation of $\mathcal F$ vanishes. Its quadratic term is negative definite at $z=0$, so near the candidate its change is bounded above by $-c\|X\|^2+C\|z\|^2$ for some $c,C>0$. The total singleton change adds $(1-W)\sum_j z_j$, which for $z_j\le0$ strictly dominates $C\|z\|^2$ when $z\ne0$ is sufficiently small. Thus for any distinct nearby feasible arrangement the singleton total decreases. By (2) the total of the three complementary pair volumes then increases, so some pair exceeds $m$.

This argument is local and centered. The exact formulas (1) and (2) cover the whole centered vertex chart, but the argument does not establish the singleton aggregate bound throughout it.
