# Exact minimum for fixed body-diagonal normals and arbitrary offsets

Use the cube $[-1/2,1/2]^3$ and fix the four oriented normals to
$$
(1,1,1),\quad (1,-1,-1),\quad (-1,1,-1),\quad (-1,-1,1).
$$
The offsets in the four plane equations are independent. The minimum of the largest piece volume in this family is
$$
M_{\rm reg}=\frac1{12}-k_*^2+\frac43 k_*^3=0.08131120853\ldots,
$$
where $k_*\in(0,1/10)$ is the unique root in that interval of
$$
168k^3-84k^2-18k+1=0.
$$
Up to reflection through the cube center, equality requires all four offsets to equal $k_*$. In particular the equal-shift minimizer here is a genuine 15-cell cut and exceeds the candidate $(1+2\sqrt2)/48$ for the unrestricted problem.

## Exact pair volume

Write $y_i=n_i\cdot x$, with $\sum_i y_i=0$, and let the offsets be $h_i$. Reflect the cube if necessary to arrange $H=\sum_i h_i\ge0$. The all-positive sign pattern is then absent. It suffices to consider arrangements whose largest cell volume is $M\le1/12$, which includes the centered configuration. The positive cap of each cut has volume in $[1-8M,7M]$. For these normals its volume at offset $h\in[-1/2,1/2]$ is $1/2-3h/4+h^3/3$, by the sum-of-three-uniforms formula. Evaluating at $-1/8$ and $1/4$ proves $-1/8<h_i<1/4$ for all four offsets.

For a pair $ab$ with complementary indices $c,d$, put $A=h_a+h_b$, $C=h_c+h_d$, and use cube coordinates
$$
\alpha=y_a+y_b,\qquad \beta=y_a+y_c,\qquad \gamma=y_a+y_d.
$$
These are the three spatial coordinates up to a permutation, sign changes, and multiplication by two. Thus they form a cube $[-1,1]^3$ with normalized volume. The coordinates $v=\beta+\gamma,\ w=\beta-\gamma$ range over $|v|+|w|\le2$ with planar Jacobian factor $1/2$. At fixed $\alpha$, the sign cell with positive set $ab$ is the rectangle
$$
2h_a-\alpha\le v\le\alpha-2h_b,\qquad
-\alpha-2h_d\le w\le\alpha+2h_c
$$
intersected with this diamond. It is nonempty only for $\alpha\ge A$, because $H=A+C\ge0$. Each corner of the rectangle can lose a right triangle with legs $2(\alpha-1+h_j-h_i)_+$, for $i\in\{a,b\},j\in\{c,d\}$. These triangles do not overlap. The rectangle stays in $[-2,2]^2$, so a point cannot violate two diamond edges. Also the linear form for each diamond edge is at most two at the two neighboring corners of the rectangle, using $|h_i+h_j|<1/2$, which gives the stated triangle legs. At $\alpha=A$ the legs are all zero. Integrating the area in $\alpha\in[A,1]$ gives the exact normalized volume
$$
P_{ab}=\frac{(1-A)^3}{12}+\frac{H(1-A)^2}{8}
-\frac1{24}\sum_{i\in\{a,b\},\,j\in\{c,d\}}(h_j-h_i)_+^3.
$$

Sort the offsets $h_1\le h_2\le h_3\le h_4$, and write $k=H/4$, $h_1+h_2=2k-d$. Here $0\le k<1/4$ and $d\ge0$. Denote the value of $P_{ab}$ at four equal offsets by $P(k)=1/12-k^2+4k^3/3$. The derivative in $A$ of the first two terms of $P_{ab}$ is $-(1-A)(1-A+H)/4$; its magnitude for $A\le2k$ is at least $(1-4k^2)/4\ge3/16$. For the pair $12$ the four differences in the cubic correction are nonnegative, have sum $4d$, and are each at most $3/8$. Consequently
$$
P_{12}\ge P(k)+\frac{21}{128}d.
$$
Put $b=k-h_1\ge0$. Because $h_2\le k+b/3$, one has $d\ge 2b/3$. Since $P_{12}\le1/12$, the bounds on $b$ include
$$
b\le\frac{64}{7}k^2(1-4k/3),\qquad b\le 3/4-3k,\qquad b\le k+1/8.
\tag{1}
$$

## A lower bound for the largest singleton

The region with $y_j\le h_j$ for $j=2,3,4$ consists of the singleton with positive sign 1 and the central tetrahedron. In coordinates $z_j=h_j-y_j$, $j=2,3,4$, it has normalized volume $T(s_2,s_3,s_4)$ for $s_j=1/2+h_j$, where $T$ is one quarter of the ordinary volume of
$$
z_j\ge0,\qquad z_j+z_l\le s_j+s_l\quad(j\ne l).
$$
Indeed the paired inequalities $|y_j+y_l|\le1$ describe the cube and their other sides are automatic here. The determinant of the three normal rows has absolute value four. The central tetrahedron has volume $H^3/24=8k^3/3$ because its coordinates satisfy $z_j\ge0$ and $z_2+z_3+z_4\le H$. It lies wholly in the cube: at each tetrahedron vertex the pair sums have absolute value at most $1/2$, using $-1/8<h_i<1/4$. This also gives the correct zero volume when $H=0$.

Write $q=1/2+k$. The three $s_j$ in $T$ are at least $q-b$ and sum to $3q+b$. The cube root of $T$ is concave in these parameters by Brunn–Minkowski, since all inequalities have affine bounds. On this parameter simplex its minimum is at one of the permutations of $(q-b,q-b,q+3b)$. For $0<u\le v$, one has $T(u,u,v)=u^2v/2$. To compute it, take the box with side lengths $(u,u,v)$ and add the three regions in which exactly one coordinate exceeds its side length. When this excess is $a$, the other coordinates are bounded by their respective lengths minus $a$. Integrating their product gives the stated formula. Thus the singleton volume is at least
$$
\frac{(q-b)^2(q+3b)}2-\frac{8}{3}k^3.
\tag{2}
$$

If $k\le1/10$, (1) implies $b/q\le1/5$. The identity
$(q-b)^2(q+3b)-q^3=b(q^2-5qb+3b^2)$ shows that (2) is at least the equal-offset singleton volume $S(k)=q^3/2-8k^3/3$. If $k\ge1/10$, (2) exceeds $1/12$. For a short rational check, partition $[1/10,1/4]$ at $3/20$ and $1/5$. On an interval $[l,u]$ use $q_0=1/2+l$ and $B=\min\{(64/7)u^2(1-4u/3),3/4-3l,u+1/8\}$. The expression $(q-b)^2(q+3b)$ increases with $q>b$, and at fixed $q$ its derivative in $b$ is $(q-b)(q-9b)$, so its minimum on $[0,B]$ is at an endpoint. The three resulting lower bounds on (2), subtracting $8u^3/3$, are respectively
$$
\frac{507441933}{5359375000},\qquad
\frac{21786321457}{257250000000},\qquad
\frac{779}{6000},
$$
all greater than $1/12$.

## Equality and a direction neighborhood

On $[0,1/10]$, $P(k)$ decreases and $S(k)$ increases. They meet at the root $k_*$ displayed above, so $M\ge\max(P(k),S(k))\ge P(k_*)$. At four equal offsets $k_*$ all six pair cells have volume $P(k_*)$, all four singleton cells have volume $S(k_*)$, and the triple cells have volume $(1/2-k_*)^3/2\le1/16$. The central tetrahedron is smaller still. If $d>0$ at $k=k_*$, the pair lower bound is strict. If $d=0$, the sorted offsets are all equal. This proves the asserted minimum and equality case.

This also excludes a neighborhood of these directions with arbitrary offsets for the unrestricted candidate value. Let the unnormalized normals change by at most $10^{-3}$ each in the sum of the absolute coordinate differences, and keep arbitrary offsets. Compare to the fixed normal case with the same offsets. Each linear score differs by at most $1/2000$ throughout the cube, and the fixed score density is at most $3/4$. A specified sign cell of the fixed arrangement loses at most the sum of four one-sided strips, at most $3/2000$. Thus the actual maximum is at least $M_{\rm reg}-3/2000>(1+2\sqrt2)/48$.
