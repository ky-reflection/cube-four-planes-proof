# Affine cuts with all three-plane intersections in the cube

Let $Q=[-1/2,1/2]^3$. Suppose each three of the four planes is independent and its intersection point belongs to $Q$. The largest part of $Q$ then has volume strictly greater than
$$
M_0=\frac1{14}+\delta,\qquad \delta=\frac{283}{100000}.
$$
The arithmetic in the proof is checked by [affine_gram_bound.py](affine_gram_bound.py).

## Signed offsets and the tetrahedron

Write the exterior sides as $u_i\cdot x>h_i$, with $u_i$ unit and $\sum_i\lambda_i u_i=0$, $\lambda_i>0$, $\sum_i\lambda_i=1$. In the nonconcurrent case the exterior sides point away from the tetrahedron $T$ spanned by the four intersection points. Put $t=\operatorname{vol}(T)$, $B=\sum_i\lambda_i h_i$, so $B>0$. In the concurrent case choose signs from the positive relation, and put $t=B=0$. The all-exterior pattern is absent in both cases.

Assume all pieces have volume at most $M_0$. Let $q_i$ be the volume of exterior cap $i$, and $e_i=1/2-q_i$. For each nonempty proper exterior set $J\subset\{1,2,3,4\}$, let $d_J=M_0-v_J$, with $v_J$ the volume of the corresponding sign cell, including zero for a missing cell. Then
$$
\sum_J d_J=14\delta+t,\qquad e_i=\sum_{J\ni i}d_J-7\delta,\qquad
-7\delta\le e_i\le t+7\delta.
$$
Write $P=\sum_i(e_i)_+$, $N=\sum_i(-e_i)_+$, and $k=\#\{i:e_i<0\}$. At most three signed offsets are negative, since $B\ge0$. Counting the incidence of the nonnegative coordinates in each $J$ gives
$$
P\le\min(3,4-k)(14\delta+t)-7\delta(4-k),\qquad N\le7k\delta.
$$
In particular
$$
P\le3t+21\delta,\quad P+N\le3t+28\delta,\quad
P^2+N^2\le(3t+21\delta)^2+(7\delta)^2.
\tag{1}
$$

The centered slab between a facet of $T$ and its antipodal plane has volume $2|e_i|$. It contains the convex hull of the facet triangle and its reflection. At relative height $s$, the area of this hull is $A_i[1+2s(1-s)]$, where $A_i$ is the facet area. Thus $2|e_i|\ge(8/3)|h_i|A_i$. The signed tetrahedron-volume identity is $t=\sum_i h_i A_i/3$, so $\sum_i|e_i|\ge4t$, and (1) implies $t\le28\delta$. This also holds when $t=0$.

The theorem of [Barthe and Koldobsky](https://doi.org/10.1016/S0001-8708(02)00055-5) gives slab volume at least $2s$ for centered half-width $s\le3/8$. Here $|e_i|\le t+7\delta<3/8$, so monotonicity first gives $|h_i|<3/8$, then $|h_i|\le |e_i|$.

There is a stronger bound on $t$ when these offsets are small. Project the facet plane along a coordinate where $|u_{ij}|=c=\max_j|u_{ij}|\ge1/\sqrt3>4/7$. In the projected square its cube section is $|R-h_i|\le c/2$, where $R$ is the sum of the other two scaled uniform coordinates and has a symmetric unimodal density. For $h=h_i>0$, this section is contained in the symmetric region $|R|\le c/2+h$. Every triangle in a centrally symmetric planar region has at most half the area. Indeed, if the center is outside the triangle, the triangle and its reflection have disjoint interiors. If the center is interior, the triangle's vertices $v_1,v_2,v_3$ and their antipodes form a simple cyclic hexagon $v_1,-v_3,v_2,-v_1,v_3,-v_2$ of twice its area by the shoelace formula. That hexagon is inside the region.

For $0\le s\le h\le c/2$, unimodality gives
$$
\Pr(|R|\le c/2+h)\le(1+2h/c)\Pr(|R|\le c/2),\qquad
\Pr(|R-s|\le c/2)\ge(1-s/c)\Pr(|R|\le c/2).
$$
The second inequality bounds the lost interval $[-c/2,-c/2+s)$ by $s/c$ of the original mass. If both projected coefficients vanish, the section is a square throughout this interval and the same estimates hold. Integrating $s$ gives $h A_i\le\beta(H)e_i$ for $h_i>0$, where $H=\max |h_i|$ and
$$
\beta(H)=\frac{1+7H/2}{2(1-7H/8)}.
$$
The projected-area factor $1/c$ cancels in the ratio. Discard the negative terms of the signed volume identity and apply (1). Starting at $t\le28\delta$, the bound $t\le\beta(t+7\delta)(t+7\delta)$ can be iterated with rational upper rounding.

## Refining the Gram matrix

Let $D_{ij}(h)$ be the affine sign disagreement probability and $D_{ij}(0)$ the probability for the same normals through the center. At most eight sign cells contribute to disagreement and at most seven contribute to agreement, since the all-exterior pattern is absent. Therefore
$$
1-7M_0\le D_{ij}(h)\le8M_0.
\tag{2}
$$
Recentring changes this probability by at most $2H$, since each individual change has volume $|e_i|\le t+7\delta$. In the checker $H$ stands for the upper bound $t+7\delta$.

For centered planes use the lower angular derivative $23/100$ for a pair in a fixed normal span, checked by [angular_moment_certificate.cpp](angular_moment_certificate.cpp). The three anchor points $(c_0,D_0)$ have the form $D\ge D_0$ at angle $\theta_*=\arccos(-c_0)$. They are
$$
(c_0,D_0)=(1/3,143/250),\quad
\left(\frac{200^2-131^2}{200^2+131^2},\frac{596}{1000}\right),\quad
\left(\frac{100^2-61^2}{100^2+61^2},\frac{616}{1000}\right).
$$
They are verified by [wedge_certificate.cpp](wedge_certificate.cpp) using its default arguments, the option **--pair 200 131 404**, and the option **--pair 100 61 384**. The latter two commands bound disagreement for the complementary acute angle, so subtract their bounds from one. For $\rho=u_i\cdot u_j\le-c$, $c\ge c_0$, convexity of the derivative of $\arcsin$ gives
$$
\theta-\theta_*\ge
\frac{c-c_0}{\sqrt{1-((c+c_0)/2)^2}}.
$$
An upper bound for $D_{ij}(0)$ therefore excludes overly negative $\rho$ by a rational squared inequality. An upper bound for $1-D_{ij}(0)$ bounds overly positive $\rho$ by applying the same argument to the reversed normal. In combination with (2) this starts bounds $-c\le\rho\le u$, using the best estimate from the three anchors.

Write $c_*=\max(c,u)<1$, $A=2+\sqrt3 c_*/\sqrt{1-c_*^2}$, $L=\sqrt3/\sqrt{1-c_*^2}$. The separated affine error bound is
$$
|D_{ij}(h)-D_{ij}(0)|
\le A(h_i^2+h_j^2)+2L|h_i h_j|.
\tag{3}
$$
For the single-offset term, project the section $u_i\cdot x=s$ orthogonally back to the central plane. The two projected sections at $s,-s$ lie between $(1-2c_i|s|)Q_i(0)$ and $(1+2c_i|s|)Q_i(0)$, $c_i=\max_j|u_{ij}|$. Their symmetric difference has area at most $8|s|$, since $c_i\operatorname{area}(Q_i(0))\le1$. On their intersection, the other centered inequality changes only in a strip of area at most $2\sqrt3|s|\,|\rho|/\sqrt{1-\rho^2}$, using the cube diameter. Central symmetry relates the signed slice integrals at $s,-s$. Taking half their difference and integrating $s$ bounds the single-offset term by $(2+\sqrt3|\rho|/\sqrt{1-\rho^2})h_i^2$. The two-offset term is $-2\mathbb E[\Delta_i\Delta_j]$, where $\Delta_i$ is the change of the $i$-th sign indicator. Its support is the intersection of two one-sided slabs, bounded by $\sqrt3 |h_i h_j|/\sqrt{1-\rho^2}$ by coarea and cube diameter. In a lower bound for $D_{ij}(h)$, this term can be dropped when $h_i,h_j$ have opposite signs.

The bound (3) now sharpens $c,u$ in (2). The positive normal relation gives the further weight and mean estimates
$$
\ell\le\lambda_i\le r,\qquad r=\frac{c}{1+c},\quad
\ell=\max\!\left(0,1-3r,\frac{\sqrt{(1-2c)/3}}{1+\sqrt{(1-2c)/3}}\right)
$$
with the square-root term used for $c<1/2$. The upper bound comes from taking the dot product with $u_i$. The square-root lower bound follows from
$$
\lambda_i^2=\left\|\sum_{j\ne i}\lambda_j u_j\right\|^2
\ge(1+c)\sum_{j\ne i}\lambda_j^2-c(1-\lambda_i)^2
\ge(1-2c)(1-\lambda_i)^2/3 .
$$
For $\bar u=\sum_i u_i/4$, write $4r\bar u=\sum_i\alpha_i u_i$, $\alpha_i=r-\lambda_i$. Then
$$
16r^2\|\bar u\|^2
\le (1-u)\sum_i\alpha_i^2+u(4r-1)^2.
$$
The maximum of $\sum_i\alpha_i^2$ with $0\le\alpha_i\le r-\ell$ and total $4r-1$ is obtained by filling as many entries as possible to $r-\ell$ and putting the remainder in one entry. Also $\sum_{i<j}(\rho_{ij}+1/3)=8\|\bar u\|^2$, so each $\rho_{ij}\le-1/3+5(c-1/3)+8\|\bar u\|^2$.

For three distinct indices $i,j,k$, the exact volume relation is
$$
t=\frac{B^3}{6\lambda_i\lambda_j\lambda_k|\det(u_i,u_j,u_k)|}.
$$
For $B>0$, the coordinates $\lambda_i(h_i-u_i\cdot x)/B$ are barycentric coordinates, giving this identity by their Jacobian. At concurrence both sides are zero. Using (1) gives $B\le r(3t+21\delta)$. The Gram determinant gives the lower bound $|\det|\ge \sqrt{1-3c_*^2-2c_*^3}$ when the radicand is positive. A lower bound for the triple weight product is found with each weight in $[\ell,r]$ and the triple sum at least $1-r$. When $3\ell\le1-r$, one may first decrease weights to give sum $1-r$, then minimize the concave log product over that polygon. A minimizer is a vertex, with at least two weights at endpoints. These substitutions shrink $t$ and hence $H$ for another iteration. The exact checker obtains
$$
t\le0.000111,\qquad -0.406609\le\rho_{ij}\le0.115924.
$$

## The final volume comparison

Put $\theta_0=\arccos(-1/3)$. The tangent to $\arccos$ at $-1/3$ is a lower bound throughout $[-1,1/3]$. This follows from convexity up to zero and the fact that the derivative on $[0,1/3]$ is at least the tangent slope. Consequently
$$
\sum_{i<j}\arccos\rho_{ij}\ge6\theta_0-6\sqrt2\|\bar u\|^2.
$$
On the unit sphere, pair disagreement is angle divided by $\pi$. Positive dependence excludes the all-same signs, and the pointwise number of differing pairs is three or four. Thus the spherical fraction $\mu$ of the centered $2+2$ directions is at least $(6\theta_0-6\sqrt2\|\bar u\|^2)/\pi-3$.

For a spherical set of fraction at least $\mu$, the least cube volume of its cone is achieved by taking the directions where the cube radial function is smallest. For the value of $\mu$ in the checker, these form the six disjoint axis caps $\max_j|\omega_j|\ge1-\mu/3>1/\sqrt2$. Integrating the cube radial function $(2\max_j|\omega_j|)^{-1}$ gives a lower bound
$$
V_2(0)\ge \frac{\pi}{4}\big[(1-\mu/3)^{-2}-1\big].
$$

For the affine arrangement, $\sum_{i<j}D_{ij}(h)=3(1-t)+V_2(h)$. Sum the lower side of (3), keeping the mixed loss only for pairs of equal-sign offsets. Since $|h_i|\le|e_i|$, (1) bounds the loss by
$$
(3A-L)(t+7\delta)(3t+28\delta)
+L[(3t+21\delta)^2+(7\delta)^2].
$$
After crediting the $+3t$ in $V_2(h)=\sum D_{ij}(h)-3+3t$, the bound on the loss is convex in $t$, with maximum at an endpoint of $[0,0.000111]$. The checker proves $V_2(h)>0.44607413$. Six $2+2$ pieces of size at most $M_0$ would have total less than $0.44555143$, a contradiction.

Any 15-piece cut has all four three-plane intersections inside $Q$. The same is true for a 14-piece cut if its two absent sign patterns differ in at least two coordinates. In those cases every three cuts exhibit all eight signs. The affine map given by their three plane equations must have full rank and its convex image contains zero in its interior, since a separating plane would otherwise exclude an open orthant. Thus the bound applies to those 14-piece cuts as well. A 14-piece cut with two absent patterns differing in just one coordinate is a separate seven-piece problem after deleting that coordinate.
