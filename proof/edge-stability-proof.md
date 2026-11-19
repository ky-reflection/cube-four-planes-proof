# A quantitative restriction on two-positive cube-edge intervals

Let $m=(1+2\sqrt2)/48$. In a unit cube, suppose four cuts have largest volume at most $m$. Orient them as exterior facets of their interior tetrahedron, with either simultaneous choice of signs at concurrence. Let $L$ be the maximum length of a cube-edge interval on which exactly two cuts are strictly positive and two strictly negative. Take $L=0$ if there is no such interval.

If $L\le 10^{-35}$, the cuts are one of the centered equality configurations of the [edge-separated family](edge-separated-family-proof.md). The numerical cutoff here is deliberately conservative. The proof uses the three-cut seven-cell bound, the affine family interval checker, and the exact local certificate.

## Uniform margins for a configuration with small maximum volume

Use the centered unit cube and unit normals $n_i$, with evaluations $\ell_i(x)=n_i\cdot x-h_i$. The three-cut reduction shows that every triple intersection lies inside the cube. In fact its distance to every face is greater than $d=3/1000$. If its distance to one face were at most $d$, remove the slab between that face and the parallel plane through the intersection. The remaining rectangular box has volume at least $1-d$, and the triple cuts it into at most seven cells. Rescale it affinely to a unit cube and apply the [three-cut bound](seven_scaled_charts_certificate.cpp). A cell of the triple then has volume greater than $(4/25)(1-d)>2m$, which cannot be reduced to pieces of volume at most $m$ by one cut.

The smallest singular value of each three-row matrix of normals exceeds $\sigma=1/6000$. Otherwise, for a unit vector $w$ attaining this singular value, orthogonally project the three normals onto $w^\perp$, keeping their offsets. The projected normals are nonzero and dependent, so again some cell for these three cuts has volume greater than $4/25$. In the original cube each corresponding oriented sign cell can lose at most
$$
\sqrt3\sum_i \frac{\sqrt3}{2}|n_i\cdot w|
\le \frac{3\sqrt3}{2}\sigma.
$$
Here the linear score of any unit normal has density at most $\sqrt3$, by integrating along its largest coordinate, and loss of a specified sign is one-sided. Since $4/25-(3\sqrt3/2)\sigma>2m$, this is a contradiction.

Write $\sum_i\lambda_i\ell_i=-B$, where the weights are positive and sum to one, and $B\ge0$. Every boundary point has at least one evaluation greater than $\mu=1/4000000$ and at least one less than $-\mu$. Indeed, $\{\ell_i\le\mu\text{ for all }i\}$ is a simplex whose vertices solve triples of equations $\ell_i=\mu$. In the evaluation coordinates, this follows from the positive dependence and the fixed weighted sum. The vertices are within $\sqrt3\mu/\sigma<d$ of the original triple points. Thus the entire simplex is inside the cube. The same argument for $\{\ell_i\ge-\mu\text{ for all }i\}$ uses triples $\ell_i=-\mu$. The latter set is empty if $B>\mu$ and is a point if $B=\mu$.

## A vertex close to any cut forces an edge interval

Put $\eta_0=10^{-12}$, $\theta=10^{-17}$, and $a_0=10^{-4}$. For a cut with $|\ell_i(v)|\le\eta_0$ at a vertex $v$, use inward coordinates based at $v$. In these coordinates the sums of its positive and absolute negative normal coefficients are each at least $d-\eta_0$. To see this, evaluate at an interior triple intersection of this cut. All three inward coordinates there belong to $[d,1-d]$, so if the sums are $P,N$, both $dP-(1-d)N$ and $dN-(1-d)P$ are at most $\eta_0$. Use $P+N\ge1$.

If two cuts are within $\eta_0$ of zero at $v$, the other two evaluations are greater than $\mu$ and less than $-\mu$. The two near cuts must have opposite inward normal coefficients, both of absolute value at least $a_0$, in some common coordinate. For otherwise, replace the smaller absolute entry in every opposing coordinate by zero and move both cuts through $v$. The modified normals still have both signs by the coefficient-sum bound and now agree weakly in every coordinate. Each new score differs uniformly from the old by at most $3a_0+\eta_0$, so the sign-change volume per cut is at most $2\sqrt3(3a_0+\eta_0)<1/900$. The two new cap volumes, in either simultaneous orientation, are in $[9/25,16/25]$, using $m<997/12500$ and the original cap bounds $[1-8m,8m]$. The [aligned corner cap bound](edge-separated-family-proof.md#excluding-vertices-on-the-cutting-planes-at-the-candidate-threshold) now gives disagreement at most $16727/38400$ for the modified cuts. The original disagreement would be at most $16727/38400+2/900<11/25$, contrary to its lower bound $1-7m$. Thus an opposite-coefficient coordinate exists. On its edge at inward distances from $\mu/4$ to $\mu/2$, the two near cuts have opposite signs, and the other signs are unchanged. This produces an interval of length $\mu/4$.

If exactly one cut has evaluation of absolute value at most $\eta_0$ and that value is at most $\theta$, choose an inward coefficient of the sign needed to make the total two positives. There is such a coefficient of magnitude at least $(d-\eta_0)/3$. At inward distances from $\eta_0/4$ to $\eta_0/2$, the near cut takes the chosen sign and the other three keep their vertex signs. This interval has length $\eta_0/4$. Consequently, $L\le10^{-35}$ forces every vertex evaluation to have magnitude greater than $\theta$.

## Recovering the vertex sign pattern when the intervals are short

At every vertex there are one or three positive signs, since a two-positive vertex would yield such a sign pattern for a length at least $\theta/2$ along an incident edge. Give a singleton vertex its positive label and a triple vertex its negative label. Along an edge between two different singleton labels, the only two changing cuts must cross at a separation at most $L$: the interval between crossings has either two positives or all negatives, with the latter forbidden. At one crossing the other changing evaluation has magnitude at most $L$ because all normals are unit. This contradicts the boundary's positive margin $\mu$. The same argument with signs reversed applies to distinct triple labels.

An edge between singleton $i$ and triple $i$ also gives an interval of length at least $\theta/2$. At the crossing of cut $i$, there are two other cuts with evaluations greater than $\mu$ and less than $-\mu$. They retain their signs for distance $\mu/2$ on both sides. Cut $i$ crosses from positive to negative and the remaining cut from negative to positive. If the remaining cut is nonpositive at this crossing, use the interval just before it; if it is positive, use the interval just after it. The crossing is at distance greater than $\theta$ from both endpoints.

For an edge between singleton $i$ and triple $j\ne i$, cuts $k,l$ with the other two labels cross at separation at most $L$. Pick one crossing point $x$. There is a point $y$ on the true intersection line of $k,l$ with $\|x-y\|_2\le e=L/\sigma$, by the two-row singular value bound. That bound follows by restricting the quadratic form of the three-row Gram matrix to the two rows. Choose a triple intersection $p$ on the line, with distance greater than $d$ from all faces. On the ray from $p$ to $y$ the line meets the cube boundary at $z=p+t(y-p)$. Since $y$ is within $e$ of a boundary point, $d/(d+e)\le t\le 1+e/(d-e)$. The lower bound follows from all six face inequalities, and the upper bound from a face containing $x$. Also $\|y-p\|_\infty\le1+e$. These bounds give $\|z-x\|_2<1000e<\theta/4$.

Two distinct edges with the same ordered pair $(i,j)$ would therefore yield two distinct boundary endpoints of the same line, since both crossings are at least $\theta$ from every vertex of their edges. Both line endpoints would have $\ell_i>0$, since the evaluation of cut $i$ is greater than $\theta$ everywhere on each edge. The line has an interior triple point with $\ell_i=0$, which is impossible. Thus every ordered pair occurs at most once. The graph argument in the [strict-vertex classification](edge-separated-family-proof.md#the-vertex-sign-pattern-follows-from-the-edge-condition) now applies without any change. The vertices have the prescribed alternating singleton and triple labels, with complementary signs at antipodes.

## Approximating the cuts by the solved family

Scale the cube to $[-1,1]^3$, keeping unit normals and denoting the new offsets by $H_i=2h_i$. Use the parity vertices $v_i$ and off-diagonal coefficients $C_{ij}=-n_i\cdot v_j$ of the family note. Then $C_{ij}>|H_i|+2\theta$. For opposite edges $v_i,-v_j$ and $v_j,-v_i$, the crossing fractions for cut $k$ are $(C_{ki}+H_k)/(C_{ki}+C_{kj})$ and $(C_{kj}+H_k)/(C_{ki}+C_{kj})$. Thus, for complementary rows $k,l$, the differences of their centered fractions and the differences of their normalized offsets each have magnitude at most $L$. The centered fractions belong to $[\theta/2,1-\theta/2]$, since each $C_{ij}\le\sqrt3$. Taking log odds bounds the corresponding disjoint minor of $\log C$ in absolute value by $4L/\theta$.

Subtract row and column constants from $\log C$, setting row 1, column 1 and entry $(2,3)$ to zero off the diagonal. The disjoint minor equations from the family note express every remaining entry as a signed sum of at most three minors. Thus $C_{ij}=r_i c_j\exp(X_{ij})$ for positive factors, with $\sum_j c_j=4$ and $|X_{ij}|\le\eta=12L/\theta$. Define $b=\sum_j(c_j/4)v_j$, $n_i^0=v_i-b$, $n_i'=n_i/r_i$, $H_i'=H_i/r_i$, and $B_i=\sum_{j\ne i}c_j$. At every cube vertex the score difference $|(n_i'-n_i^0)\cdot v|$ is at most $2\eta B_i$, by the off-diagonal bounds and the zero sum of the four parity vertices. Hence it obeys the same bound on the cube.

Choose $k$ with largest $c_k$. There is an index other than $k$ attaining $c_{\min}$, so $|H_k'|\le e^\eta c_{\min}$. For any $i\ne k$ use complementary columns $j,l$, and put $A_i=(C_{ij}+C_{il})/r_i$, with $A_k$ defined in the same way. The normalized offset difference bound says $|H_i'/A_i-H_k'/A_k|\le L$, and $e^{-2\eta}\le A_i/A_k\le e^{2\eta}$. It follows that
$$
|H_i'-H_k'|\le A_i L+ |H_k'||A_i/A_k-1|\le(2L+8\eta)B_i.
$$
Clamp $H_k'$ to $[-c_{\min},c_{\min}]$ to obtain a feasible common family offset $H$; this adds at most $2\eta B_i$ to the bound for each row. The inequalities use $\eta\le 1/10$. The resulting score difference between the actual and approximating planes is at most $14\eta B_i$ everywhere. Since $n_i^0\cdot v_i=B_i$, its largest component has absolute value at least $B_i/3$. Integrating along that coordinate in the side-two cube bounds the sign-change volume by $42\eta$ per cut. Thus each normalized volume in the approximating family differs from its actual counterpart by at most $168\eta<10^{-14}$.

## Returning to the exact local certificate

Set $\varepsilon=10^{-14}$. The approximating family has largest volume at most $m+\varepsilon$. The [rational interval checker](edge_separated_affine_certificate.py) excludes this outside the analytic region $x\le23/50$, $y,z\ge9/10$, $r\le x/5$ in the notation of the family note, with $x\le y\le z\le1$. Here $r=|H|/(4w_{\max})$. Put $x_0=\sqrt2-1$ and $s_0(x)=(5+x)/(48(1+x))$. The largest-weight singleton has volume at least
$$
s_0(x)+(2-y-z)/100+r/6.
$$
In the centered singleton formula the absolute partial derivatives in $y,z$ are at least $1/96$. For the affine term the linear gain is at least $3r/16$, while the cubic loss is at most $r\, x(1+x+y+z)/(300yz)<r/100$, using $r\le x/5$. This proves the displayed bound. If $x\le x_0$, use $s_0(x)-m\ge(x_0-x)/48$. If $x\ge x_0$, the centered pair bound for the split $ad\mid bc$ is $(1-(1-x)^2/[4(1+x)^2])/12$, which grows at slope at least $1/100$ on this interval. The larger affine pair volume is at least the centered value plus
$$
\frac{r S\,(|x-yz|-rS)}{4D},\qquad
S=1+x+y+z,\quad D=(x+y)(x+z)(1+y)(1+z).
$$
Here $S\ge1$, $D\le16$, and $|x-yz|-rS\ge 81/100-23/50-(23/250)(173/50)>64/2500$, so the gain is at least $r/2500$. Combining the pair and singleton inequalities, in both cases
$$
|x-x_0|\le100\varepsilon,\qquad 2-y-z\le1000\varepsilon,\qquad r\le2500\varepsilon.
$$
For the middle inequality when $x\ge x_0$, use $|s_0'(x)|\le1/12$.

The unnormalized weights $(x,y,z,1)$ differ from the candidate weights in sum of absolute differences by at most $1100\varepsilon$. Normalization to sum one at most doubles this, so each coefficient of $b$ differs from its candidate value by at most $2200\varepsilon$. Also $|H|\le4r\le10000\varepsilon$. Returning to the actual planes adds at most $8\eta$ to each normal coefficient, by the parity-vertex tight frame, and at most $48\eta$ to the offsets. Finally rescale the body-diagonal row to have coefficient sum 3 and the other rows to have diagonal component $\sqrt2$, and return the cube side to one. At the candidate the body row has each coefficient greater than 1, and the other rows have diagonal coefficients greater than 1. These rescalings multiply the maximum coefficient error by at most 4 and the scaled-cube offset magnitude by at most 1 for the stated perturbations. Thus all twelve parameters are within $1/2500000$ of the candidate in the [local certificate chart](local_certificate.py). The [strict local minimum](local-candidate-proof.md) there forces the actual planes to be the equality configuration.
