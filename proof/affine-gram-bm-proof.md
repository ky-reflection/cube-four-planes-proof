# Affine pair bounds for the four-plane partition

The accompanying [rational certificate](affine_gram_bm_certificate.py) verifies a lower bound when every triple of planes is independent and its intersection belongs to the unit cube. We use the tetrahedron notation, the signed offsets, and the deficit identities of the [affine Gram proof](affine-gram-bound-proof.md). In particular, let $M=1/14+\delta$ be a proposed upper bound on every block volume, $t$ the tetrahedron volume, $d_J=M-v_J$ for each nonempty proper exterior sign set $J$, and $e_i=1/2-q_i$ for each exterior cap. Then

$$
\sum_Jd_J=14\delta+t,\qquad e_i=\sum_{J\ni i}d_J-7\delta,\qquad -7\delta\le e_i\le t+7\delta.
$$

For the spherical comparison we use the [four-line radial bound](four-line-radial-proof.md), whose arithmetic is verified by the rational certificate.

The default certificate takes $\delta=428/100000$. It bounds the affine $2+2$ mass from below by $0.454296852068\ldots$ and obtains a contradiction margin greater than $0.00004542349$ over the proposed upper bound $6M$.

## Passing to arbitrary cuts

If some three of the cuts give at most six positive-volume cells, one of those cells has volume at least $1/6$. If they give seven, the [three-cut certificate](seven_scaled_charts_certificate.cpp) gives a cell of volume greater than $4/25$. The fourth plane splits any such cell into at most two, so in either case a final cell exceeds $2/25$.

Otherwise every triple realizes all eight open orthants under the affine map given by its three cut equations. Its convex image of the cube contains the origin in its interior. Indeed, any separating or supporting hyperplane through the origin would exclude the open orthant whose coordinate signs oppose the signs of its normal. The affine map therefore has full rank, and the triple intersection lies in the cube interior. For nonconcurrent cuts the four triple intersections form a tetrahedron. For concurrent cuts, orient the same normals as the facets of a nearby tetrahedron by perturbing the offsets. This gives the same positive relation among the normals, and the deficit estimates below apply directly to the original concurrent cuts with $t=0$.

The maximum cell volume attains its infimum. Normalize each normal to unit length and truncate the offsets to $[-\sqrt3/2,\sqrt3/2]$, which preserves the volume of each sign cell of the centered unit cube. The parameter space is compact and all sign-cell volumes are continuous by dominated convergence off the cutting planes. Thus the pointwise strict lower bound is also strict for $C$.

## A sharper facet estimate

Use the same projection as in the earlier tetrahedron bound. Divide the two projected coefficients by the largest absolute coefficient of the normal. With $a,b\in[0,1]$ and independent uniform $X,Y\in[-1/2,1/2]$, put $R=aX+bY$, $p_0=\Pr(|R|\le 1/2)$, and $U(s)=\Pr(1/2<R\le 1/2+s)$. For $0\le s\le\eta\le 1/2$, symmetric unimodality bounds the lost tail of the shifted band by $s p_0$, and $U(s)\ge(s/\eta)U(\eta)$ by concavity. Hence the mean shifted band has mass at least $p_0(1-\eta/2)+U(\eta)/2$. The facet triangle lies in the symmetric band of mass $p_0+2U(\eta)$ and occupies at most half that mass.

Let $S=a+b-1$. If $S\le0$, then $U(\eta)=0$. Otherwise the two-uniform triangular tail gives

$$
\frac{U(\eta)}{p_0}=\frac{S^2-(S-2\eta)_+^2}{8ab-2S^2}\le\frac{2\eta(1-\eta)}3.
$$

For the inequality use $ab\ge S$ and $S\le1$. If $S\ge2\eta$, the ratio is at most $2\eta(1-\eta/S)/(4-S)$. If $S\le2\eta$, it is at most $\eta/(4-2\eta)$, which satisfies the same bound for $\eta\le1/2$. Both the projected width and this tail ratio increase up to their bounds at $\eta=7H/4$, where $H=t+7\delta$ bounds the absolute offsets. Therefore one may use

$$
\beta(H)=\frac{1+2u}{2-\eta+u},\qquad u=\frac{2\eta(1-\eta)}3
$$

in the inequality $t\le \beta(H)(t+7\delta)$.

## Pair disagreement and cap deficits

For a pair of cuts, let $D$ denote its affine disagreement mass and $D_0$ its centered disagreement mass. Write $p_0=(1-D_0)/2$. First suppose both signed offsets are nonnegative, and vary them along the common multiplier $s\in[0,1]$. The cube root of the positive overlap volume is concave while that overlap is nonempty by Brunn–Minkowski. Its volume derivative at $s=0$ is at most $-(e_i+e_j)/2$. Indeed, each centered section is halved by the other centered plane, and the central marginal density is maximal. The two strips overlap to order $s^2$ for independent normals. Thus, if $e_i+e_j\le6p_0$, the affine overlap is at most $p_0[1-(e_i+e_j)/(6p_0)]^3$. The inequality also holds if the final overlap is empty.

Use the identity $D=1-e_i-e_j-2p$, with $p$ the affine positive overlap. For two negative offsets reflect the cube. For offsets of opposite signs, expand the disagreement in the two sign indicators. The mixed change is nonnegative, so only the two single-offset losses are needed. With $x_+=\max(x,0)$, $x_-=\max(-x,0)$, the resulting lower bound is

$$
D\ge D_0-E((e_i)_++(e_j)_+,1-D_0)-E((e_i)_-+(e_j)_-,1-D_0),\qquad
E(s,q)=\frac{s^2}{3q}-\frac{s^3}{27q^2}.
$$

The range is $0\le s\le3q$. In that range $E$ increases in $s$ and decreases in $q$. For agreement, reverse one cut, which replaces $D,D_0$ by $1-D,1-D_0$ and negates its signed cap deficit. Starting with the elementary recentering estimate $|D-D_0|\le2H$, substitute the pair estimate repeatedly to tighten the upper bounds on centered disagreement and agreement. The spherical angular anchors and the derivative interval $[.23,.51]$ then bound all pair dots in $[-c,u]$.

## Coupling the six pairs through sign-cell deficits

Let $U$ be the uniform centered-disagreement upper bound from the previous step and $T$ the current tetrahedron upper bound. A uniform pair error is $\epsilon=E(2(T+7\delta),1-U)$, so $1-D_0\ge1-D-\epsilon$. For a graph $G$ on the four cuts, the affine pair sum is

$$
\sum_{ij\in G}D_{ij}=8|G|M-\sum_J d_J |\partial_G J|.
$$

On $s/q\in[0,3]$, the polynomial $z^2/3-z^3/27$ is convex. Its perspective $E(s,q)$ is convex and increasing in $s$. Consequently the sum of pair-error upper bounds, using $q=1-D-\epsilon$, is convex in the deficits. The same holds after adding the affine pair sum above. Maximize either expression by placing all deficit mass $14\delta+t$ in a single proper nonempty sign set $J$. After doing so, maximize in $t\in[0,T]$ at an endpoint, again by perspective convexity. The checker verifies the range $s\le3q$ using $2(T+7\delta)\le3(1-8M-\epsilon)$. For the final six-pair error, include the linear credit $-3t$.

For an obtuse angular anchor $(c_0,D_*)$, with $x=-u_i\cdot u_j\le c$, the upper angular derivative gives

$$
D_0\ge D_* -K\big((c_0-c)_++c-x\big),\qquad K=\frac{.51}{\sqrt{1-\max(c_0,u)^2}}.
$$

For $x\le c_0$, integrate the angular derivative bound and use $x\ge-u$. For $x\ge c_0$, monotonicity suffices. Combining this with a graph upper bound for the centered pair sum gives a lower bound on $\sum_{ij\in G}(u_i\cdot u_j+c)$.

Take $G$ first to be each three-edge star. If $g$ is the resulting lower bound, $r=c/(1+c)$, $\alpha_i=r-\lambda_i$, and $\lambda_i\ge\ell$, the positive relation gives $(1+c)\alpha_i=\sum_{j\ne i}\lambda_j(u_i\cdot u_j+c)\ge\ell g$. Thus every $\alpha_i\ge a_0=\ell g/(1+c)$. Write $A=4r-1$. When $\ell\ge r/2$, one also has $\sum_i\alpha_i^2\le A^2/2$. If $a=\max_i\alpha_i\le A/2$, use $\alpha_i^2\le a\alpha_i$. Otherwise use the nonnegative edge weights $\lambda_i\lambda_j(u_i\cdot u_j+c)$. Their degree at vertex $i$ is $(1+c)(r-\alpha_i)\alpha_i$. The degree of the largest-$\alpha$ vertex is at most the sum of the other degrees, which gives $\sum_i\alpha_i^2\le r(A-2a)+2a^2\le A^2/2$ for $A/2\le a\le r/2$.

Set $b_i=\alpha_i-a_0$. Then $\sum_i b_i u_i=4(r-a_0)\bar u$, $\sum_i b_i=A-4a_0$, and $0\le b_i\le r-\ell-a_0$. Bound $\sum b_i^2$ by filling as many entries as possible to this cap. The degree estimate, when applicable, also gives $\sum b_i^2\le A^2/2-2a_0A+4a_0^2$. Taking the squared norm and applying the upper pair dot $u$ now bounds $\|\bar u\|^2$.

For the second graph use $K_4$ minus any one edge $ij$. If its graph lower bound is $g_5$, then $u_i\cdot u_j\le-2+5c+8\|\bar u\|^2-g_5$, because the sum of all six pair dots is $-2+8\|\bar u\|^2$. This refines the upper pair dot in the next iteration.

## Spherical comparison

The centered spherical $2+2$ area fraction is $\sum_{i<j}\arccos(u_i\cdot u_j)/\pi-3$. The checker optionally strengthens the tangent at $x=-1/3$ by a quadratic term $q(x+1/3)^2$. The difference of $\arccos x$ and this quadratic lower candidate vanishes with its derivative at $-1/3$. Its second derivative is decreasing and initially positive. Hence it suffices to check the upper endpoint on the right, as done with rational cosine Taylor bounds. The rank-three Gram matrix gives $\sum_{i<j}(u_i\cdot u_j)^2\ge2/3$ by Cauchy–Schwarz on its nonzero eigenvalues. It follows that $\sum_{i<j}(u_i\cdot u_j+1/3)^2\ge(16/3)\|\bar u\|^2$. This improves the spherical fraction used in the radial cube bound. The affine $2+2$ mass is at least the resulting centered cube mass minus the computed six-pair error after the $3t$ credit. A value greater than $6M$ is the contradiction.

## Reproducing the arithmetic

Compile and run the angular derivative checker in both modes. The pair checker uses the listed triples in `LEGACY_PAIRS` with `--pair`, and those in `SHARP_PAIRS` with `--sharp-pair`. With no arguments it checks the anchor at $c_0=1/3$. The sharp mode uses the derivative upper bound $0.51$. A box of half-width $1/d$ in the quaternion chart changes either normal angle by at most $2\sqrt3/d$, giving an error less than $3.6/d$ in the disagreement of two moved normals.

```sh
c++ -std=c++20 -O3 angular_moment_certificate.cpp -o /tmp/cube-angular
/tmp/cube-angular
/tmp/cube-angular --upper
c++ -std=c++20 -O3 wedge_certificate.cpp -o /tmp/cube-wedge
/tmp/cube-wedge
python3 - <<'PY'
from subprocess import run
from affine_gram_bm_certificate import LEGACY_PAIRS, SHARP_PAIRS

for mode, anchors in (("--pair", LEGACY_PAIRS), ("--sharp-pair", SHARP_PAIRS)):
    for anchor in anchors:
        run(["/tmp/cube-wedge", mode, *map(str, anchor)], check=True)
PY
python3 affine_gram_bm_certificate.py
```

## Why two charts cover the seven-cell case

Orient the three cuts so that the missing open sign cell is $+++$. Write their forms as $f_i(x)=u_i\cdot x-h_i$ on the centered unit cube $K$. Its affine image $(f_1,f_2,f_3)(K)$ misses the open positive orthant. Separation gives a nonnegative vector $\lambda$ with $\sum\lambda_i=1$ such that
$$
g(\lambda)=\max_{x\in K}\sum_i\lambda_i f_i(x)
=\tfrac12\Big\|\sum_i\lambda_i u_i\Big\|_1-\sum_i\lambda_i h_i\le0.
$$
On the boundary of the weight triangle at most two weights are positive. Every pair of positive halfspaces meets the cube interior, so $g>0$ there. The three equations $(\sum\lambda_i u_i)_k=0$ subdivide the triangle into polygons on each of which $g$ is affine. Its minimum occurs at a polygon vertex, necessarily in the interior of the triangle. At least two coordinate equations hold there. Identically zero coordinates already vanish everywhere and can simply be counted as zero. Thus all three weights are positive and, after a cube symmetry, the resultant is parallel to the $z$ axis (it may vanish).

Absorb the weights into the cut equations, then divide all equations by their largest absolute normal coefficient. If this coefficient is horizontal, reflect and permute coordinates to obtain the first row $(1,a,z_1)$ with $0\le a\le1$. The other two $x$ coefficients sum to $-1$ and both have absolute value at most one. Put the larger first, so it lies in $[-1/2,0]$. Reflect $z$ to make the resultant nonnegative. These are exactly the six intervals of the horizontal root. If the largest coefficient is vertical, choose the first row $(a,b,1)$ with $0\le b\le a\le1$ and order the other two rows by descending $z$ component. The other horizontal components in the last row are determined by the zero sums. This gives the vertical root. Ties satisfy at least one of the two descriptions.

For the offsets, under the hypothetical maximum $4/25$, the seven masses sum to one, so every mass is at least $1/25$. Each positive halfspace is a union of three cells; its complement has four. Thus its deficit from $1/2$ lies between $1/50$ and $7/50$. The thin-slab bound translates these into offset bounds. The centered-cone estimate and its derivation from the three pair volumes appear in the opening derivation of the [seven-cell Gram certificate](seven_gram_certificate.py). The two full runs of the [integer verifier](seven_scaled_charts_certificate.cpp) use these inequalities to cover all feasible offsets at each direction box. A run succeeds only after exhausting its stack without an unresolved node.
