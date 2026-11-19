# A continuous weight cover for the four-plane problem

The [rational certificate](affine_weight_graph_certificate.py) proves that if every triple of four cutting planes is independent and intersects inside the unit cube, the maximum block volume is strictly greater than

$$
M=\frac1{14}+\frac{44}{10000}=0.075828571428\ldots.
$$

The result combines with the [seven-block certificate](seven_scaled_charts_certificate.cpp) and the counting argument in the [affine Gram proof](affine-gram-bound-proof.md). It applies to every four-plane cut. The argument uses the notation $t,c,u$ from the [height and Gram bound](affine-frame-height-proof.md). Here $-c\le\rho_{ij}\le u$ for every pair dot $\rho_{ij}=u_i\cdot u_j$, and $t$ is the volume of the central tetrahedron.

## A lower envelope for centered pair disagreement

Let $x=-\rho$ for a centered pair of planes, and let $D_0(x)$ be its disagreement mass. This notation also depends on the pair's orientation in the cube. The argument below bounds it uniformly over those orientations. A verified pair anchor $(x_0,D_*)$ is a bound $D_0(x_0)\ge D_*$ for all orientations. The angular derivative bounds give

$$
\frac{23}{100}\le \frac{d}{dx}D_0(x)
\le \frac{51}{100\sqrt{1-x^2}}
\le \frac{572}{1000}
$$

on the interval used here, including the adjacent anchors. In the first two inequalities one varies the angle while keeping the pair in a fixed plane. The last inequality is a rational squared comparison since all anchor points used in the interpolation have $|x|\le0.445$. Thus between two consecutive anchors use the larger of the lower lines from the left anchor with slope $0.23$ and from the right anchor with slope $0.572$. To get a globally valid affine bound $D_0(x)\ge a+b x$, subtract $bx$ and minimize that maximum on each interval. It suffices to evaluate the endpoints and, when inside, the intersection of the two lines. The resulting calculation proves, for example,

$$
D_0(x)\ge0.452+0.348x,\qquad -0.12\le x\le0.445.
$$

The checker computes several such lines with rational slopes. Each provides a valid bound in the full dot domain of this certificate.

## Coupling all six pair dots

Let $G$ be any nonempty graph on the four normals. The capped-deficit calculation from the height and Gram bound gives an upper bound $U_G$ on the sum of the centered disagreements for edges of $G$. For each lower affine line $D_0(x)\ge a+b x$, this implies

$$
\sum_{ij\in G}(\rho_{ij}+c)\ge\frac{|G|(a+bc)-U_G}{b}.
$$

Use the maximum of these lower bounds and the earlier single-anchor graph lower bound. The checker uses all 63 nonempty graphs. It also uses the star graph lower bound to tighten the weight intervals. If $g$ bounds a star sum from below, $r=c/(1+c)$, and $\ell$ bounds every weight from below, the positive relation implies $(1+c)(r-\lambda_i)\ge\ell g$. The weight sum then raises $\ell$ if the new upper weight bound is strong enough.

Write $z_{ij}=\lambda_i\lambda_j(\rho_{ij}+c)\ge0$, with $\lambda_i$ the positive normalized weights in the normal relation. Their weighted degrees are

$$
b_i=\sum_{j\ne i}z_{ij}=\lambda_i[c-(1+c)\lambda_i].
$$

Index the four normals from 0 to 3. Set $x=z_{01},y=z_{02}$ and $d=(b_0+b_1+b_2-b_3)/2$. The remaining four edges are

$$
z_{03}=b_0-x-y,\quad z_{12}=d-x-y,\quad z_{13}=b_1-d+y,\quad z_{23}=b_2-d+x.
$$

This expresses every edge as a constant plus a fixed linear form in two variables. The constants depend on the weights. For a continuous box of weights, compute rational lower and upper bounds $k_{ij}^-\le k_{ij}\le k_{ij}^+$ on these constants and a rational upper bound $L_{ij}$ on $1/(\lambda_i\lambda_j)$. The inequalities $0\le z_{ij}\le(c+u)\lambda_i\lambda_j$ yield linear constraints in $x,y$ after outward rounding. The graph lower bound $g_G$ yields another necessary inequality

$$
\sum_{ij\in G}L_{ij}(k_{ij}^++s_{ij}x+t_{ij}y)\ge g_G.
$$

Indeed, the true edges are nonnegative and $L_{ij}$ and $k_{ij}^+$ are upper bounds. The checker adds the same type of inequality for the hypothesized lower bound $m_2$ on $\|\bar u\|^2$, since

$$
\sum_{i<j}(\rho_{ij}+c)=8\|\bar u\|^2-2+6c.
$$

This gives a two-variable linear feasibility problem for each weight box. The checker encloses the coefficients with integers on the grid $2^{-22}$. It clips the initial rectangle in $(x,y)$ by every half-plane using exact intersection arithmetic. A vertex stores the constraint defining its outgoing edge, so new vertices are intersections of two original constraint lines. The weight boxes are bisected whenever the polygon is nonempty. Sorting the weights and imposing their sum equal to one narrows boxes while preserving the full parameter domain. Every terminal box is certified empty or eliminated by those weight constraints. Consequently the default run proves $\|\bar u\|^2<0.011$.

## Final comparison

Use this mean bound in the angular quadratic estimate from the [affine pair proof](affine-gram-bm-proof.md), and then use the four-line radial correction. The centered $2+2$ cube volume exceeds $0.4690$, and the affine $2+2$ mass exceeds $0.45530$. This is incompatible with an upper bound $6M<0.45498$ on its six pieces. The checker carries out the comparison with fractions.

## Reproduction

Run the angular, pair and height checks listed in the [height and Gram proof](affine-frame-height-proof.md). The additional pair anchors use the `--sharp-pair-any-angle` mode. It takes the tangent of half the complementary angle, so it also permits a tangent greater than one. Run the following from the proof directory.

```sh
c++ -std=c++20 -O3 wedge_certificate.cpp -o /tmp/cube-wedge
python3 - <<'PY'
from subprocess import run
from affine_weight_graph_certificate import PAIR_LINE_ANCHORS

for sine, upper in PAIR_LINE_ANCHORS:
    run(["/tmp/cube-wedge", "--sharp-pair-any-angle", "100", str(sine), str(upper)], check=True)
PY
python3 affine_weight_graph_certificate.py
```
