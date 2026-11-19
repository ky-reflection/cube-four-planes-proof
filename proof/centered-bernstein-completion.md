# Completion of the centered alternating vertex chart

Let $m=(1+2\sqrt2)/48$. Assume the four cutting planes pass through the center and are oriented so that at each of the four vertices of one parity there is a distinct, unique positive sign. Then the largest block volume is at least $m$. Equality holds precisely at the centered candidate, up to cube symmetries and relabeling. This covers all six maximal-column choices in this chart, including ties. It leaves other vertex patterns and affine offsets to be treated.

The [first cover](centered-bernstein-proof.md) deals with $1000,1032,1230$. This note closes $1001,1002,1200$ by stronger joint bounds on parameter boxes and an explicit centered neighborhood theorem. Use the [exact face formulas](vertex-chart-overlap-proof.md) and retain the convention that $s_i$ and $p_{ab}$ are multiplied by twelve.

## Joint bounds with a shared affine model

The existing row test uses a common denominator for the three triangle terms. The strengthened test additionally represents each triangle on its two-variable box as an affine form plus a rigorously bounded remainder. The affine coefficients may be proposed by midpoint differentiation. Write the proposed affine form as $L$ in box coordinates $u,v\in[0,1]$. The residual $N-DL$ for the rational term $N/D$ has degree at most $(3,3)$. The exact coefficient formula, after raising the numerator and denominator to that degree, uses the four terms with indices $(i-a,j-b)$, $a,b\in\{0,1\}$, weighted by
$$
\frac{(a i+(1-a)(3-i))(b j+(1-b)(3-j))}{9}.
$$
In each term use $N_{i-a,j-b}-D_{i-a,j-b} L(a,b)$. Using the same weights on the denominator gives an enclosing interval for the remainder by the minimum and maximum ratios of the nonzero-denominator coefficients. A zero denominator coefficient again has zero numerator. This is the Bernstein product rule for degrees two and one, using only integer operations on the proposed coefficients.

The overlap model can use a constant enclosure or another affine form. On a positive box the two branches can be written together as
$$
o=\frac{(\delta_+)^2}{d_+}+\frac{(({-\delta})_+)^2}{d_-}.
$$
The first derivatives match at $\delta=0$. The code encloses the derivatives of each rational branch by interval quotient differentiation. On the positive support of a branch the denominator is at least one. Outside its positive support the value and first derivatives are zero. For each branch this gives an enclosing interval throughout the box. Let $[g_j^-,g_j^+]$ be the resulting bound for derivative $j$ in normalized coordinates $[-1,1]^8$. If the proposed affine coefficient is $g_j$, the mean value estimate for the remainder adds at most $\max(|g_j^--g_j|,|g_j^+-g_j|)$, summed over the coordinates. The center value is also enclosed. Every bound used to accept an overlap model is outward rounded with integers.

On the same normalized coordinates, each weighted sum of the seven volume expressions has an affine slope vector. Its minimum affine contribution is minus the sum of the absolute slopes. A linear program proposes the nonnegative weights, and integer arithmetic checks the resulting inequality at the candidate threshold. The check combines both the slopes and the signed remainder coefficients before taking their bounds. This is important because different cells can contain the same triangle or overlap with opposite coefficients. A failure to find a weight only leaves the box unresolved. The interval overlap test is enabled only on boxes with every entry at least $1/8$; elsewhere the original enclosures remain available.

## An explicit neighborhood of the candidate

Relabel so the small singleton at the candidate has index zero. Put $x=\sqrt2-1$, $b=1/(3+x)$, and $w=(xb,b,b,b)$. After rescaling the rows, write the entries as $C_{ij}=w_j(1+\eta_{ij})$ with every row sum of $\eta$ zero. Suppose the pivot-normalized entries are within $1/200$ in sup norm of the corresponding candidate. Before this rescaling, the relative errors have magnitude at most $r/x$ with $r=1/200$, and their squared norm is at most $r^2(3/x^2+5)$. There are three free entries in the small column and five other free entries. The fixed entry in each row has error zero. Subtracting the row mean contracts its Euclidean norm and leaves each coordinate bounded by $r/x$. Dividing by one plus the row mean costs at most a factor $(1-2r/(3x))^{-1}$. Hence
$$
\|\eta\|_\infty\le e=1/80,\qquad \|\eta\|\le1/40,
$$
where unadorned norms are Euclidean.

Decompose $\eta_{ij}=X_{ij}+y_j+y_i/3$ with $\sum_j y_j=0$ and $X$ having zero row and column sums. This is orthogonal, with $\|\eta\|^2=\|X\|^2+(8/3)\|y\|^2$. For the complementary rows $k,l$ of an edge $ij$, the signed linear minor is $d_{ij}=X_{ki}+X_{lj}-X_{kj}-X_{li}$. The [minor estimate](vertex-chart-overlap-proof.md#the-logarithmic-split-directions) gives $\sum d_{ij}^2\ge4\|X\|^2$; it is an inequality on the matrix of linear differences and applies in these relative coordinates as well.

Let $W=2x/(1+x)$, $W_0=1$, $W_i=W$ for $i>0$, and $\mathcal F=\sum_i W_i s_i$. In each row the weighted triangle sum has cubic numerator $N=\sum_j W_j C_{ij}^2(C_{ik}+C_{il})$ and cubic denominator $K$ from the first cover. For its change from the candidate, $N-\mathcal F_{i,*}K$ has zero constant and linear terms, where $\mathcal F_{i,*}$ is that row's candidate value. The quadratic term divided by $K_*$ is at most $(129/2000)\|\eta_i\|^2$ by the [row quadratic formulas](vertex-chart-overlap-proof.md#a-local-bound-from-the-singleton-caps). The [exact checker](centered_candidate_neighborhood.py) also verifies those formulas on a basis and its sum.

To bound the cubic term, each numerator monomial contains a repeated entry. Its relative cubic perturbation $\eta_j^2\eta_k$ is bounded by $e\eta_j^2$, and $W_j t_{ij,*}<0.21$. The analogous denominator coefficient for a repeated entry is $t_{ij,*}<0.354$. The only remaining denominator monomial is twice the product of the three entries. Its fraction of $K_*$ is at most $1/4$ by the arithmetic-geometric mean inequality, and the absolute product of the three relative perturbations is at most $e\|\eta_i\|^2/3$. Since $\mathcal F_{i,*}\le1$, the total cubic error divided by $K_*$ is at most $0.65 e\|\eta_i\|^2$. Also $K\ge(1-e)^3K_*$.

For the overlap terms the disjoint minor is $w_iw_j(d_{ij}+\epsilon_{ij})$, where $\epsilon_{ij}=\eta_{ki}\eta_{lj}-\eta_{kj}\eta_{li}$. Each squared entry contributes to two such error bounds, so $\sum\epsilon_{ij}^2\le2e^2\|\eta\|^2$, by bounding the square of the difference of products by twice the sum of the squared products. The overlap denominator has positive degree-four terms, so it is at most $(1+e)^4$ times its candidate value $(w_i+w_j)^3$. For both edge types $(W_i+W_j)w_i^2w_j^2/(w_i+w_j)^3>0.028$. Using $(d+\epsilon)^2\ge(15/16)d^2-15\epsilon^2$ and the preceding row bounds gives
$$
\mathcal F-\mathcal F_*\le 0.076\|\eta\|^2-0.099\|X\|^2. \tag{1}
$$

We now use the high singleton bounds. Let $z_j=s_j-12m$ for $j=1,2,3$, and write $z=D\eta+R$ using the derivative at the candidate. One share in a triangle has the exact form $f=f_0(1+\eta_A)/(1+v)$, where $v=f_0\eta_A+(1-f_0)\eta_B$. Its difference from $f_0$ is $f_0(1-f_0)(\eta_A-\eta_B)/(1+v)$. The first-order remainder is obtained by multiplying that difference by $-v$. Thus the difference has magnitude at most $0.36\|\eta_i\|$ and the first-order error at most $0.36\|\eta_i\|^2$ in row $i$. Both reference shares in a triangle are at most $3/4$, so the triangle first-order error is at most $0.7\|\eta_i\|^2$. The same exact share differences give an overlap upper bound of $0.26(\|\eta_k\|^2+\|\eta_l\|^2)$, using the two complementary rows and the overlap denominator at least one. Each other row enters twice among overlaps incident to a given singleton. Therefore $\|R\|\le2.2\|\eta\|^2$.

For clarity, the linear map has a particularly useful split. Put
$$
A=\frac{1+2x}{2(1+x)^2},\quad B=\frac{1-x}{4(1+x)^2},\quad
G=\frac38+\frac{5A-B}3,\quad K=\frac{-8A+16B}3.
$$
The residual part sends $X$ to $((3/8-A)X_{0j}-B X_{j0})_{j>0}$, a vector of sum zero. Its norm is at most $d_0\|X\|$ with $d_0^2=(3/8-A)^2+B^2$. The edge-weight part sends $y$ to the vector with entries $G y_j+(1/8-A/3+5B/3)y_0$. It has singular values $G,G,|K|/2$ with respect to the norm of $y$, by separating the zero-sum high coordinates from their common-value direction. Thus the exact constants checked above give $\|D\eta\|\le0.7\|\eta\|$, and, using the fact that the residual part has sum zero,
$$
\sqrt{8/3}\,\|y\|\le \tfrac16\|X\|+4\|z\|+9\|\eta\|^2.
$$
Since $\|\eta\|\le1/40$, the first and last terms together are at most $(47/120)\|\eta\|$. Squaring with $(a+b)^2\le(6/5)a^2+6b^2$ and using the orthogonal decomposition yields
$$
\|\eta\|^2\le (5/4)\|X\|^2+120\|z\|^2,\qquad \|z\|<0.019.
$$

If all the high singletons are at most $12m$, then $z_j\le0$. The change in the total singleton volume is $\mathcal F-\mathcal F_*+(1-W)\sum z_j$. By (1) and $\sum z_j\le-\|z\|$, this is strictly negative unless both $X=0$ and $z=0$. In the latter case the preceding norm bound gives $\eta=0$. The seven representative volumes sum to six. Thus a strict decrease of the singleton total forces one of the three pair volumes to exceed $12m$, proving the neighborhood theorem.

## Independent affine offsets in the chart with two reciprocal pairs

Use the cube $[-1,1]^3$ with normalized volumes, and write the four cuts as $n_i\cdot x=H_i$, with the normals still in the alternating vertex chart. Suppose the pivot word is $1032$. The four offsets may vary independently subject to
$$
|H_i|\le \eta\min_{j\ne i} C_{ij},\qquad \eta=3/25.
$$
Then the largest block is strictly greater than $m$. This proof works directly with the four two-positive sign sets $02,13,03,12$.

At scale $\rho\in(\eta,1]$, the face coordinates of the smaller concentric cube retain the alternating vertex signs, since the effective offsets are $H_i/\rho$. On each face with positive parity corners $v_a,v_b$, the sign set $ab$ occupies the portion outside the union of its four corner triangles. Its area is therefore at least the full face area minus the sum of the triangles. Apply this to the four faces with the selected pairs and discard any contributions from the other two faces. On opposite faces the triangle sizes associated with one entry $C_{ij}$ have numerators $(C_{ij}+H_i/\rho)^2$ and $(C_{ij}-H_i/\rho)^2$ over the same denominator $12(C_{ij}+C_{ik})(C_{ij}+C_{il})$, in the face-pyramid volume normalization. In this pair selection every row contributes its two nonpivot entries, each on two opposite faces. The centered triangle terms in each row sum to at most $1/2$ in the scale multiplied by twelve, by the row inequality in the [first cover](centered-bernstein-proof.md#the-three-covered-charts). The terms linear in the offsets cancel. The additional quadratic loss is at most $\eta^2/(3\rho^2)$ after summing over the four rows, since $|H_i|\le\eta C_{ij}$. Thus the four selected cells have a radial face-pyramid contribution of at least $1/3-\eta^2/(3\rho^2)$ at this scale. Integrating over the outer shell with weight $3\rho^2\,d\rho$ gives their actual combined volume at least
$$
\frac13-\eta^2+\frac23\eta^3>4m.
$$
The strict comparison at $\eta=3/25$ is verified by the exact neighborhood checker. This includes nonconcurrent cuts and imposes no relation among their four offsets.

## Reproduction

From the directory of this note, run

```sh
c++ -O3 -std=c++20 centered_chart_cover.cpp -o /tmp/centered-chart-cover
python3 -B centered_candidate_neighborhood.py
/tmp/centered-chart-cover 1001 88 5000000 16 24
/tmp/centered-chart-cover 1002 88 5000000 16 24
/tmp/centered-chart-cover 1200 88 7000000 16 24
```

The outputs have zero pending boxes, with the following handoffs to the neighborhood theorem.

| Chart | Visited | Excluded by strengthened joint bounds | Unresolved | Max distance to a candidate |
| --- | ---: | ---: | ---: | ---: |
| 1001 | 91,863 | 12,560 | 179 | 112/32768 |
| 1002 | 270,257 | 37,420 | 519 | 128/32768 |
| 1200 | 1,521,057 | 177,229 | 2,715 | 112/32768 |

The box program returns a nonzero status for these handoffs because it verifies only the boxes it excludes. Each distance in the table is less than $1/200$, so the local theorem covers every remaining box. The possible candidate points in a chart are exactly the choices of small label omitted from the pivot word, with all nonpivot entries equal to one except for the entries in the small column, which equal $\sqrt2-1$. The distance calculation encloses that irrational entry on the integer grid and uses the farthest endpoint in each coordinate.
