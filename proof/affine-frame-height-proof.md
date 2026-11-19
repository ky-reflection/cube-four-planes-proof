# Height and Gram bounds for the central tetrahedron

Suppose each triple of the four cutting planes is independent and its intersection lies in the centered unit cube. Use the outward unit normals $u_i$, signed offsets $h_i$, positive relation $\sum_i\lambda_i u_i=0$ normalized by $\sum_i\lambda_i=1$, and tetrahedron volume $t$ as in the [affine Gram proof](affine-gram-bound-proof.md). The [rational checker](affine_frame_height_certificate.py) proves that if all blocks have volume at most $M=1/14+0.0062$, then

$$
t\le0.001818,\qquad |h_i|\le0.045218,\qquad -0.490800\le u_i\cdot u_j\le0.290961.
$$

This is a restriction on the parameters of a possible cut with that maximum block volume. The parameter restriction leaves cases to check before deriving a contradiction at this threshold.

The `--candidate` mode also covers the proposed optimum $m=(1+2\sqrt2)/48$ by using the rational upper bound $0.079759>m$. Under the assumption that all blocks have volume at most $m$, it gives

$$
|h_i|\le0.138072,\qquad -0.740529\le u_i\cdot u_j\le0.385202,\qquad
\lambda_{\min}(N^\mathsf{T}N)>0.049451.
$$

Here $N$ is the matrix with the four unit normals as rows. These conditions apply to 15-block cuts and to 14-block cuts whose two absent sign patterns differ in at least two coordinates. For the other 14-block cuts, deleting the coordinate in which the patterns differ gives the seven-block case.

## Height bound and the four minors

Write $B=\sum_i\lambda_i h_i\ge0$, $a=\min_i\lambda_i$, and $H=t+7\delta$ for an upper bound on the absolute offsets, where $\delta=M-1/14$. Let $v_i$ be the intersection of the three facets other than $i$. The positive relation gives

$$
B/\lambda_i=h_i-u_i\cdot v_i\le H+\sqrt3/2<H+7/8.
$$

The equality follows by evaluating the relation at $v_i$. This gives $B\le (H+7/8)a$, which offsets the small denominator when a weight approaches zero.

Let $N$ be the matrix with the unit normals as rows. Its four signed minors form a vector in the one-dimensional kernel of $N^\mathsf{T}$, so their absolute values are $K\lambda_i$ for a common positive $K$. By Cauchy–Binet,

$$
K^2=\frac{\det(N^\mathsf{T}N)}{\sum_i\lambda_i^2},\qquad
 t=\frac{B^3}{6K\prod_i\lambda_i}.
$$

The second formula is the usual barycentric volume formula, substituting the corresponding minor. It holds by continuity at concurrence, where $t=B=0$. The checker partitions the allowed range of $a$ into 64 intervals. On $[a_-,a_+]$, the three remaining weights lie between $a_-$ and the available upper weight bound, and their sum is at least $1-a_+$. Lower their sum to the smallest value allowed by these constraints and minimize the product at a vertex, by concavity of the logarithm. If $B_*$ is any upper bound on $B$ in the interval, use

$$
 B^3/a\le\min\big((H+7/8)^3a_+^2,\ B_*^3/a_-\big),
$$

omitting the second expression when $a_-=0$. This bounds $t$ in each weight interval. We also retain the bound obtained from the triple-minor formula in the earlier proof.

## Bounds on the frame determinant

Put $Q_2=\sum_{i<j}(u_i\cdot u_j)^2$. The six pair dots lie in $[-c,u]$, and their sum lies between $-2$ and $-2+8\|\bar u\|^2$, where $\bar u=\sum_i u_i/4$. Evaluate the convex function $Q_2$ at the vertices of this box cut by two sum inequalities. The eigenvalues of $N^\mathsf{T}N$ sum to four. Their deviations $x_i$ from $4/3$ have sum zero and sum of squares $s=2Q_2-4/3$. For a zero-sum triple one has $|x_1x_2x_3|\le s^{3/2}/(3\sqrt6)$. For $s>0$, the extrema on the circle defined by the sum and squared sum have two equal coordinates, as follows from Lagrange multipliers. The bound is immediate at $s=0$. Expanding the eigenvalue product yields

$$
\det(N^\mathsf{T}N)\ge(4/3)^3-(4/3)s/2-(s/3)\sqrt{s/6}.
$$

A complementary estimate uses projection onto the plane perpendicular to an eigenvector for the smallest eigenvalue $e$. The four projected normals still have a positive relation. If $e\le\eta<1$, then after normalizing these projections their pair dots lie in $[-c',u']$, where

$$
c'=\frac{c+\eta/2}{\sqrt{1-\eta}},\qquad
u'=\frac{\max(0,u)+\eta/2}{\sqrt{1-\eta}}.
$$

Indeed, the squared discarded components sum to at most $\eta$, and each pair's discarded product has magnitude at most $\eta/2$. For planar unit vectors with positive relation and pair dots in $[-c',u']\subset(-1,1)$, the four cyclic gaps lie in $[\alpha,\beta]$, with $\alpha=\arccos u'$ and $\beta=\arccos(-c')$. For each of the two opposite pairs, one of the two complementary two-gap sums is at most $\beta$. The chosen two-gap sums overlap in one gap of size at least $\alpha$. The fourth gap is at most $\beta$, so $3\beta\ge2\pi+\alpha$. This is impossible for $c'<1/2$. Otherwise it implies $3c'-4c'^3\le u'$ by the cosine triple-angle identity. The polynomial test

$$
(3(c+\eta/2)-\max(0,u)-\eta/2)(1-\eta)>4(c+\eta/2)^3
$$

therefore proves $e>\eta$. The strict inequality also ensures $c',u'<1$, because $u'\ge0$ and $3c'-4c'^3\le1$ for $c'\ge0$.

If the least eigenvalue is $e$, the two others have sum $4-e$ and squared sum $4+2Q_2-e^2$. Hence the determinant is at least $e(6-Q_2-4e+e^2)$. It is also at least $e^2(4-2e)$, by minimizing the product of the other two eigenvalues subject to their sum and their lower bound $e$. The first bound is concave for $0\le e\le4/3$. The second is increasing in that interval. The checker applies the two bounds on either side of a few rational split points.

## Capped deficits and an extended pair bound

Use the sign-cell deficits $d_J=M-v_J$, cap deficits $e_i$, and the pair-error formula $E(s,q)$ from the [affine pair proof](affine-gram-bm-proof.md). Besides $\sum_Jd_J=14\delta+t$, the deficit coordinates satisfy $0\le d_J\le M$. For each value of $t$, the capped simplex vertices have as many coordinates equal to $M$ as fit, at most one residual coordinate, and zeros elsewhere. In computing graph disagreement bounds, maximize convex functions at these vertices. Along each interval of $t$ for which the number of full coordinates is fixed, the vertex expressions are convex in $t$. Evaluate the endpoints, including the changes at $14\delta+t=kM$.

One also has $h_i\le g(e_i)$, where $g(e)=e$ for $e\ge0$ and $g(e)=4e/7$ for $e<0$. The positive side follows from the thin-slab estimate. On the negative side, condition on the coordinate with largest normal component, whose magnitude is at least $1/\sqrt3>4/7$, to bound the marginal density. At fixed weights, maximize the convex function $\sum_i\lambda_i g(e_i)$ at a capped-simplex vertex. Its maximum over $t\le T$ is attained at $t=T$, because the cap deficits are nondecreasing in every $d_J$. For each candidate cap vector maximize over the bounded weight simplex by assigning available weight to the largest entries first. This supplies an upper bound on $B$. A separate bound used in each height bin is $B\le(1-a)H-4aG/7$, where $G=7\delta$. Indeed, maximize the same convex function over the larger uncapped deficit simplex. A vertex puts all deficits at one proper sign set, so at least one coordinate has $e_i=-G$, while the others have $e_i\le H$. Its weight is at least $a$, giving the stated bound.

The pair error extends to $s\ge3q$ by setting $E(s,q)=s-q$. The Brunn–Minkowski tangent then forces the two positive shifted halfspaces to have empty overlap, and the resulting disagreement is $1-s$. The extension matches the cubic branch in value and first derivative. It remains a convex perspective and is nondecreasing in $s$ and nonincreasing in $q$. Moreover, $E(s_1,q)+E(s_2,q)\le E(s_1+s_2,q)$ by convexity and $E(0,q)=0$. Thus all offset sign cases satisfy $D\ge D_0-E(2H,1-D_0)$. The right-hand side is nondecreasing in $D_0$: on the cubic branch its derivative is $(1-z/3)^2(1+2z/3)$ with $z=2H/(1-D_0)$, and on the linear branch it is constant. Invert it by rational bisection. The checker applies the graph, Gram, and tetrahedron bounds repeatedly, beginning with $t\le\min(M,28\delta)$.

## Reproduction

Run the angular and pair checks listed in the affine pair proof. Verify the additional triples in `WIDE_SHARP_PAIRS` using `wedge_certificate --sharp-pair A B U`. Then run `python3 affine_frame_height_certificate.py` and `python3 affine_frame_height_certificate.py --candidate` from the proof directory.

## 用中央四面体的棱长收紧候选阈值

取有理阈值 $q=79759/10^6$，仍以边长一的立方体计量。若所有块都不超过 $q$，进一步有

$$
t\le 8452/10^6,\qquad |h_i|\le66765/10^6,\qquad
-573370/10^6\le u_i\cdot u_j\le143814/10^6\quad(i\ne j).
$$

这里的法向都朝向中央四面体外侧。任意两张面的交棱长度记为 $L_{ij}$，两面的面积为 $A_i,A_j$。以交棱为底，两面三角形在各自面内的高为 $2A_i/L_{ij},2A_j/L_{ij}$。用三个棱向量计算四面体体积便得
$$
6t=4A_i A_j\sin\theta_{ij}/L_{ij}.
$$
由 $A_i=3t\lambda_i/B$（面积向量平衡及带符号的棱锥体积公式），得到
$$
t=\frac{L_{ij} B^2}{6\lambda_i\lambda_j\sin\theta_{ij}}.
$$
共点时直接取 $t=0$。非共点时四个顶点均在立方体内，所以 $L_{ij}\le\sqrt3$。选择权重最大的两个面，把这条二次上界和前文由四个子式得到的三次上界一起用于同一组权重。即使最小权重趋近零，二次上界仍然有效。

核验时把四个权重依次排成 $a\le b\le c\le d$，总和为一。前三个各用步长 $1/128$ 的闭区间覆盖，用总和和次序界定第四个。对每组权重，按前述截断单纯形顶点求出四个有序的高度上界 $g_1\le g_2\le g_3\le g_4$。由重排不等式，
$$
B\le g_4-(g_4-g_3)(a+b+c)-(g_3-g_2)(a+b)-(g_2-g_1)a.
$$
此外 $B\le (7/8+t_{\max}+7\delta)a$ 来自最小权重所对的三面交点。将两个 $B$ 上界以及正弦、框架行列式的下界代入体积公式，向上取整后回代。程序仅删除被另一个有序上界逐坐标支配的高度向量。以 $1/10^6$ 为舍入单位，二十轮得到上述数值，全部使用整数和分数。[有理数核验程序](affine_frame_height_certificate.py)

在每轮迭代中，两把向外定向的刀的同号区域由八个符号块组成，其中一块为空，一块为中央四面体，余下六块各至多为 $q$。因此同号体积的上界可以取 $\min(7q,6q+t_{\max})$。将此值代入前文的居中两刀比较，得到上述加强界。任何单刀两侧的体积均属于 $[0.433235,0.566765]$。这是因为外侧含七个非空的普通块，内侧含七个普通块和中央四面体。可在本目录运行以下命令复核。

~~~sh
python3 -B affine_frame_height_certificate.py --candidate --edges
~~~
