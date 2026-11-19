# A cell cannot contain opposite cube vertices

Suppose four planes cut the unit cube and the largest cell volume is at most $2/25$. Then no closed sign cell contains opposite cube vertices. Consequently every closed sign cell contains at most three cube vertices. By the three-cut reduction, the volume assumption implies that every triple intersection lies in the open cube and has independent normals. [Three-cut reduction](affine-gram-bm-proof.md)

Center the unit cube at the origin. Orient the cuts as the exterior sides $u_i\cdot x>h_i$ of the tetrahedron of triple intersections, with $\|u_i\|=1$. At concurrence choose an orientation with a positive dependence. In both cases the normals satisfy a positive relation. The all-negative cell lies in the open cube, and the all-positive cell has no points on the cube boundary. In the concurrent case both weak sign regions reduce to the concurrence point.

## Pair angles under an antipodal sign constraint

Let $M_*=2/25$, $G=7M_*-1/2$, and $H=8M_*-1/2$. The positive cap of cut $i$ contains at most seven cells, and the negative side at most eight. Its signed cap deficit $e_i=1/2-\operatorname{vol}\{u_i\cdot x>h_i\}$ therefore satisfies $-G\le e_i\le H$. The centered slab bound gives $|h_i|\le |e_i|$. [Affine Gram proof](affine-gram-bound-proof.md)

Let the two vertices be $\pm(\sqrt3/2)w$ with $w$ unit. For a cell containing both vertices, call its positive signs type $N$ and its negative signs type $P$. The type $N$ cuts have $h_i\le 0$, $|e_i|\le G$, and the type $P$ cuts have $h_i\ge 0$, $|e_i|\le H$. For either type,

$$
|s_i|=|u_i\cdot w|\le \frac{2|h_i|}{\sqrt3}\le \frac{2|e_i|}{\sqrt3}.
$$

In particular, orthogonal projection into $w^\perp$ gives four nonzero vectors. Denote their unit directions by $\widehat u_i$. These still satisfy a positive dependence. Their dots satisfy

$$
\widehat u_i\cdot\widehat u_j
=\frac{u_i\cdot u_j-s_i s_j}{\sqrt{(1-s_i^2)(1-s_j^2)}}.
$$

For each pair of cuts the affine disagreement mass is at most $8M_*$ and the affine agreement mass is at most $7M_*$. The Brunn–Minkowski pair bound relates the former to the centered disagreement $D_0$. Its loss at $q=1-D_0$ is the sum of the errors

$$
E(s,q)=\frac{s^2}{3q}-\frac{s^3}{27q^2}
$$

for the positive and negative deficit totals separately, where each total $s$ lies in $[0,3q]$. For $s>3q$ the same bound holds with $E(s,q)=s-q$, since the Brunn–Minkowski tangent forces the overlap to be empty. The extended error is increasing in $s$ and bounded by $s$. Reversing one cut gives the analogous bound for agreement. [Brunn–Minkowski pair bound and angular certificates](affine-gram-bm-proof.md)

Write $A$ for the affine mass upper bound and $s_1,s_2$ for upper bounds on the two deficit totals. Here $A<1-s_1-s_2$. If $q<\max(s_1,s_2)/3$, the extended bound would force the affine mass to be at least $1-s_1-s_2$, a contradiction. Thus the centered mass is at most $1-\max(s_1,s_2)/3$. It is also at most $A+s_1+s_2$ by recentering the two cuts one at a time. Starting with the lesser of these upper bounds $U$, we repeatedly improve it by $A+\sum_k E(s_k,1-U)$ using the polynomial expression. The certified angular anchors bound the centered dots from these mass bounds. The projection formula and rational Taylor bounds for the cosine give the following angle intervals in radians. [Exact rational checker](antipodal_vertex_certificate.py)

| Pair | Lower angle | Upper angle |
| --- | ---: | ---: |
| $NN$ | $1.235$ | $2.220$ |
| $NP$ | $1.150$ | $2.275$ |
| $PP$ | $1.128$ | $2.625$ |

## Exclusion in the projection plane

There is at least one cut of each type because the common cell contains cube boundary points. There are three possibilities for the number of $N$ cuts.

With one $N$ cut, the three $P$ directions must lie in two arcs at angular distances between $1.150$ and $2.275$ from the $N$ direction. Each arc has length $1.125$. Two $P$ directions must lie in one arc, violating their pair separation of at least $1.128$.

For the remaining cases use the four cyclic gaps $g_1,g_2,g_3,g_4$ between the directions. Each is less than $\pi$ because the directions have a positive dependence and are pairwise neither equal nor opposite. Suppose all gaps lie in $[\alpha,a]$, and the two diagonal-pair angles are at most $b,c$, respectively. Then

$$
a+b+c\ge 2\pi+\alpha.
$$

Indeed, $g_1+g_2$ is either at most $b$ or at least $2\pi-b$, and $g_2+g_3$ is either at most $c$ or at least $2\pi-c$. If both are small, $g_4=2\pi-(g_1+g_2)-(g_2+g_3)+g_2\ge 2\pi-b-c+\alpha$. If both are large, the same argument applies to the complementary sums. If one is small and one large, $|g_1-g_3|\ge 2\pi-b-c$, whereas $|g_1-g_3|\le a-\alpha$.

With two $N$ cuts the only pair that might exceed $2.275$ is the pair of $P$ cuts. Thus $a+b+c\le 2(2.275)+2.625<2\pi+1.128$. With three $N$ cuts every pair is at most $2.275$, giving the same contradiction.

## The number of cube vertices in a cell

If the closure contains four affinely independent cube vertices, it contains their convex hull, whose volume is an integer multiple of $1/6$. In the absence of antipodal pairs, four vertices consist of one choice from each of the four opposite pairs. Use as representatives the vertices $v_1,\ldots,v_4$ of one parity class. A selection $\{\sigma_i v_i\}$ with all $\sigma_i$ equal or one differing is affinely independent. A selection with two of each sign is a full cube face, as follows by permuting the four representatives by cube symmetries.

To exclude the face case, rescale the cube to $[-1,1]^3$ and take the face $x=1$. For a cell containing this face, each cut's oriented side containing it is $x\ge a_i y+b_i z+c_i$ with $c_i+|a_i|+|b_i|\le1$. Its volume fraction is at least $1-8M_*$ since the opposite side contains at most eight cells. Conditional on $(y,z)$ uniform in the face square, the cap fraction is bounded above by $(1-a_i y-b_i z-c_i)/2$. Hence $c_i\le16M_*-1$. Put $R=\max(|y|,|z|)$. The maximum of the four graphs is at most $16M_*-1+(2-16M_*)R$. Since $\mathbb E R=2/3$, the cell containing the full face has volume fraction at least $(1-8M_*)/3>M_*$.
