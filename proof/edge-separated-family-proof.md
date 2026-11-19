# Four cuts whose two-positive cells avoid the cube edges

Work in $Q=[-1,1]^3$ and normalize volumes by $8$. Let $v_1,\dots,v_4$ be the four cube vertices whose coordinate product is positive. They satisfy $\sum_i v_i=0$ and $v_i\cdot v_j=4\mathbf 1_{i=j}-1$.

Assume the planes are labeled and oriented so that the only positive sign at $v_i$ is $i$, and the signs at $-v_i$ are the complements. Assume also that no relatively open interval of a cube edge belongs to a two-positive cell. The centered cuts in this chart are exactly the cuts with normals

$$
n_i=v_i-b,\qquad b=\sum_i w_i v_i,\qquad w_i>0,\qquad \sum_i w_i=1,
$$

up to positive rescaling of plane equations. Within this family the smallest possible largest cell is

$$
m=\frac{1+2\sqrt2}{48},
$$

attained only when the weights are proportional to $(\sqrt2-1,1,1,1)$, up to permutation. The same minimum holds in the affine extension of this chart, and equality forces all offsets to vanish.

## Characterizing the plane equations

Let the cut equations be $n_i\cdot x-h_i=0$, and set $C_{ij}=-n_i\cdot v_j>0$ for $i\ne j$. The assumed signs at the two vertices also give $C_{ij}>|h_i|$. Every edge of the cube has endpoints $v_i,-v_j$ for some $i\ne j$. The two planes with indices $k,l\notin\{i,j\}$ both cross this edge, and the other two planes have constant signs there. If the crossings occurred at distinct points, the interval between them would have two positive signs. Thus the crossings coincide. The same argument applies to the opposite edge $v_j,-v_i$. Their two crossing equations are

$$
\begin{aligned}
(C_{ki}+h_k)(C_{lj}-h_l)&=(C_{li}+h_l)(C_{kj}-h_k),\\
(C_{kj}+h_k)(C_{li}-h_l)&=(C_{lj}+h_l)(C_{ki}-h_k).
\end{aligned}
$$

Subtracting and adding gives $C_{ki}C_{lj}=C_{li}C_{kj}$ and $h_k(C_{li}+C_{lj})=h_l(C_{ki}+C_{kj})$. The six disjoint $2\times2$ minors imply a factorization $C_{ij}=r_i c_j$ for $i\ne j$. To see this, take logarithms and subtract row and column constants so that all off-diagonal entries of row $1$ and column $1$ vanish, as does the entry $(2,3)$. The disjoint-minor equations for row pairs $(1,2),(1,3),(1,4)$ give respectively $X_{24}=X_{23}=0$, $X_{32}=X_{34}$, and $X_{42}=X_{43}$. The equations for row pairs $(2,3),(2,4),(3,4)$ then give $X_{34}=X_{43}=X_{42}=X_{32}=0$. Hence the original logarithms are additive in the row and column.

Since $\sum_j n_i\cdot v_j=0$, the diagonal entry is $r_i\sum_{j\ne i}c_j$. Put $w_i=c_i/\sum_jc_j$. Rescaling each equation makes its normal $v_i-b$, since $(v_i-b)\cdot v_j=4(\mathbf 1_{i=j}-w_j)$. The offset relations above give $h_i/r_i$ constant, so the rescaled offsets all equal one number $h$. The assumed vertex signs hold whenever $|h|<4\min_i w_i$. The equations and volume formulas extend continuously to the endpoints.

## Centered volumes

On the face $x=1$, relabel so that $v_i=(1,1,1),v_j=(1,-1,-1)$. Let the weights of $i,j$ be $a,b$ and the other weights be $c,d$. Put $u=y+z$, $v=y-z$. The face becomes the diamond $|u|+|v|\le2$, and the scores $v_r\cdot(1,y,z)$ of $i,j$ are $1\pm u$. The other two scores are $-1\pm v$. Thus $i,j$ are the two highest scores everywhere on the face. Subtracting $b\cdot(1,y,z)$, with $b$ here the vector from the normal formula, leaves the order unchanged.

Consider the quadrants with $u\ge0$, so that $i$ is the highest score. Write $U=u$, $V=|v|$. The line where the second score $n_j\cdot(1,y,z)$ equals zero meets the $U$-axis at $2(c+d)/(1+a-b)$. It meets the hypotenuse $U+V=2$ at the point with $V=2a/(a+q)$, where $q$ is either $c$ or $d$, according to the quadrant. The singleton $i$ occupies the triangle between this line and the vertex $(U,V)=(2,0)$, since the highest score is positive and the other two scores are below the second. Its $U,V$ area is $4a^2/[(1+a-b)(a+q)]$.

The conversion from $U,V$ area to normalized cone volume is $1/48$. Summing the two quadrants gives the singleton contribution from this face

$$
\frac{a^2}{12(1+a-b)}
\left(\frac1{a+c}+\frac1{a+d}\right)
=\frac{a^2}{12(a+c)(a+d)}.
$$

The three faces where $i$ is one of the two highest scores give the full singleton volume

$$
S_i=\frac{w_i^2(1+2w_i)}{12\prod_{j\ne i}(w_i+w_j)}.
\tag{1}
$$

The two-positive cell with indices $i,j$ is entirely in the pyramid from the origin to their common cube face. Indeed, comparison of their scores with the two lower scores gives $x\ge |y|,|z|$ in the coordinates above. On that face there are four other cells, namely the singletons $i,j$ and the antipodes of the singletons $k,l$. Subtracting their triangle volumes from $1/6$ gives

$$
P_{ij}=\frac1{12}
\left[1-\frac{(ab-cd)^2}{(a+c)(a+d)(b+c)(b+d)}\right].
\tag{2}
$$

For the simplification, write the four weights as a $2\times2$ table with opposite corners $a,b$ and the other two corners $c,d$. The four triangle terms, multiplied by $12$, are its squared entries divided by the corresponding products of row and column sums. Each entry differs from that product by $\pm(ab-cd)$. Expanding the squares gives a sum of $1+(ab-cd)^2/[(a+c)(a+d)(b+c)(b+d)]$.

The face-pyramid argument also proves the converse in the characterization. The two positive scores strictly exceed the two negative scores, giving $x>|y|,|z|$. Hence the two-positive cell has no open interval on a cube edge. A common affine offset preserves this order argument.

## The centered minimum

Sort the weights as $a\le b\le c\le d$, and put $x=a/d,\ y=b/d,\ z=c/d$. By (1) the largest-weight singleton satisfies

$$
12S_d=\frac{3+x+y+z}{(1+x)(1+y)(1+z)}
\ge \frac{5+x}{4(1+x)}.
$$

The left side decreases strictly in $y,z$, so equality requires $y=z=1$. If $x\le\sqrt2-1$, the right side is at least $12m$, with equality only at $x=\sqrt2-1$.

If $x\ge\sqrt2-1$, use the split $\{a,d\}\mid\{b,c\}$. The square root of the correction in (2) satisfies

$$
\frac{|x-yz|}{\sqrt{(x+y)(x+z)(1+y)(1+z)}}
\le \frac{1-x}{2(1+x)}
\le \frac{\sqrt2-1}{2}.
$$

For $x\ge yz$, bound the numerator by $x-x^2$ and the denominator below by $2x(1+x)$. For $yz>x$, the ratio increases in $y,z$ because, for example, $z/(yz-x)>\tfrac12[(x+y)^{-1}+(1+y)^{-1}]$. Evaluating at $y=z=1$ proves the first inequality. Since $12m=1-(\sqrt2-1)^2/4$, equation (2) gives $P_{ad}\ge m$. Equality requires $x=\sqrt2-1$, then the singleton inequality requires $y=z=1$. At those weights the three large singletons and all two-positive cells have volume $m$. The small singleton has volume $(3-2\sqrt2)/8$. The normals in this case are the body diagonal and the three cyclic permutations of $(\sqrt2,-1,-1)$, up to scaling and cube symmetries.

## Exact affine volumes in this chart

Reflect the cube and reverse all signs if needed to arrange $0\le h\le4\min_i w_i$. The four vertices of the interior tetrahedron are $-h v_i/(4w_i)$. Its normalized volume is

$$
t=\frac{h^3}{768\prod_i w_i}.
\tag{3}
$$

For the determinant computation, use the four rows $(v_i,1)$, whose determinant has absolute value $16$. The determinant with last column $w_i$ has absolute value $4$, since its four cofactors in that column agree and $\sum_i w_i=1$. Rescaling each row and dividing the spatial determinant by $6\cdot8$ yields (3).

The triple-positive cell complementary to $i$ is the cone cut out by the three positive inequalities with indices other than $i$, with apex $-h v_i/(4w_i)$. The corresponding central cone exits only the three faces in the directions of $-v_i$. The displaced apex is closer to all three by the same factor $1-h/(4w_i)$, and farther from the other faces. Its volume is therefore $S_i(1-h/(4w_i))^3$.

For singleton $i$, split each assigned face pyramid into the four simplex cones defined by the signs of $u,v$. In a quadrant where $i$ has the highest score, the scores of $i,j$ on the four simplex vertices are, respectively,

$$
(0,\ 2Q,\ 4q,\ 4(1-a)),\qquad (0,\ 2Q,\ 4q,\ -4a),\qquad Q=c+d.
$$

The singleton is the portion where the first score exceeds $h$ and the second is below $h$. For a linear score with values $g_r$ at a simplex's four vertices, the lower-tail fraction at $h$ is $\sum_{r:g_r<h}(h-g_r)^3/\prod_{s\ne r}(g_s-g_r)$. This follows by slicing a simplex, or by the divided-difference formula for the integral of a truncated power. For these two lists and $h\le4\min_iw_i$, subtracting the first lower-tail fraction from the second gives

$$
\frac{(h+4a)^3}{(4a)(4a+2Q)(4a+4q)}
-\frac{h^3}{(2Q)(4q)}\left(\frac1{4a}+\frac1{4(1-a)}\right).
$$

Each simplex cone has normalized volume $1/24$. Summing the first term gives $S_i(1+h/(4w_i))^3$. The second term gives $t$, using $(1/c+1/d)/(c+d)=1/(cd)$. Thus the affine singleton volume is $S_i(1+h/(4w_i))^3-t$.

For the two-positive cell with weights $a,b$ and complementary weights $c,d$, put $\tau=h/4$ and $K=(a+c)(a+d)(b+c)(b+d)$. The cell lies entirely in the pyramid over the common face of its two highest scores, as in the centered case. Its volume is

$$
P_{ab}(h)=
\frac{d^2(1-\tau/d)^3/[(a+d)(b+d)]-
      c^2(1-\tau/c)^3/[(a+c)(b+c)]}{12(d-c)}.
\tag{4}
$$

For a derivation, let $z_r$ be the tetrahedron vertices' coordinates perpendicular to the face, with the face at height $1$. The two high coordinates are $z_c=\tau/c,\ z_d=\tau/d$, and the two low coordinates are $z_a=-\tau/a,\ z_b=-\tau/b$. In barycentric coordinates let the two negative coordinates of the cell be $-A,-B$, where $A,B\ge0$. The two other coordinates range over the segment with sum $1+A+B$. The endpoints' heights are $z_r+A(z_r-z_a)+B(z_r-z_b)$ for $r=c,d$. Integrating the segment length below height $1$, first over the segment and then over $A,B\ge0$, gives

$$
-t\sum_{r=c,d}\frac{(1-z_r)^3}{\prod_{s\ne r}(z_r-z_s)}.
$$

Substitution of (3) gives (4). At $c=d$ take the divided-difference limit, and at $h=0$ use continuity. In particular (4) is equivalently

$$
P_{ab}(h)=
\frac{(a+b)cd+ab(c+d)+3\tau(cd-ab)-3\tau^2
      +\tau^3(c+d-cd+ab)/(cd)}{12K}.
$$

The opposite cell swaps $(a,b)$ with $(c,d)$.

## The affine minimum

Keep the weights sorted $a\le b\le c\le d$ and set

$$
x=a/d,\qquad y=b/d,\qquad z=c/d,\qquad
r=\frac{h}{4d}\in[0,x].
$$

The affine singleton associated with $d$, multiplied by $12$, is

$$
F(r)=L(1+r)^3-Cr^3,\qquad
L=\frac{3+x+y+z}{(1+x)(1+y)(1+z)},\qquad
C=\frac{1+x+y+z}{xyz}.
\tag{5}
$$

Its derivative has positive constant and linear coefficients. If the quadratic coefficient is nonnegative, $F$ is increasing for $r\ge0$. Otherwise the derivative has just one positive root and changes from positive to negative there. In either case $F$ takes its minimum on any closed subinterval of $[0,x]$ at an endpoint.

At $r=x$,

$$
F(x)=\frac{(3+x+y+z)(1+x)^2}{(1+y)(1+z)}
-\frac{x^2(1+x+y+z)}{yz}.
$$

Holding either $y$ or $z$ fixed, this has the form $\mathrm{const}+A/(1+s)-B/s$ with $A,B>0$. Its derivative can change sign only from positive to negative, since $(s/(1+s))^2$ is increasing. Thus the minimum on $[x,1]^2$ occurs at a corner. The three corner values up to symmetry are $2,\ 2+x-x^2,\ (5+11x-5x^2-3x^3)/4$, all at least $5/4$. If $h>0$ and all cells have volume at most $m$, it follows that $L<12m$. The centered singleton bound then gives $x>\sqrt2-1>2/5$.

We can also assume $r<7/80$. For fixed $r=7/80$, the same form $\mathrm{const}+A/(1+s)-B/s$ holds in each of $x,y,z$ in (5), so its minimum on $[2/5,1]^3$ occurs at a corner. The values, ordered by the number of coordinates equal to $2/5$, are

$$
\frac{1970021}{2048000},\quad
\frac{17697947}{14336000},\quad
\frac{15686423}{10035200},\quad
\frac{3123839}{1605632}.
$$

All exceed $959/1000>12m$. The endpoint bound at $x$ and the derivative of $F$ therefore give $F(r)>12m$ for $7/80\le r\le x$.

For $r<7/80$, the function in (5) decreases strictly in $y,z$. In the derivative with respect to $y$, the ratio of the positive term to the magnitude of the negative term is

$$
\frac{r^3}{(1+r)^3}
\frac{1+x+z}{2+x+z}\frac{(1+x)(1+z)}{xz}
\left(\frac{1+y}{y}\right)^2
\le \left(\frac{7}{80}\right)^3\left(\frac72\right)^4<1.
$$

The same comparison holds for $z$. Evaluating at $y=z=1$ gives

$$
F(r)\ge L_{\mathrm{eq}}(1+r)^3-r^3(3+x)/x
\ge L_{\mathrm{eq}}(1+3r),\qquad
L_{\mathrm{eq}}=\frac{5+x}{4(1+x)}.
$$

For the second inequality use $L_{\mathrm{eq}}\ge3/4$ and $r(3+x)/x\le (7/80)(17/2)<9/4$. Put $u=(1-x)/(1+x)$, $u_0=\sqrt2-1$. Then $0\le u<u_0<1/2$, $L_{\mathrm{eq}}=(3+2u)/4$ and $12m=(3+2u_0)/4$. Since $F(r)\le12m$, we obtain

$$
r\le\frac{2(u_0-u)}{3(3+2u)}.
\tag{6}
$$

Now use the two cells in the split $\{a,d\}\mid\{b,c\}$. In their cubic formulas the linear terms in $\tau$ are opposites and both cubic coefficients are positive. Thus their maximum is at least the centered value minus $\tau^2/(4K)$. By the centered pair bound that value is at least $(1-u^2/4)/12$. Put $S=1+x+y+z$ and $D=(x+y)(x+z)(1+y)(1+z)$, so $\tau^2/K=r^2S^2/D$. The expression $S^2/D$ decreases in $y,z$ since $2/S<1/(x+y)+1/(1+y)$, with the corresponding inequality for $z$. At $y=z=x$ it equals $(2-u)^2(1+u)^2/[4(1-u)^2]$.

Combining these estimates with (6), the ratio of the pair loss to the centered gap $(u_0^2-u^2)/48$ is at most $4Q(u)/3$, where

$$
Q(u)=\frac{u_0-u}{u_0+u}
\left[\frac{(2-u)(1+u)}{(3+2u)(1-u)}\right]^2.
$$

The logarithmic derivative of $Q$ is

$$
-\frac1{u_0-u}-\frac1{u_0+u}-\frac2{2-u}
+\frac2{1+u}-\frac4{3+2u}+\frac2{1-u}<0.
$$

The first two terms sum to less than $-4$. The positive terms are at most $2$ and $4$, and the remaining negative terms are at most $-1$ each. Thus $Q(u)\le Q(0)=4/9<3/4$. One of the two affine pair cells is strictly larger than $m$, a contradiction. This excludes every positive $h$, and reflection excludes negative $h$. The centered equality conditions finish the proof.

## The vertex sign pattern follows from the edge condition

Assume every triple of cuts is independent and its intersection lies in the open cube. Assume no cut contains a cube vertex. Orient the cuts as the exterior sides of their interior tetrahedron, or with a positive dependence of the normals if they are concurrent. Then no point of the cube boundary has all four evaluations nonpositive or all four nonnegative. The nonpositive region is the closed tetrahedron, or the point of concurrence. For the other claim, a positive linear combination of the four oriented affine evaluations is a constant that is strictly negative for a tetrahedron and zero at concurrence. In the zero case all four nonnegative evaluations would have to vanish.

Suppose no relatively open interval of a cube edge has two strictly positive and two strictly negative evaluations. Every vertex is then of one of two types. It has either a unique positive label $i$, or a unique negative label $j$. Call these the singleton and triple types respectively.

Adjacent vertices of the same type must have the same label. If two distinct singleton labels occurred, just those two cuts would change sign along the edge. Distinct crossings give either an all-negative or a two-positive interval. Coincident crossings give a boundary point with all evaluations nonpositive. The triple case follows by reversing all signs in this argument.

For an edge from singleton $i$ to triple $j$, the labels must be different. If $i=j$, all four cuts switch signs. At the crossing of cut $i$, let $u$ other cuts have crossed already and let $v$ cross simultaneously with it. The positive counts immediately before and after are $1+u$ and $u+v$, both in $\{1,3\}$. Thus $(u,v)$ is one of $(0,1),(0,3),(2,1)$. The crossing point has all evaluations nonpositive in the first two cases, and all nonnegative in the last case, each a contradiction.

When $i\ne j$, the other two cuts must cross at the same point to avoid a two-positive interval. There is at most one edge for each **ordered** pair $(i,j)$. Otherwise the two coincidence points are distinct boundary points of the line common to the other two cuts. Their segment is the entire line section of the cube because the line meets the open cube at its triple intersection with cut $i$. But cut $i$'s evaluation is strictly positive at both coincidence points and hence on their whole segment, a contradiction.

Consider any connected component of the subgraph induced by the singleton vertices, or by the triple vertices. Its vertices share one label. Its edge boundary has at most three edges by the preceding ordered-pair bound. There are both types of vertices: if all vertices had one type then connectivity would force a single label, and some cut could not meet the open cube. A nonempty proper vertex set of the cube graph with boundary at most three has one or seven vertices. To see this, take the smaller of the set and its complement, of size $s\le4$. At $s=2,3,4$, there are respectively at most $1,2,4$ induced edges, by simplicity and bipartiteness, so the boundary count $3s-2e$ exceeds three. A component of size seven is also impossible: its complementary vertex has opposite type and a different label. In the singleton case the cut with the component's label would be positive at all eight vertices. In the triple case it would be negative at all eight. Both contradict its open-cube triple intersections.

Hence all components are single vertices, so the two types are the two parity classes. Each singleton vertex has three incident edges to triple vertices with distinct labels, all different from its own label. It uses all three allowed ordered pairs starting at that label, so no two singleton vertices share a label. The same holds for the triple vertices. The only triple vertex without an edge to singleton $i$ is its antipode, so that vertex has label $i$. This is precisely the vertex sign pattern at the start of the proof, up to a cube reflection and relabeling.

Thus under the triple-intersection and strict-vertex assumptions, the edge condition alone implies largest volume at least $m$, with the same equality characterization.

## Excluding vertices on the cutting planes at the candidate threshold

Suppose the same edge condition holds and the largest volume is $M\le m<2/25$. Every vertex has at least one strictly positive and one strictly negative evaluation by the first observation in the preceding section. If exactly one plane passes through a vertex, the other three evaluations there have either one or two positive signs. In the first case the zero plane's directional derivatives along the three inward edges must all be nonpositive, since a strictly positive derivative would create a two-positive interval. In the second case they must all be nonnegative. In either case the plane supports the cube and cannot have open-cube triple intersections.

It remains to consider a vertex on two planes. The other two cuts have opposite signs. Consequently the two zero planes cannot have opposite strict directional signs along any inward edge. Write $X,Y,Z\in[0,1]$ for the inward coordinates from the vertex. The two normal vectors in these coordinates have coordinatewise products nonnegative, and each has both a positive and a negative component because its plane meets the open cube. Their union of positive coordinate indices and union of negative indices are disjoint nonempty sets. One of the unions is a singleton, and its coordinate is nonzero in both vectors. Reverse both planes if necessary and call that coordinate $Z$, so the two chosen halfspaces have forms
$$
Z<aX+bY,\qquad Z<cX+dY,\qquad a,b,c,d\ge0.
$$
Reversing both does not change their sign disagreement volume $D$.

Each chosen cap has volume in $[1-8M,8M]\subset [9/25,16/25]$, since either halfspace is split into at most eight blocks. Let $s=a+b$, $t=c+d$, and $\phi(u)=\min(1,u)$ on $u\ge0$. The two cap volumes are $\mathbb E\phi(aX+bY)$ and $\mathbb E\phi(cX+dY)$, hence are bounded above by $s/2$ and $t/2$ respectively. In particular $s,t\ge18/25$.

Their common cap volume is $\mathbb E\min(1,aX+bY,cX+dY)$. The integrand is concave in the four coefficients. For given $s,t$, interpolate each pair of coefficients between its two endpoints, for example $(s,0),(0,s)$ for $(a,b)$. The common volume is at least the minimum of the four endpoint values. When the two cap slopes use different coordinates the overlap equals
$$
I(s,t)=\int_0^1 (1-u/s)_+(1-u/t)_+\,du.
$$
When they use the same coordinate, the integrand is the smaller of the two survival probabilities and is at least their product. Thus in all cases
$$
D\le \min(16/25,s/2)+\min(16/25,t/2)-2I(s,t).
$$
Since $I$ is increasing in both arguments, for this upper bound one may replace each argument above $K=32/25$ by $K$. By symmetry assume $18/25\le s\le t\le K$. For $s\le1$ the bound becomes $(t-s)/2+s^2/(3t)$, increasing in $t$ and convex in $s$. Its maximum is at $t=K$ and $s=18/25$ or $1$. For $s,t\ge1$ the bound is
$$
\frac{s+t}{2}-2+\frac1s+\frac1t-\frac{2}{3st},
$$
which is convex in each variable separately and so is maximized at a corner of $[1,K]^2$. Evaluating all the endpoints gives $D\le16727/38400<11/25$. On the other hand, with the original exterior orientation at least one of the eight agreement patterns is absent, namely the all-positive pattern. Consequently $D\ge1-7M>11/25$, a contradiction. This argument applies even when one of the planes contains an entire cube edge.

Combining this with the strict-vertex classification proves that **any** four-cut configuration with all triple intersections inside the cube, no two-positive edge interval, and largest volume at most $m$ must be one of the equality configurations above. By the three-cut reduction in the [global lower-bound proof](affine-gram-bm-proof.md), the triple-intersection hypothesis follows from the volume bound $M\le m$.

## A uniform positive edge interval for any strict counterexample

Let $L$ be the greatest length, in a unit cube, of a relatively open edge interval with two strictly positive and two strictly negative signs, and put $L=0$ if there is none. There is a single constant $\epsilon>0$ such that any configuration with largest volume strictly below $m$ has $L\ge\epsilon$.

Indeed, use the compact parameter space of normalized plane normals and bounded offsets described in the global lower-bound proof. Restrict to $M\le m$ and to exterior orientations. This is closed: the orientation condition may be written as the existence of $\lambda_i\ge0$ with $\sum_i\lambda_i=1$, $\sum_i\lambda_i n_i=0$ and $\sum_i\lambda_i h_i\ge0$. All triple intersections are interior and all triples are independent on $M\le m$, so every such dependence has strictly positive weights. At concurrence either simultaneous sign convention has the same two-positive edge lengths.

Remove the open neighborhoods of all equality configurations supplied by the [local proof](local-candidate-proof.md). The remaining parameter set is compact and has $L>0$ everywhere by the result just proved. The function $L$ is lower semicontinuous, since any closed subinterval of a strict two-positive interval keeps all its signs under small perturbations. It therefore attains a positive minimum on the remaining set if that set is nonempty. Every strict counterexample belongs to this set. A [quantitative version](edge-stability-proof.md) gives an explicit cutoff.
