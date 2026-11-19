# Four two-vertex cap sides cannot form a perfect matching

Suppose four planes cut the unit cube and the largest block volume is at most $2/25$. Choose one closed side of each plane. It is impossible for each chosen side to contain exactly two strictly signed cube vertices and for these pairs of vertices to form a perfect matching. Here the six remaining vertices for each side may lie on its boundary or strictly outside it. In particular, if all vertex signs are strict, the eight vertices cannot all have one exterior sign, and they cannot all have three exterior signs. For those two applications, orient the cuts as the exterior sides of the tetrahedron of triple intersections, or by the positive dependence when they are concurrent. [Three-cut reduction](affine-gram-bm-proof.md)

## Geometry of a two-vertex cap

Every chosen side has volume at least $9/25$, since the other side is split by the remaining three cuts into at most eight blocks. It contains at least two strict vertices. A side with no strict vertex has volume zero by convexity. If its only strict vertex is $(0,0,0)$, write the side as $a x+b y+c z\le d$. The three neighboring vertices give $a,b,c\ge d>0$, so the cap lies in the corner simplex $x+y+z\le1$, of volume $1/6$. Consequently, if every vertex has just one exterior sign, the four exterior sides have at least two vertices each and exactly eight incidences in total. They must select two apiece and form a perfect matching. If every vertex has three exterior signs, the same argument applies to the four interior sides.

The two strict vertices in a chosen side must be adjacent. If they differ in two coordinates, the two other vertices of their face have the same midpoint as the given pair, which would be both strictly in the side and in the complementary closed side. If they differ in three coordinates, another antipodal pair gives the same contradiction.

Work in coordinates $x,y,z\in[0,1]$. For a chosen side with strict vertices $(0,0,0)$ and $(1,0,0)$, write its halfspace as $a y+b z\le \ell(x)$, with $\ell$ affine and $a,b>0$. Indeed, each endpoint is strict while its neighbors at $y=1$ or $z=1$ are outside or on the plane. This gives $0<\ell(x)\le\min(a,b)$ for all $x\in[0,1]$. The cap lies in the reference prism $y+z\le1$.

Two chosen sides at opposite parallel edges would have disjoint interiors, as their reference prisms are complementary. Their symmetric difference would have volume at least $18/25$, whereas it is a union of at most eight blocks and so has volume at most $16/25$.

For a perfect matching the number of edges in each coordinate direction is even. To see this, one side of the corresponding midplane contains four vertices, and its internal matched edges use an even number of them. Thus the matching uses either four parallel edges or two edges in each of two directions. The former contains opposite parallel edges. In the latter, the two parallel edges must belong to one face, since they cannot be opposite. They use all four vertices of that face, and the remaining edges match the opposite face in the other direction.

After relabeling coordinates, the first two chosen sides have edges at $y=0,z=0$ and $y=1,z=0$, and the other two have edges at $x=0,z=1$ and $x=1,z=1$. The last two caps lie in the reference prisms $z\ge x$ and $z\ge 1-x$. Therefore in the wedge $0\le z<\min(x,1-x)$ they are absent. The part of the union of the first two caps in this wedge lies in only three blocks, corresponding to the two singleton patterns and the pair pattern for these caps.

## A volume bound in the wedge

At fixed $(x,z)$, the first cap covers an interval starting at $y=0$ and the second an interval ending at $y=1$. We bound the length for the first cap; reflection in $y$ gives the same bound for the second.

Divide $a y+b z\le\ell(x)$ by $a$ and write the upper endpoint as

$$
L=u+t(2x-1)-wz,\qquad w=b/a.
$$

Let $v=|t|$ and $a_0=\min(1,w)$. The vertex inequalities give $v\le u$ and $u+v\le a_0$. The full cap has volume $(u^2+v^2/3)/(2w)$ by integrating its triangular sections at each $x$. Put $k=18/25$. The cap volume implies $u^2+v^2/3\ge kw$, and hence $k\le w\le 1/k$, since $u^2+v^2/3\le (u+v)^2\le a_0^2$.

We have

$$
(a_0-v)^2+v^2/3\ge kw,\qquad
0\le v\le a_0/2.
$$

The left side decreases on this interval, so

$$
v\le v_*=\frac{3a_0-\sqrt{12kw-3a_0^2}}4.
$$

In fact $v_*\le a_0/2$, since $kw/a_0^2\ge k>1/3$. The function $\sqrt{kw-v^2/3}-v|2x-1|$ decreases in $v$ on $[0,v_*]$. Using $u\ge\sqrt{kw-v^2/3}$ at the given $v$, and the defining equality at $v_*$, gives

$$
L\ge a_0-(1+|2x-1|)v_*-wz.
$$

The right side is concave in $w$ separately on $[k,1]$ and $[1,1/k]$. For $w\le1$ it is affine in $w$ plus a positive multiple of $\sqrt{12kw-3w^2}$, the square root of a positive concave quadratic. For $w\ge1$ it is affine in $w$ plus a positive multiple of $\sqrt{12kw-3}$. Its minimum is therefore at $w=k,1$, or $1/k$. The respective endpoint values are $k(1-z)$, $1-(1+|2x-1|)(3-\sqrt{12k-3})/4-z$, and $1-z/k$.

Since $(3-\sqrt{12k-3})/4<4/25$ and $1+|2x-1|=2(1-\min(x,1-x))$, the middle value is at least $(17/25)(1-z)$ when $z\le\min(x,1-x)$. The first value also exceeds this. Each interval thus has length at least

$$
L_0(z)=\min\left(\frac{17}{25}(1-z),1-\frac{25}{18}z\right).
$$

Here $L_0(z)>0$ for $z\in[0,1/2]$, so taking the nonnegative part of the cap endpoint preserves the bound. The union of the two intervals covers at least $\min(1,2L_0(z))$ of each horizontal segment. Integrating over $x,z$ gives a volume of at least

$$
\int_0^{1/2}(1-2z)\min(1,2L_0(z))\,dz
=\frac{1550224607}{6352328664}>\frac6{25}.
$$

The integral is obtained by splitting at $9/34$, where $2L_0$ first falls below one, and $144/319$, where the two branches of $L_0$ meet. The volume exceeds that of any three blocks under the assumed cap, a contradiction.

## Three cap edges incident to a common vertex

There is also a constraint on three chosen two-vertex sides with edges incident to the same cube vertex. The three sides cannot all be exterior in the tetrahedron orientation. Regardless of orientation, they cannot occur if one side of the fourth cut contains the common vertex and its three edge neighbors.

To see this, normalize the $y$-leg of one of the caps as $\ell(x)/a$. Let its minimum and maximum be $p$ and $P$, and let $w=b/a$. Its volume is

$$
\frac{p^2+pP+P^2}{6w}\le \frac{p^2+p+1}{6}.
$$

For $w\ge1$ use $P\le1$. For $w\le1$ use $P\le w$ and the fact that $p^2/w+w$ increases for $w\ge p$. Put $s=24/35$, which satisfies $(s^2+s+1)/6<9/25$. The cap volume lower bound forces $p>s$. Interchanging $y,z$ gives the same bound on the other leg. Thus the cap contains the entire reference prism $y+z\le s$.

Move the common vertex to the origin. The intersection of the three caps contains the region

$$
K=\{(x,y,z)\in[0,1]^3:x+y\le s,\ x+z\le s,\ y+z\le s\}.
$$

Its volume is $s^3/4$ by integrating over the six coordinate orders. If the three chosen sides are exterior, their intersection gives just one block because the all-exterior sign pattern is absent. That block would have volume at least $s^3/4=3456/42875>2/25$.

If one side of the fourth cut contains the common vertex and all three neighbors, that side contains the corner simplex $x+y+z\le1$. At a fixed sum $S>1$, the region $K$ requires each coordinate to be at least $S-s$. Subtracting $S-s$ from each coordinate gives a nonnegative triangle with coordinate sum $3s-2S$, whose area in the volume integral is $(3s-2S)^2/2$. Hence

$$
\operatorname{vol}(K\cap\{x+y+z\le1\})
=\frac{s^3}4-\frac{(3s-2)^3}{12}
=\frac{10366}{128625}>\frac2{25}.
$$

This subset lies in a single block, again a contradiction. In particular, suppose the common vertex is a strict triple, missing label $j$, and the other seven vertices are strict singletons, four of label $j$ and one of each other label. The three non-$j$ cuts would select exterior edges incident to the common vertex, so this pattern is excluded. The pattern obtained by reversing every vertex sign is also excluded, because the fourth cut then has constant sign on the common vertex and its three neighbors.

## One or two cap edges with the other cuts fixed on the corner simplex

More generally, suppose one or two chosen sides of a cell are two-vertex caps whose edges meet at a common vertex. Assume each remaining chosen side contains that vertex and its three edge neighbors. Each remaining side then contains the corner simplex $x+y+z\le1$. If there is just one two-vertex cap, its reference prism $y+z\le s$ intersects this simplex in volume

$$
\int_0^s r(1-r)\,dr=\frac{s^2}{2}-\frac{s^3}{3}>\frac2{25}.
$$

Two caps cannot select the same edge: their union lies in the same reference prism of volume $1/2$, and their individual volumes are at least $9/25$. Their symmetric difference would be at most $1-18/25$, yet it is at least $1-16/25$ because its complement contains at most eight blocks. If the two incident edges are distinct, their reference prisms can be written as $x+z\le s$ and $y+z\le s$. Their intersection has volume $s^3/3$. The part outside the corner simplex has triangular sections of area $(2s-1-z)^2/2$ for $0\le z\le2s-1$, and is empty at larger $z$. The single cell therefore contains volume at least

$$
\frac{s^3}{3}-\frac{(2s-1)^3}{6}>\frac2{25}.
$$

This rules out, for example, a strict two-positive vertex if its two positive sides are two-vertex caps and both negative sides contain all three neighbors of the vertex.
