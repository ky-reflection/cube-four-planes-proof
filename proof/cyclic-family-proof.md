# Three cyclic central cuts and an arbitrary fourth cut

Let the first three planes through the center of a unit cube have normals obtained by cyclically permuting one vector. The minimum possible largest block, with an arbitrary fourth plane, is

$$
m=\frac{1+2\sqrt2}{48}.
$$

Equality occurs when the first three normals, up to cube symmetries, relabeling and sign and scale of plane equations, are the permutations of $(\sqrt2,-1,-1)$, and the fourth plane is $x+y+z=0$ in centered coordinates.

## Reduction of the fourth plane

In a first-three sign cone and its antipode, both of volume $W$, let $S(h)$ be the first cone's volume on the positive side of the fourth plane at offset $h$. For $h\ge0$, the two negative-side volumes are $W-S(h)\ge W-S(0)$ and $S(-h)\ge S(0)$. Reflection covers $h<0$. Recentring the fourth plane therefore cannot increase the largest block.

A centered arrangement with largest block below $1/12$ must have fourteen nonempty sign cells. The number is even by antipodality, and four concurrent planes produce at most fourteen cells. The first three normals must be independent, since three dependent central normals give at most six cells before the last cut.

There are four antipodal pairs of first-three cones. Exactly one pair remains unsplit. This pair has volume at most $M$ per cone if the largest final block is $M$, while each split cone has volume at most $2M$. Two pairs of volume at most $M$ would give total volume at most $6M<1/2$ for one representative from each pair. Cyclic symmetry permutes the three mixed-sign pairs transitively. Hence the unsplit pair must be the allsame pair.

Write the first three linear forms as $f_i$. Up to sign, a fourth linear form leaving the allsame pair unsplit has the form $\sum_i\lambda_i f_i$ with $\lambda_i>0$. In the cone $f_j>0$, $f_k<0$ for $k\ne j$, the cut is determined by

$$
\lambda_j f_j>\sum_{k\ne j}\lambda_k|f_k|.
$$

For the equal-weight cut the two sides have the same volumes across all three such cones, by cyclic symmetry. For arbitrary weights, choosing $j$ of largest weight enlarges the displayed side compared with the equal-weight cut. Choosing $j$ of smallest weight enlarges the complementary side. The equal-weight cut therefore minimizes the largest block among directions leaving the allsame pair unsplit. Its plane is $x+y+z=0$ whenever the first three normals are independent.

## Computing the body-diagonal cut

Use normals cyclically permuting $(a,-\alpha,-\beta)$, with

$$
a>0,\qquad \alpha=1-d,\quad \beta=1+d,\qquad 0\le d\le1.
$$

This normalization follows from sign changes of all three equations together, a coordinate permutation, and scaling, when the generating normal has both positive and negative coordinates. If its coordinates all have the same weak sign, an entire coordinate octant lies in an allsame cone. Its volume is at least $1/8$, ruling out the required allsame unsplit pair.

Let $u=x+y+z$. Orient the diagonal by $u$ for $a<2$, and by $-u$ for $a>2$, so all four oriented normals have a positive dependence. At $a=2$ the first three normals are dependent and the largest block is at least $1/12$.

Let $p$ be the sign disagreement probability between the oriented diagonal and any one of the three other planes. Let $q$ be the sign disagreement probability between any two of the other planes, and $T$ the volume of one allsame first-three cone. The triple is either allsame or has exactly two pair disagreements, so $q=2(1-2T)/3$. The three types of four-cut blocks have volumes

$$
T,\qquad A=\frac12-\frac{2p+q}{4},\qquad B=\frac{p+q-1}{2}.
$$

There are two blocks of the first type and six each of the other two types. For a one-positive or one-negative pattern, the three equal-sign entries disagree pairwise exactly twice when they are mixed, giving the formula for $A$. For a two-positive pattern, sum the four disagreements across its two-by-two split of the indices. This count is four on the pattern and its opposite, and two on every other feasible sign pattern, giving the formula for $B$. In particular $A+B=q/4$.

### Face polygons for $a<2$

Scale coordinates to $[-1,1]^3$ while retaining relative volume. A polygon on a cube face of area $E$ contributes $E/24$ to the normalized volume of its cone from the center. On the face $x=1$ write the other coordinates as $Y,Z$. Then

$$
f_1=a-\alpha Y-\beta Z,\quad
f_2=-\beta+aY-\alpha Z,\quad
f_3=-\alpha-\beta Y+aZ.
$$

Put $c=2-a$, $r=1-Y$, $s=1-Z$, and $L=\beta^2+\alpha a$. The allnegative polygon is bounded by $\alpha r+\beta s\le c$, $\beta r\le c+as$, and the nonnegative coordinate axes. The remaining inequality $\alpha s\le c+ar$ follows from the first one. Its vertices are

$$
(0,0),\quad(c/\beta,0),\quad
\bigl(c(a+\beta)/L,c(\beta-\alpha)/L\bigr),\quad(0,c/\beta).
$$

These points lie in $[0,2]^2$. For $a\le\beta$, the third vertex's first coordinate is at most $c$. For $a>\beta\ge1$, $c<1$, $L\ge2$ and $a+\beta\le4$ give the same square containment. On $x=-1$, the sum $f_2+f_3=2+(a-\beta)Y+(a-\alpha)Z$ is nonnegative throughout the square, since $|a-\beta|+|a-\alpha|\le2$. Thus only the three positive faces contribute to this cone. Taking the polygon area divided by eight gives

$$
T(a,d)=\frac{(2-a)^2(a+1+3d)}
 {16\beta(\beta^2+\alpha a)}
=\frac{(2-a)^2}{16(1+R)},\qquad
R=\frac{d^2(3+d-a)}{a+1+3d}.
$$

Here $R$ increases in $d$, as is seen by logarithmic differentiation for $d>0$. Also $-R_a/(1+R)\le4d^2(1+d)/(1+3d)^2\le1/2$. The logarithmic derivative of $T$ in $a$ is therefore strictly negative. Thus $T$ decreases in $a$ and $d$, and $q(a,d)\ge q(a,0)=(4+4a-a^2)/12$.

For $p$, let $N_1,P_2,P_3$ be the face areas of $f_1<0,f_2>0,f_3>0$, respectively, and $E$ the area of $f_3>0,u<0$. The triangle $u<0$ on this face has area $1/2$, lies in $f_1>0$, and has zero overlap with $f_2>0$. Adding both equal-sign types on this face over all three forms and then using cyclic and antipodal symmetry gives

$$
p=\frac58+\frac{N_1-P_2-P_3+2E}{12}.
$$

Intersecting each line with the square edges, or with the edges of the triangle $u<0$ for $E$, gives the following areas.

| Quantity | $a\ge2d$ | $a\le2d$ |
|---|---|---|
| $N_1$ | $c^2/(2\alpha\beta)$ | $2(1-a/\beta)$ |
| $P_2$ | $(a-2d)^2/(2a\alpha)$ | $0$ |
| $P_3$ | $(a+2d)^2/(2a\beta)$ | $4d/\beta$ |
| $E$ | $2d^2/[a(a+\beta)]$ | $1/2-(a+\alpha)^2/[2\beta(a+\beta)]$ |

Consequently, with $p_0(a)=(19-6a+a^2)/24$,

$$
\begin{aligned}
p(a,d)&=p_0(a)+\frac{d^2}{12a}
\left[-\frac{(2-a)^3}{2(1-d^2)}+\frac4{a+\beta}\right] &&(a\ge2d),\\
p(a,d)&=\frac{13}{24}
+\frac{-3a^2-4+8\beta-\beta^2}{12\beta(a+\beta)} &&(a\le2d).
\end{aligned}
$$

The formulas agree at the common boundary. At $d=1,a<2$ use the second formula, which has no singularity.

### Lower bound for $a<2$

If $\sqrt2\le a<2$, then $\max(A,B)\ge q/8\ge q(a,0)/8\ge m$, since $q(a,0)$ increases up to $a=2$.

For $2/3\le a\le\sqrt2$, we have $p(a,d)\ge p_0(a)$. In the first formula the bracket is at least $[4-2(2-a)^2]/(a+2)>0$, using $d\le a/2$. In the second formula the derivative in $\beta$ has the sign of $(3a^2+4)(a+2\beta)-(8+a)\beta^2$, which changes at most once from positive to negative. The minimum is thus at $d=a/2$ or $d=1$. The first endpoint was already treated. At the second, the difference from $p_0(a)$ is $[3a-a^2-4/(a+2)]/24>0$, since the numerator is positive at $a=2/3$ and increases up to $a=\sqrt2$. Hence

$$
B\ge B(a,0)=\frac{3+2a-a^2}{48}\ge m,
$$

where the last inequality follows because $2-\sqrt2<2/3\le a\le\sqrt2$.

For the remaining range, it suffices to use $m<2/25$.

- For $a\le12/25$, monotonicity gives $T\ge T(12/25,1)=2527/31250>2/25$.
- For $12/25\le a\le2/3$ and $d\le3/5$, it gives $T\ge T(2/3,3/5)=325/3816>2/25$.
- For $12/25\le a\le2/3$ and $d\ge3/5$, the second formula gives $p\ge\min(p(a,a/2),p(a,1))$. Both endpoints decrease in $a$ on this interval. In fact $p(a,1)=(19-3a-4/(a+2))/24$, and $p(a,a/2)=13/24+(12+12a-13a^2)/[12(a+2)(3a+2)]$ has a positive decreasing numerator in its last fraction and an increasing denominator. Therefore $p\ge p_0(2/3)=139/216$. Also $q\ge q(12/25,3/5)=\tfrac23(1-2\cdot14801/137600)>523/1000$, so $B>(139/216+523/1000-1)/2>2/25$.

### The range $a>2$ and equality

For $a>2$, the allpositive triple cone on $x=1$ is a quadrilateral with vertices

$$
(1,1),\quad(2/a,1),\quad
\left(\frac{a\beta+\alpha^2}{a^2-\alpha\beta},
      \frac{a\alpha+\beta^2}{a^2-\alpha\beta}\right),\quad(1,2/a).
$$

Its area is $(a-2)^2(a+1)/[a(a^2-\alpha\beta)]$, and the cone has no patch on the negative faces. Therefore

$$
q=\frac23-\frac{(a-2)^2(a+1)}{6a(a^2-\alpha\beta)}
\ge q(a,0).
$$

The face-area computation for the other pair uses $N_1=0$, $P_2=2(1-\beta/a)$, $P_3=2(1-\alpha/a)$, and the same $E=2d^2/[a(a+\beta)]$. Reversing the diagonal sign gives

$$
p=\frac{17}{24}-\frac1{3a}-\frac{d^2}{3a(a+\beta)}
\ge p(a,0)-\frac1{3a(a+2)}.
$$

At $d=0$, $B(a,0)\ge5/48$, since
$(a^2-2a+2)/[6a(a-1)]\le1/6$. Thus $B\ge5/48-1/48=1/12>m$ throughout this range.

Equality in the lower bound is possible only at $a=\sqrt2,d=0$, by the strict inequalities above and strict decrease of $T(a,d)$ in $d>0$. There, $A=B=m$ and $T=(3-2\sqrt2)/8<m$, so the diagonal cut attains the bound. If positive fourth weights are unequal, choosing a largest weight makes the inclusion used in the direction comparison strict. Indeed the first-three normals are independent, so their coordinates realize a neighborhood of zero and hence an open set of positive-measure points between the two splitting lines within the cone. Finally, at this solution the two sides of each split cone both have volume $m$, and their dividing plane crosses the cone interior. Any nonzero fourth offset strictly increases one of the antipodal paired volumes. This proves the equality statement for arbitrary fourth planes.
