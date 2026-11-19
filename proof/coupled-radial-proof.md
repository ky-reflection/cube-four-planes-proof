# Radial rearrangement using the three linked depths of each plane

Take four planes through the center of a unit cube. Let $S$ be any set of directions whose membership is constant on each sign cell of the planes. Suppose its spherical area fraction is at least $\mu$, and define
$$
r_0^2=(1-\mu/3)^{-2}-1,
\qquad \kappa=\frac{r_0^2}{1+r_0^2}.
$$
For $8/25\le\kappa$ and $r_0\le4/5$, the [rational checker](coupled_radial_certificate.py) gives a lower bound on the cube volume of $S$ above the disk rearrangement value $\pi r_0^2/4$. It uses the fact that the same normal determines one line depth on each of three pairs of faces.

## Asymmetric annuli

Use the face coordinates and weight $w(r)=(1+r^2)^{-3/2}$ of the [four-line proof](four-line-radial-proof.md). In every face square compare the selected region with the disk $r\le r_0$. As before, the planar area excess is at least the sum of the mismatch integrals with weight $|1-w(r)/w(r_0)|$. For an angular ray, parameterize the inner and outer endpoints by
$$
r_-^2=r_0^2(1-u),\qquad r_+^2=r_0^2(1+v(u)),\qquad
v(u)=\frac{u}{1-2u/5}.
$$
Use a horizon $U\le3/5$ such that $r_+(U)\le1$. The whole annulus lies in the square. With $v$ temporarily an independent variable, the mismatch densities per $du$ and $dv$, respectively, are
$$
\frac{r_0^2}{2}[(1-\kappa u)^{-3/2}-1],\qquad
\frac{r_0^2}{2}[1-(1+\kappa v)^{-3/2}].
$$
For $v_\kappa(u)=u/(1-\kappa u)$, the outer density times $v_\kappa'(u)$ equals the inner density times $(1-\kappa u)^{-1/2}$. Since $\kappa\le 16/41<2/5$, our $v$ and its derivative are at least $v_\kappa$ and its derivative, respectively. Thus the outer density per $du$ is at least the inner one. With $\kappa_*=8/25$, the first five nonzero terms of the inner density give the common lower bound
$$
\frac{3r_0^4}{4(1+r_0^2)}\,p(u),\qquad
p(u)=u+\frac54\kappa_*u^2+\frac{35}{24}\kappa_*^2u^3
+\frac{105}{64}\kappa_*^3u^4+\frac{231}{128}\kappa_*^4u^5.
$$

At any angle where no line intersects the annulus up to a given $u$, one side of the disk boundary up to that $u$ has the wrong membership throughout. The nested annuli imply that the angular integral of the mismatch cost is at least the integral of the common lower density times the angular measure avoiding all line crossings. For one line at actual distance $d$ from the face center the angular measure of crossings is
$$
U_d(u)=2\left[\arccos\min(1,d/r_+)-\arccos\min(1,d/r_-)\right].
$$
The same expression is zero for a line through the face center, ignoring the two rays along it. Use a union bound for the four lines. Write $I(d)=\int_0^U p(u)U_d(u)\,du$ and $W=\int_0^U p(u)\,du$.

## Coupling opposite face pairs

For a unit normal $n$, the line depth on the two faces with axis $j$ is $d_j=|n_j|/\sqrt{1-n_j^2}$. When the denominator is zero the plane misses both faces and contributes zero. The depths are linked by $\sum_j n_j^2=1$. If $I_*$ bounds $\sum_j I(d_j)$ for any such normal, the total over six faces and four cutting planes is at most $8 I_*$. After dividing by 24 to convert areas to cube volume, the excess bound is
$$
\frac{3r_0^4}{16(1+r_0^2)}\left(2\pi W-\frac43 I_*\right).
$$
The argument only uses the normals and remains valid when a face line misses the square.

## Exact evaluation

On each rectangle of the $u$ and $d^2/r_0^2$ grids, the checker uses the upper endpoint of $u$, the lower endpoint of the squared depth for the outer angle, and the upper endpoint of the squared depth for the inner angle. These choices enlarge $U_d$. The cosine comparisons use outward-rounded Taylor bounds on the square of the positive part of the cosine. Integrate the polynomial $p$ exactly. The checker next partitions the squared components $n_j^2$ into bins. It bounds the resulting depth range for each component bin and checks all sorted triples of component bins that can have sum one. It takes the best valid excess for horizons $9/20,1/2,11/20,3/5$, or zero if all the estimates are negative. At $r_0=3/4$ this certifies excess at least $0.006858273546\ldots$. At $r_0=77/100$ it certifies excess at least $0.006123057541\ldots$.
