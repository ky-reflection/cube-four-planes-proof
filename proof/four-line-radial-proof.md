# A four-line improvement to radial rearrangement

Let $S\subset\mathbb S^2$ be the directions of the $2+2$ regions of four planes through the center of a unit cube. Suppose $S$ has spherical area fraction at least $\mu$, and put
$$
r_0=\sqrt{(1-\mu/3)^{-2}-1}.
$$
If $3/4\le r_0\le4/5$, then the volume of the cube in these directions is at least
$$
\frac{\pi r_0^2}{4}+
\frac38\frac{r_0^4}{1+r_0^2}\frac{619}{10000}
\;>\;\frac{\pi r_0^2}{4}+\frac{47}{10000}.
$$
The same statement holds for any set of directions whose membership on each of the six face charts is constant on every cell of an arrangement of at most four straight lines.

## Face coordinates and the mismatch cost

Use the face charts of $[-1,1]^3$, writing $(x,y)\in[-1,1]^2$ for coordinates in a face and $r=\sqrt{x^2+y^2}$. In these coordinates the spherical area element is $w(r)\,dx\,dy$, where $w(r)=(1+r^2)^{-3/2}$. The volume of the corresponding radial cone in the **unit** cube is the face area divided by 24. Summing over six faces gives the total volume.

The disk $B=\{r\le r_0\}$ in each face gives the six axial spherical caps whose total area fraction is $\mu=3(1-1/\sqrt{1+r_0^2})$. With $w_0=w(r_0)$, the planar area of the portions of $S$, minus the combined disk area, is at least the sum over the faces of
$$
\int_{S_{\rm face}\triangle B}|1-w(r)/w_0|\,dx\,dy.
\tag{1}
$$
Indeed, subtract the disk area and add and subtract the spherical area difference divided by $w_0$. The spherical area difference is nonnegative. The remaining integrand has the claimed sign on both sides of the disk boundary.

Along a ray at polar angle $\theta$, write $r=r_0(1\pm y)$, with $0\le y\le1/4$. Both radial coordinates stay in the square because $r_0\le4/5$. Put $\kappa=r_0^2/(1+r_0^2)\le2/5$. The weighted cost per $dy$ of a mismatch, on either side of the disk boundary, is at least
$$
\frac32\frac{r_0^4}{1+r_0^2}\,p(y),
\qquad p(y)=y(1-y)(2-y)(1+3y/4).
\tag{2}
$$
For the inner side, put $z_-=\kappa(2y-y^2)$ and use
$$
(1-z_-)^{-3/2}-1\ge\tfrac32z_-+\tfrac{15}{8}z_-^2.
$$
The second term gives the factor $1+5z_-/4\ge1+3y/4$, since $\kappa\ge9/25$ and $y\le1/4$. Include the polar factor $r\,dr=r_0^2(1-y)\,dy$. For the outer side, let $z=\kappa(2y+y^2)$. Use
$$
1-(1+z)^{-3/2}\ge\tfrac32 z(1+z)^{-5/2}.
$$
The polar factor is now $r_0^2(1+y)\,dy$. The same bound (2) follows because
$$
\log\frac{(1+y)(2+y)}{(1-y)(2-y)}-\tfrac52\log(1+z)
\ge 3y-(2y+y^2)\ge3y/4.
$$
Exponentiation supplies at least the factor $1+3y/4$.

## How much of a circle can four lines cover?

If none of the four lines crosses the radial interval $[r_0(1-y),r_0(1+y)]$ at angle $\theta$, membership in $S_{\rm face}$ is constant throughout that interval. On such a ray, up to the first crossing or $y=1/4$, one entire side of the disk boundary is therefore a mismatch. Integrate (2) over that side and then over $\theta$.

For a line at distance $d r_0$ from the face center, the angular measure of rays crossing it in this annulus is
$$
U(d,y)=2\left[
\arccos\!\left(\min(1,d/(1+y))\right)
-\arccos\!\left(\min(1,d/(1-y))\right)
\right].
$$
The rays along a line through the center form a measure-zero exception. The formula follows by writing the intersection radius as $d r_0/\cos\theta$ after rotating the line to have perpendicular direction at angle zero.

By a union bound over the four lines and integration in $y$, the contribution to (1) on one face is at least
$$
\frac32\frac{r_0^4}{1+r_0^2}\left[
2\pi\int_0^{1/4}p(y)\,dy
-4\sup_{d\ge0}\int_0^{1/4}p(y)U(d,y)\,dy
\right].
\tag{3}
$$
The bound remains valid if the bracket is negative. The checker confirms it is positive.

## Exact one-line quadrature

The code [four_line_radial_certificate.py](four_line_radial_certificate.py) partitions both $y\in[0,1/4]$ and $d\in[0,5/4]$ in units of $1/1024$. Beyond $d=5/4$, the line misses every annulus. On each rectangle, $U(d,y)$ is bounded above by taking the upper endpoint of $y$, the lower endpoint of $d$ in the first angle, and the upper endpoint of $d$ in the second angle. The first angle is rounded up to the grid $1/2000$, and the second down, by comparison with alternating Taylor bounds for the cosine. The exact integrated weight on each $y$-interval is the difference of $y^2-y^3/2-5y^4/16+3y^5/20$ at its endpoints.

The resulting one-line upper bound is $0.068728526021$, so the bracket in (3) is at least $0.061938142011>619/10000$. Apply (3) on all six faces and divide by 24 for the radius-dependent excess stated at the beginning. At $r_0=3/4$ the directly verified excess is $0.004703427659$. Using the simpler bracket bound $619/10000$ still gives more than $47/10000$.
