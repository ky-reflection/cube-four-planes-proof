# Arbitrary independent offsets for column-weight directions

Let $m=(1+2\sqrt2)/48$, and use the side-two cube with normalized volumes. Consider normals whose off-diagonal scores satisfy $C_{ij}=c_j$ after positive row rescaling, with $c_j>0$. Equivalently their directions are $v_i-b$, where $b=\sum c_i v_i/\sum c_i$ lies inside the parity tetrahedron. The offsets $H_i$ can vary independently over all values that preserve the alternating vertex signs. In this entire seven-dimensional family the least possible largest volume is $m$. The only equality configurations are the centered candidate and its cube symmetries. In particular, all nonconcurrent cuts in this family have a strictly larger block.

This family extends the [common-offset result](edge-separated-family-proof.md). Independent offsets usually separate the crossings on a cube edge. The proof uses the exact volume recurrence from the [affine face formulas](affine-vertex-volumes-proof.md) and verifies boxes directly at $m$.

## A useful fact for the entire alternating affine chart

For any positive off-diagonal scores $C$, including nonfactorable ones, every triple of planes intersects inside the cube whenever the offsets preserve the strict alternating vertex signs. Each triple is independent by the rank and positive relation of the evaluation matrix. At zero offsets its intersection is the center. Join zero offsets to the specified offsets by a segment. The vertex signs stay strict along it. At a cube face, two of the three cuts have the same type, both adding a positive triangle or both removing one. Their corner triangles lie at opposite corners, and their hypotenuses cannot meet on the face. In fact each triangle lies strictly on its own side of the diagonal between the two other corners, since its legs have length strictly less than a full edge. Thus the triple point cannot meet the cube boundary as the offsets vary, so it remains interior.

After the reflection that makes the positive weighted offset sum nonnegative, a nonzero such sum places each $p_j$ strictly on the negative side of its omitted cut. All eight choices of signs of the three meeting cuts occur locally around $p_j$. Hence all fifteen proper sign sets have positive volume. This makes the face recurrence applicable throughout the strict affine chart, with no additional condition on the offsets.

## Stable formulas for box evaluation

Use cube symmetries and common scaling to sort $c=(x,y,z,1)$, $0<x\le y\le z\le1$, and parametrize $(x,y,z)=(suv,su,s)$ with $0<s,u,v\le1$. Write $d_0=y$, $d_i=x$ for $i>0$, and $H_i=d_i t_i$, $-1<t_i<1$. This is the full set of offsets with the prescribed signs. All ratios of a smaller column to a larger one are products of a subset of $s,u,v$, so common vanishing factors can be cancelled before interval evaluation.

Put $S=\sum c_i$, $q=\sum H_i v_i/S$, and $\bar H=\sum c_i H_i/S$. After reflection take $\bar H\ge0$, or equivalently $y(t_0+t_1)+zt_2+t_3\ge0$. The point $p_j=q-(\bar H/c_j)v_j$ is the intersection of the three planes other than $j$. Indeed $n_i\cdot q=H_i-\bar H$, and $n_i\cdot v_j=-c_j$ when $i\ne j$. On each face, cancel the contribution of $H_j$ before evaluating its height. For a face with default pair $B$, set $s_i=+1$ for $i\in B$ and $-1$ otherwise. Then
$$
\nu_F\cdot p_j=\frac1S\sum_{i\ne j}\frac{d_i}{c_j}(s_i c_j-s_j c_i)t_i.
$$
All ratios in this formula have numerator at most denominator. The difference between two heights uses $\bar H( s_i/c_i-s_j/c_j)$ without introducing independent errors for the common translation $q$.

For each pair of corner triangle types on a face, the overlap depends only on the corresponding offset difference. If the adding and removing cuts are $j$ and $i$, and the corner labels are $a,b$, its area fraction is
$$
\frac{(H_i-H_j)_+^2}{2S(c_a+c_b)}.
$$
The sum of the longitudinal legs minus one is $(H_i-H_j)/(c_a+c_b)$ and the sum of their ratios to the transverse legs is $S/(c_a+c_b)$. Factor out $\max(d_i,d_j)$ in the numerator before evaluating the quotient to keep the ratios bounded even as a column vanishes. The triangle legs are $(c_a\pm H_i)/(c_a+c_b)$, which can be evaluated as the centered share times $1\pm(d_i/c_a)t_i$.

For a cap feasibility test, write the three off-diagonal entries of row $i$ as $a\le b\le c$ and use $H_i=a t_i$. The same sum-of-three-uniforms expansion as in the [rescaled case](affine-vertex-volumes-proof.md#all-offsets-for-coordinate-rescaled-tetrahedral-directions) gives the signed deficit
$$
e_i=\frac{a}{b+c}\frac{(1+a/b+a/c)t_i-(a/b)(a/c)t_i^3/3}{(1+a/b)(1+a/c)}.
$$
At a hypothetical maximum of $m$, it must lie between $1/2-7m$ and $8m-1/2$. The checker also uses the volume of every three-cut cell before splitting by the omitted plane. The latter is one pyramid sum $\Phi_{p_j}(A_J+A_{J\cup\{j\}})$, bounded above by $2m$ and by $m$ when the second cell is the absent all-positive cell.

## Joint lower bounds on the seven-variable boxes

The [checker](column_weight_affine_cover.cpp) first evaluates scalar interval bounds on the preceding formulas. In boxes left over, it gives the cell volumes a common affine dependence on the seven normalized box coordinates. It keeps an outward-rounded interval remainder for every computed quantity. Multiplication uses the product of the two deviations from the centers for the nonlinear term. For division, if $q+L$ is the proposed affine part of $a/b$, the exact remainder is $[a-(q+L)b]/b$. The constant and linear coefficients cancel up to known dyadic rounding errors, and the remaining numerator uses the input remainders and the product of $b-b_0$ with $L$. The positive-part square has a continuous derivative of Lipschitz constant two, so its first-order error at an input center is between zero and the squared input displacement. These rules are implemented in [fixed-point affine arithmetic](affine_box_linear.hpp).

The fourteen volume forms are combined with nonnegative weights proposed by a linear program. For fixed weights the minimum of their combined affine part on the normalized box is its center minus the sum of the absolute combined slopes. The check accepts only after the resulting lower inequality, including the weighted remainders, is verified by integer arithmetic. Thus the floating-point weight search does not enter the validity of an exclusion. The checker compares with a dyadic upper bound on $m$. The only exceptions are boxes wholly inside the following exact local neighborhood.

## Exact local exclusion with independent offsets

At $c=(\sqrt2-1,1,1,1)$, $H=0$, twelve cells have volume $m$. Write $r$ for the sup norm of the change in $(c_0,c_1,c_2,H_0,H_1,H_2,H_3)$. The [local certificate](column_weight_local_certificate.py) differentiates the exact pyramid recurrence at the candidate. The overlaps and their first derivatives vanish there. For each of the fourteen faces of the unit displacement box, it checks an exact nonnegative combination of the twelve active derivatives that is greater than $1/80$ everywhere on that face. A nonzero displacement therefore has at least one active first-order gain greater than $r/80$.

Here are uniform bounds for the remainder. Work first at $r\le1/200$. Every pair of columns sums to at least $7/5$, the least column is at least $2/5$, the other columns are at least $99/100$, and $S\ge27/8$. The exact triangle share $f=(A\pm H)/(A+B)$, differentiated along a unit sup-norm direction, satisfies $|f'|\le2/(A+B)$ because $0\le f\le1$, and $|f''|\le4|f'|/(A+B)$. A triangle, half a product of shares, thus has first-order error at most $(150/49)r^2$. An overlap is at most $r^2/2$. A pair cell appears as the default on one face and can leak onto four others. A singleton or triple appears as a single triangle minus overlaps on three faces. Summing these contributions bounds the absolute face-area remainders of any proper nonempty cell by $17r^2$, and the sum of its absolute area changes by $6r$. Its total reference face area is less than $1/2$, since its reference volume is at most $m$ and its heights from the center equal one.

For the individual point heights, the exact bound is $|\nu_F\cdot p_j|\le(13/4)r$, using the three noncancelled offset coefficients in the formula above. The difference from their linear part at the candidate is at most $12r^2$. To see the latter, in $q-(\bar H/c_j)v_j$ the change of the factor $1/S$ in $q$ costs at most $12r^2/(27/8)^2$. The normalized column vector changes in sum of absolute values by at most $6r/(27/8)$, and the change of $1/c_j$ is at most $r/(2/5)^2$. Together these give the asserted error. Similarly a difference of two heights is at most $(18/5)r$ and has first-order error at most $15r^2$, using that the two reciprocal columns include at most one small column. Expand (1) as a base face-area sum with one height multiplier and at most two corrections with height differences. The total first-order error is at most
$$
\left(\frac{17+12/2+(13/4)6}{6}+2\frac{15/2+(18/5)6}{6}\right)r^2<18r^2.
$$
Consequently for $0<r\le1/1600$ one active cell exceeds $m$, since $r/80-18r^2>0$. The candidate at $r=0$ has maximum exactly $m$. The code encloses $\sqrt2-1$ by adjacent dyadic grid points and hands off only boxes wholly in this neighborhood.

## Reproduction

From this directory run

```sh
python3 -B column_weight_local_certificate.py
c++ -O3 -std=c++20 column_weight_affine_cover.cpp -o /tmp/column-cover
/tmp/column-cover 98 5000000 24
```

The run visits 852,217 nodes. It hands off 16 boxes to the local proof and has zero unresolved and zero pending boxes. This covers all seven parameters of the column-weight family. Both the scalar and affine volume enclosures in the joint test were cross-checked against exact volume evaluation on 281 rational configurations in 64 boxes. The general alternating affine chart still has five additional direction parameters beyond these column weights.
