# Direct certification at the candidate value

Let $m=(1+2\sqrt2)/48$. Work in the [centered alternating vertex chart](vertex-chart-overlap-proof.md). In each row of the off-diagonal score matrix $C$, choose a maximal entry, rescale it to one, and let the other two entries vary independently in $(0,1]$. A four-digit word specifies the maximal column chosen in each row, numbered from zero. Simultaneous permutation of the four parity vertices reduces the 81 choices to the six words $1000,1001,1002,1032,1200,1230$. Tied maxima can use any of the applicable choices.

The tests below prove the bound $M\ge m$ over the entire eight-dimensional charts $1000,1032,1230$, with strict inequality in the last two. Chart $1000$ attains $m$ at the known candidate. These results use the candidate threshold directly. The [follow-up proof](centered-bernstein-completion.md) covers the other three centered charts and gives a partial affine extension.

## A small rational Bernstein bound for a row

All volume expressions here are multiplied by twelve. A row of $C$ has entries $1,x,y$ in pivot order. The three triangle terms share the positive denominator
$$
K=(1+x)(1+y)(x+y).
$$
Their numerators in the same order are $x+y, x^2(1+y), y^2(1+x)$. For a subset of the seven cell volumes, first add their expressions. A triangle occurs with coefficient $+1$ in one singleton and $-1$ in one complementary pair. Thus the combined coefficient on each triangle is $-1,0$, or $1$. This combines the triangle constraints before bounding them.

For any such row combination $N/K$, use a tensor Bernstein representation of degree $(2,2)$ on its two coordinate intervals. If $X_i$, $X_i^{(2)}$ are the degree-two control coefficients for $x$, $x^2$, respectively, then
$$
X=(x_-,(x_-+x_+)/2,x_+),\qquad
X^{(2)}=(x_-^2,x_-x_+,x_+^2).
$$
Use $Y,Y^{(2)}$ for $y$. At index $(i,j)$, the denominator control coefficient is
$$
K_{ij}=X_i+Y_j+X_i^{(2)}(1+Y_j)+Y_j^{(2)}(1+X_i)+2X_iY_j.
$$
All these coefficients are nonnegative. If one is zero, the numerator coefficient there is also zero, since it is the corresponding signed combination of the first three numerator terms. On the strictly positive chart, $N/K$ is a convex combination of the ratios of its nonzero-denominator controls. Their minimum bounds the row sum on the whole box. The [C++ checker](centered_chart_cover.cpp) uses a dyadic coordinate grid and integer outward rounding for each ratio.

The six overlap terms are bounded separately. On one branch the overlap is $\delta^2/d$ for $\delta>0$, with $\delta$ a difference of two edge fractions. The denominator $d$ is a sum of two terms of the form $(a+e)/(a+b)$. This ratio is increasing in $e$, decreasing in $b$, and monotone in $a$ with the other two fixed, so endpoint tests suffice. In addition $d\ge1+\delta$, since its terms are the edge fractions divided by transverse fractions at most one. This allows a finite upper bound even at box boundaries where a denominator might approach zero. The opposite branch reverses the sign of $\delta$ and exchanges the relevant two rows. Combining these overlap bounds with the row bounds gives a lower bound for the sum of any subset of cells. If that lower bound exceeds the number of cells times $12m$, one of the cells exceeds $m$ throughout the box.

## The three covered charts

In $1032$, the single-box witness is the sum of the two pair volumes for splits $02\mid13$ and $03\mid12$. This witness also gives a simple analytic improvement within that chart. Each row contributes its two nonpivot triangle terms to the sum being subtracted, and these two terms sum to at most $1/2$. Indeed, with the notation above, the difference between $K$ and twice the sum of those two numerators is
$$
x(1-x)(1+y)+y(1-y)(1+x)\ge0.
$$
The overlaps in the pair volumes are nonnegative. The two scaled pair volumes therefore sum to at least two, giving $M\ge1/12$. Equality in the triangle inequality for strictly positive $x,y$ forces $x=y=1$ in every row. Uniform rows attain $1/12$, so this also gives the exact restricted minimum.

In $1230$, exact subdivision verifies all boxes directly. In $1000$, the only remaining boxes lie within $20/32768$ in sup norm of one of the two candidate parameter points. Their small singleton labels are 2 or 3, because either is absent from the pivot word. At either point, all entries in its column are $\sqrt2-1$, while the other entries equal one. Interchanging labels 2 and 3 preserves the chart.

The [local checker](centered_bernstein_certificate.py) treats the box around the candidate with small label 2. The six active volumes are the three other singletons and the three pair volumes, all equal to $12m$ at the candidate. Overlaps have zero first derivative there. For a normalized displacement of sup norm one, fix a coordinate attaining its norm and its sign. Each coordinate whose reference value is one is constrained to be nonpositive. In each of the eleven possibilities for the fixed coordinate and sign, the checker verifies a nonnegative combination of the six volume gradients whose value over the entire displacement face is greater than $19/1000$. The minimum of this linear form is evaluated at endpoints using arithmetic in $\mathbb{Q}(\sqrt2)$.

For a displacement of radius $r\le1/1000$, any two entries in the same row sum to at least $7/5$ throughout the segment. Along a unit sup-norm direction a single edge fraction $f=A/(A+B)$ has $|f'|\le 1/(A+B)$, by differentiating the ratio. Also $f''=-2(A'+B')f'/(A+B)$. The second derivative of a triangle term, a product of two fractions, is thus bounded in magnitude by $250/49$. The Taylor error for the triangle is at most $(125/49)r^2$. An overlap is at most $(100/49)r^2$, because its edge-fraction difference starts at zero, has derivative bounded by $10/7$, and its denominator is at least one on its actual branch. For either type of active volume, these errors total at most $14r^2$ in the lower bound. Thus some active volume is greater than $12m+(19/1000)r-14r^2>12m$ for every nonzero displacement in these neighborhoods. This closes the $1000$ cover without assuming any bound on the other charts.

## A flexible linear-programming witness

On a strictly positive parameter box, the [Python verifier](centered_bernstein_certificate.py) writes each rational triangle or overlap term as an affine function $L$ plus a remainder. For a rational term $N/D$, choose $L$ by differentiating at the midpoint and use the control coefficients of $N-DL$ and $D$ at the same degrees to bound the remainder. If both overlap branches are possible, the same affine approximation is checked against both. Each atom depends on at most four variables. The resulting certified affine model retains the dependence on the eight shared coordinates.

A linear program searches for nonnegative weights on the seven cell volumes. Its objective is the lower bound on the weighted sum minus $12m$ times the sum of the weights. The weights are combined before bounding the affine terms and signed remainders. The verifier rounds the proposed weights and rechecks the entire inequality with rational arithmetic. An LP failure leaves the box unresolved, and an LP result by itself never excludes a box.

For the following concrete box in chart $1231$, every individual lower enclosure falls below the threshold. The displayed combination of four cell constraints nonetheless excludes the whole box, with a certified margin greater than $0.00015365$ in scaled volume. The best individual lower enclosure has margin about $-0.00923$ in the same units. This example exercises the use of nonuniform weights as well as the independent exact verifier.

## Reproduction

Run these commands from the directory containing this note.

```sh
c++ -O3 -std=c++20 centered_chart_cover.cpp -o /tmp/centered-chart-cover
/tmp/centered-chart-cover 1032 0
/tmp/centered-chart-cover 1230 35
/tmp/centered-chart-cover 1000 96
python3 -B centered_bernstein_certificate.py --local-1000
python3 -B centered_bernstein_certificate.py 1231 \
  241/250:1 19/20:493/500 407/1000:443/1000 119/125:247/250 \
  203/500:221/500 122/125:1 191/500:209/500 939/1000:39/40 \
  --weights 0 28800837 0 0 32682813 19213600 19302749
```

The three chart runs visit respectively 1, 929,569 and 12,763 boxes. The $1000$ run reports `incomplete` because it hands off 20 boxes to the local lemma, with `open-radius=20/32768` and zero pending boxes. The other two runs finish with zero unresolved or pending boxes. Removing the final `--weights` argument lets the LP rediscover a certificate for the last box.
