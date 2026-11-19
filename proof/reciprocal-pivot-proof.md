# Reciprocal maximum scores and a two-cut agreement bound

Use the strict alternating vertex chart, with normalized cube volume and each row's maximum score rescaled to one. Put $M=0.079759>(1+2\sqrt2)/48$. Orient by center reflection so the all-positive cell is absent. The cap deficits $e_i$, central tetrahedron volume $t_{\rm cen}$, and $G=7M-1/2$ then satisfy the [deficit identities](affine-pyramid-offset-proof.md#cap-deficits-and-tetrahedron-bounds) if all block volumes are at most $M$. In particular
$$
-G\le e_i\le 8M-1/2,\qquad \sum_i e_i\le 2G+3t_{\rm cen}.
$$

Call two rows reciprocal if each has its score at the other row's label equal to one. This includes ties for the maximum. The proof excludes the following domains, with all offsets independently free to preserve the vertex signs.

| Direction condition | Exclusion |
| --- | --- |
| One reciprocal pair | The maxima of its two rows' nonpivot scores cannot be bounded by $(8/15,8/15)$, $(0.23,1)$, $(0.25,0.95)$, or $(0.30,0.85)$, in either order. The other two cuts need no direction restriction. |
| Two disjoint reciprocal pairs, i.e. pattern $1032$ up to relabeling | The eight nonpivot scores cannot all be at most $347/500=0.694$. They may approach zero independently in this exclusion. |

## The two-cut inequality

Write the cube coordinates as $X,Y,Z\in[-1,1]$, and put $U=(Y+Z)/2$, $V=(Y-Z)/2$. Thus $|U|+|V|\le1$ and the normalized volume density in $X,U,V$ is $1/4$. For a reciprocal pair relabeled $0,1$, divide each score and offset by $1+(x_i+y_i)/2$, where $x_i,y_i$ are its two nonpivot scores. The oriented affine evaluations become
$$
U+s_0 X+\tau_0 V-r_0,\qquad -U+s_1 X+\tau_1 V-r_1,
$$
where
$$
s_i=\frac{x_i+y_i}{2+x_i+y_i},\qquad
|\tau_i|=\frac{|x_i-y_i|}{2+x_i+y_i},\qquad
r_i=\frac{2H_i}{2+x_i+y_i}.
$$
These formulas follow from $n_i=\sum_{j\ne i}C_{ij}(v_i-v_j)/4$. If $x_i,y_i\le b_i$, then
$$
|\tau_i|+|r_i|\le s_i,\qquad
(1+b_i)s_i+|\tau_i|\le b_i.
\tag{1}
$$
The first uses $|H_i|<\min(x_i,y_i)$, and the second follows from
$s_i+|\tau_i|=2\max(x_i,y_i)/(2+x_i+y_i)$.

Let $S=s_0+s_1$, $R=r_0+r_1$ and $T=\tau_0+\tau_1$. At fixed $X,V$, the two cuts have equal signs exactly when $U$ lies between the two thresholds. The interval length before clipping to the diamond is $|R-SX-TV|$. Its clipped length is at most
$\min(|R-SX-TV|,2(1-|V|))$.

Take $X,V$ independently uniform on $[-1,1]$ for the following integral. By (1), $|R-TV|\le S$. Direct integration in $X$ gives
$$
\mathbb E|R-SX-TV|=S/2+R^2/(2S)+T^2/(6S).
$$
The clipping loss is at least the expectation of $(|R-SX|-2(1-|V|))_+$, by averaging the two signs of $V$ and using the convexity of $w\mapsto(|w|-c)_+$. For fixed $X$, $|R-SX|\le2$ by (1), so the integral in $|V|$ is $(R-SX)^2/4$. The agreement probability is therefore at most
$$
S/2-S^2/12+(1/(2S)-1/4)R^2+T^2/(6S).
$$
Set $s=S/2$, $r=|R|/2$, $\tau=|T|/2$. The bound becomes
$$
F(s,\tau,r)=s-s^2/3+(1/s-1)r^2+\tau^2/(3s).
\tag{2}
$$
The eight possible sign cells in which the pair disagrees have total volume at most $8M$. Hence $F\ge1-8M>0.3619$ is necessary.

## Checking the two-cut bounds exactly

In the coordinates of (1), the signed cap deficit of a row is
$$
e_i=\frac{r_i(2s_i-s_i^2-\tau_i^2-r_i^2/3)}{2s_i(1-\tau_i^2)}.
\tag{3}
$$
This is the [cap formula](affine-pyramid-offset-proof.md#cap-deficits-and-tetrahedron-bounds) with $x_i+y_i=2s_i/(1-s_i)$, $x_i y_i=(s_i^2-\tau_i^2)/(1-s_i)^2$, and $H_i=r_i/(1-s_i)$. At fixed $s_i,\tau_i$, the magnitude of (3) increases with $|r_i|$ on the domain (1), since the numerator in its derivative contains $2s_i-s_i^2-\tau_i^2-r_i^2>0$. At fixed $s_i,r_i$, the magnitude decreases with $|\tau_i|$. To certify $|r_i|\le u$ under a cap $|e_i|\le c$, evaluate (3) at $|r_i|=u$ and $|\tau_i|=\min(s_i-u,b_i-(1+b_i)s_i)$ for $u\le s_i\le b_i/(1+b_i)$. Cross-multiplication gives two cubic polynomial comparisons with $c$. If both are nonnegative, any $|r_i|>u$ violates the cap by strict monotonicity.

For pair row bounds $a\le b$, the mean variables in (2) satisfy $\tau,r\ge0,\ \tau+r\le s$ and $s\le \tfrac12(a/(1+a)+b/(1+b))$. Write $B=(a+b)/2$ and $s_a=a/(1+a)$. Further bounds from (1) are
$$
\tau\le B-(1+a)s,\qquad
\tau\le B+(b-a)s_a/2-(1+b)s.
$$
Indeed $a s_0+b s_1\ge2a s$ and $a s_0+b s_1=2bs-(b-a)s_0\ge 2bs-(b-a)s_a$. The certified individual caps give an upper bound on $r$ by the average of the two radii. For fixed $s$, (2) is convex in $\tau,r$, so its maximum on this polygon is at a vertex. The intersections of two bounding lines give affine vertex formulas in $s$, with feasible intervals determined by the other linear inequalities. The checker verifies the cubic $s(0.3619-F)$ on each interval by rational Bernstein coefficients. This proves the first line of the table.

## Two reciprocal pairs

Suppose the four rows have two reciprocal pairs and all nonpivot scores at most $b=347/500$. Label the pairs $01$ and $23$. Formula (3) and the cap bounds give
$$
-0.0736\le r_i\le 0.1765.
$$
The two-cut checker proves $F\le0.3619$ if the average $s$ in either pair is at most $0.303$, or if the average $r$ is at most $0.07525$. Hence each pair has sum of its two $s_i$'s greater than $0.606$, and absolute sum of its two $r_i$'s greater than $0.1505$. The latter sums must be positive, since two negative offsets can have total magnitude at most $0.1472$.

In $X,U,V$ coordinates, the four rows of the normalized linear scores are
$$
(s_0,1,\tau_0),\quad (s_1,-1,\tau_1),\quad
(-s_2,\tau_2,1),\quad (-s_3,\tau_3,-1).
$$
Write $D_i$ for the absolute determinant omitting row $i$. The $D_i$ give a positive left dependence. Let $A$ be the sum of $s_i$ in the pair whose two rows occur in a minor. Put $s_*=b/(1+b)$ and $\eta=2b-0.606(1+b)$. In each pair $\sum|\tau_i|\le2b-(1+b)A\le\eta$. For example, in the minor omitting row 3, add row 0 to row 1. Expanding the determinant gives magnitude differing from $A$ by at most
$$
|\tau_0+\tau_1|s_2+|\tau_2|\,|s_0\tau_1-s_1\tau_0|
\le s_*(1+\eta)\bigl(2b-(1+b)A\bigr).
$$
The same estimate holds for all four minors by swapping pairs and labels within a pair. The function $A+s_*(1+\eta)(2b-(1+b)A)$ is increasing in $A$, so the upper bound is $2s_*$ at $A=2s_*$. The lower bound follows at $A=0.606$. The rational comparisons give $0.404\le D_i\le0.82$.

As in the [central volume formula](affine-pyramid-offset-proof.md#cap-deficits-and-tetrahedron-bounds), but now with density $1/4$ in $X,U,V$,
$$
t_{\rm cen}=\frac{(\sum_i D_i r_i)^3}{24\prod_i D_i}.
$$
If all $r_i$ are nonnegative, bound each by $0.1765$. If any $r_i<0$, replace one by zero and bound the other three by $0.1765$. Then maximize over endpoint choices for the $D_i$. This is valid by the same one-weight derivative argument as in the center-pyramid proof. When one offset is replaced by zero, its weight can be chosen as the smallest for each fixed number of upper endpoints.

For the opposing bound, (3) gives $e_i/r_i\le1$ and
$$
\frac{e_i}{r_i}\ \ge\ k_*:=1-\frac{s_*}{2}-\frac{0.1765}{6(1-s_*^2)}.
$$
For the latter, the centered slope is
$1-(x_i^2+y_i^2)/(2(x_i+y_i)(1+x_i)(1+y_i))\ge1-s_*/2$.
The cubic correction to the slope is at most $0.1765/[6(1-s_*^2)]$, using $|r_i|\le s_i$ and $|\tau_i|\le s_*$. The sum of the four $r_i$ is greater than $0.301$. If a negative offset occurs, its magnitude is at most $0.1765-0.1505=0.026$, since its pair has positive sum greater than $0.1505$. There are at most two negative offsets. Thus $\sum_i e_i$ is bounded below by $k_*\sum_i r_i-(1-k_*)\sum_{r_i<0}|r_i|$. The deficit condition $\sum_i e_i\le2G+3t_{\rm cen}$ forces a lower bound on $t_{\rm cen}$ in either sign case. The exact rational comparisons are

| Offset signs | Upper bound on $t_{\rm cen}$ | Required lower bound |
| --- | ---: | ---: |
| All nonnegative | $0.0362928854$ | $0.0373589034$ |
| At least one negative | $0.0182821868$ | $0.0331956043$ |

Both comparisons are impossible, proving the second line of the first table.

## Reproduction

From this directory run

~~~sh
python3 -B reciprocal_pivot_certificate.py
~~~

The [certificate](reciprocal_pivot_certificate.py) uses exact rational arithmetic, including the checks at the cubic Bernstein coefficients and the central volume comparisons. The cap and central volume formulas and the agreement bound were also cross-checked against exact cell volumes on 320 rational arrangements. A potential counterexample in pattern $1032$ must have a nonpivot score above $0.694$ and another below $0.70$. If the largest nonpivot is at most $0.93$, some nonpivot must be below $0.45$ by the [earlier center-pyramid result](affine-pyramid-offset-proof.md).

## 两对互选的完整仿射证明

设 $q=79759/10^6$。对于最大列排列 $1032$，八个非最大分数在 $(0,1]$ 内独立取值，四个偏移也可以在严格交替顶点类型允许的范围内独立变化。此时总有 **$M>q>m$**。

假设 $M\le q$。以全正块为空块定向，记 $t$ 为全负块的体积。十四个真符号块的亏量为 $d_J=q-V_J$，单刀亏量为 $e_i=1/2-\Pr(i+)$，并记 $G=7q-1/2$。有
$$
0\le t\le T=0.008452,\qquad 0\le d_J\le q,\qquad
\sum_Jd_J=2G+t,\qquad e_i=\sum_{J\ni i}d_J-G.
$$
令 $D_{ij}$ 为两刀的分歧体积，$D^0_{ij}$ 为将两刀分别平移过中心所得的分歧体积，则 $D_{ij}=8q-\sum_{|J\cap\{i,j\}|=1}d_J$。单位法向与偏移的全局证书给出 $D^0_{ij}\le U=0.654538$。另一方面 $D_{ij}\ge1-6q-T$，将两平面移过中心的总变化至多 $|e_i|+|e_j|\le2(G+T)$，所以 $D^0_{ij}\ge L=0.379464$。[全局界](affine-frame-height-proof.md#用中央四面体的棱长收紧候选阈值)

将两刀的 $++$ 区域和反射后的 $--$ 区域取 Minkowski 平均。所得集合包含于居中的 $++$ 区域，由 Brunn–Minkowski 不等式及两个原区域的体积差 $e_i+e_j$，得到下面的上界。将一把刀反向得到下界。[立方根不等式的推导](affine-triple-radial-proof.md#相反符号区域的中点)
$$
D_{ij}-\frac{(e_i-e_j)^2}{3L}\ \le\ D^0_{ij}\ \le\
D_{ij}+\frac{(e_i+e_j)^2}{3(1-U)}.
$$

四刀平移过中心后，每一对对顶块的体积相等。选其中七个代表，记四个单正号块为 $s_i$，三个双正号块为 $p_{0j}$。直接在十六个符号字上计数得
$$
s_i=\frac12-\frac14\sum_{a<b,\ a,b\ne i}D^0_{ab},\qquad
p_{0j}=-\frac12+\frac14\sum_{|\{a,b\}\cap\{0,j\}|=1}D^0_{ab}.
$$
对任一代表子集，先合并六个 $D^0_{ab}$ 的系数，再按符号取上述上下界。所得上界为亏量的线性函数加若干个非负倍数的平方，故为凸函数。由于 $q<2G\le2G+T<2q$，在固定的 $t$ 下，截断单纯形的每个顶点恰有一项为 $q$、另一项为 $2G+t-q$，其余为零。每对固定的非零指标从 $t=0$ 到 $T$ 给出一条线段，凸函数只需检查两个端点。由此得到 128 个有理上界，分别约束居中切法的七个代表体积的子集和，适用于原切法的全部独立偏移。[阈值的整数生成器](affine_centered_moment_certificate.py)

在 $1032$ 的八维方向盒上，用原居中证明的共享分母 Bernstein 系数计算同一子集和的下界。若下界严格超过新上界，即排除该盒的所有偏移。全部非最大分数至少 $7/10$ 的方向盒沿用六棱锥定理。完整二叉树共 **683 个节点、342 个叶子、最大深度 28**。其中 70 个叶子由单个代表、271 个叶子由多个代表的和排除，最后 1 个叶子由六棱锥定理排除。余下参数盒为零。[居中的 Bernstein 下界](centered-bernstein-proof.md#a-small-rational-bernstein-bound-for-a-row) · [六棱锥定理](affine-pyramid-offset-proof.md#all-pivot-patterns-near-the-regular-directions)

~~~sh
c++ -O3 -std=c++20 centered_chart_cover.cpp -o /tmp/cube-centered-moments
python3 -B affine_frame_height_certificate.py --candidate --edges
python3 -B affine_pyramid_matching_certificate.py
python3 -B affine_centered_moment_certificate.py | /tmp/cube-centered-moments --affine-moments 1032 32 1000000
~~~

这个排除还可用于按**平移过中心后的方向**分类。假定平移后的方向在 $1032$ 的严格角星域内，而实际偏移任意。在阈值假设下，所有外法向的正线性关系唯一；角星方向已有一个严格正关系，因此与外法向的定向一致，至多同时反向。中心反射可取规定的中央四面体朝向。上面的 341 个体积叶子只用居中体积及外法向的亏量约束，故可原样使用。余下的六棱锥叶子满足所有分数至少 $7/10$。行归一化后 $|n_{ik}|$ 是两个至多为一的分数的均值，故 $\|n_i\|_2\le\sqrt{3}$。在边长为二的坐标下，由物理边长为一时的距离界，
$$
|H_i|\le 2\sqrt{3}\,(0.066765)<7/10.
$$
故这个叶子内的实际四刀仍满足严格交替顶点条件，六棱锥定理适用。于是仅由平移过中心后所属的 $1032$ 方向图，已经能排除所有偏移。[平面到中心的距离界](affine-frame-height-proof.md#用中央四面体的棱长收紧候选阈值)

<a id="affine-four-cycle"></a>

## 四个最大项首尾相接的完整仿射证明

若四刀的实际顶点符号严格交替，最大项可以依次取在第 $1,2,3,0$ 列，则整个方向域及四个允许偏移满足 $M>q$，其中 $q=79759/10^6$。各行中并列最大的项也允许选为主项。此处使用实际顶点条件 $|H_i|<\min_{j\ne i}C_{ij}$。

在前节的亏量端点界以外，可按方向收紧两刀平移的误差。对刀对 $i,j$，在第 $i$ 行记 $g=C_{ij}$，其余两项为 $b,d$，并令
$$
x_i=H_i/(b+d),\qquad s_i=(b+d)e_i/H_i,\qquad B_i=\frac{g(b+d)}{(g+b)(g+d)}.
$$
零偏移处的比值按极限取值。单刀体积公式给出
$$
s_i=1-\frac{g^2+H_i^2/3}{(g+b)(g+d)}.
$$
第 $j$ 行同理。在[两刀的精确分歧式](affine-vertex-volumes-proof.md#两个交替角星平移后的精确分歧)中，两个投影绝对斜率和记为 $A\ge B>0$。分歧的平移差可以写成
$$
\Delta=D^0-D=B_i x_i^2+B_j x_j^2-k(x_i-x_j)^2 .
$$
这里 $k$ 是 $[J(u)-J(0)]/u^2$，$u=x_i-x_j$。切割条件给出 $2|u|\le A+B$。对梯形密度作积分有 $2/(3A)\le k\le1/A$，且 $k\ge(1-|u|/(3B))/A$。若 $2|u|+B\le A$，则 $k=1/A$。这些不等式也可直接代入 $J$ 的两段三次公式验证。

选择 $k$ 下界与 $B_i,B_j$ 上界，得到 $\Delta$ 的二次型上界 $Q(x)$，其对角项为 $a_i=B_i-k,a_j=B_j-k$，非对角项为 $k$。为验证 $Q(x)\le K(s_i x_i+s_j x_j)^2$，检查两阶对称矩阵的两个对角元及行列式非负。行列式中 $K$ 的系数为
$$
-a_j s_i^2-a_i s_j^2+2k s_i s_j .
$$
在 $s_i,s_j>0$ 的矩形中，依次固定一个变量后，另一个变量或者给出凹函数，或者给出单调递增函数。因此该系数的下界只需检查四个顶点。用 $k$ 上界和 $B_i,B_j$ 下界，对 $-\Delta$ 与 $K(e_i-e_j)^2$ 作同样比较。若这些方向系数带来的界较弱，保留中点不等式。由所得六对系数重新核验前节同一组亏量端点。

对还有可能的方向盒，再计算[三刀体积差](affine-vertex-volumes-proof.md)给出的七个 $O_r$、中央四面体 $t$ 及六个 $\Delta_{ij}$ 的向外区间。七组约束 $Y_r+|O_r|\le2q$ 和四个边际等式 $\sum_r\sigma_{ri}O_r=t-2e_i$ 同时使用。对非负权重加权的十四个块，先合并居中体积的共享行系数及偏移系数，再以有理数检查下界；权重的数值选取只影响核验效率。偏移盒的细分覆盖整个可行区间。方向达到设定深度后仍有未决盒时，程序报告未完成。[方向和偏移核验](affine_alternating_chart_cover.cpp) · [两刀和三刀区间公式](affine_alternating_moments.hpp)

循环重排四刀后，可取第零行指向循环对面一列的分数最小，即同时要求 $C_{02}\le C_{13},C_{20},C_{31}$。这个对称性保留了所有几何构型。方向覆盖共检查 **13,625,249 个节点**，最大深度 44。普通方向及居中体积比较结束 5,851,346 个叶子，联合偏移约束结束 961,279 个叶子，合计 6,812,625 个。检查了 77,662,732 个偏移节点。不存在余盒。这个类型不能进入候选局部邻域，因为候选方向中有一列在其他三行都严格小于该行的最大值，四循环的最大列却各出现一次。

~~~sh
c++ -O3 -std=c++20 affine_alternating_chart_cover.cpp -o /tmp/cube-four-cycle
python3 -B affine_centered_moment_certificate.py | /tmp/cube-four-cycle
~~~

<a id="affine-candidate-charts"></a>

## 含候选方向的最大列排列

对于下表中的四种排列，任意允许的八个方向参数和四个独立偏移均有 $M\ge m$。只有候选构型及其立方体对称、刀号重排可以取得等号。这里须把含等号的方向单独保留。令 $q_*=\lceil 2^{40}m\rceil/2^{40}$，可由 $(48m-1)^2=8$ 与 $48m>1$ 用整数算出。方向盒完全进入[半径 $1/35$ 的局部定理](local-candidate-proof.md#angular-radius-35)时直接采用该定理；其余成功叶子均排除了 $M\le q_*$。用于全局必要条件的较大有理阈值仍为 $q=79759/10^6$。

**居中体积的方向插值。** 记 $D^0(n,n')$ 为两个居中刀的分歧体积。立方体中心与六个面的锥体给出其精确值：在每个面裁出第一刀为正、第二刀为负的多边形，面积求和再除以三。归一到单位立方体，各面至中心的距离为 $1/2$，相反符号块等体积，因而得到这个因子。线性地改变一把刀 $n\mapsto n+zv$。选择 $n_k\ne0$，把其余两个坐标记为 $X$，并置 $p=n_\perp/n_k$，$a=n'_\perp-n'_kp$，$b=v_\perp-v_kp$。沿主轴积分给出
$$
 \partial_z D^0=-\frac1{|n_k|}\int_{[-1/2,1/2]^2,\ |p\cdot X|\le1/2}
          (b\cdot X)\operatorname{sgn}(a\cdot X)\,dX.
$$
在整个正方形上的积分由 $\mathbb E|a\cdot X|$ 的导数求得。超出中心截面的两个角三角形互为反射，在三角形顶点检查符号；同号时线性函数的积分等于三顶点平均乘以面积。

对固定 $z$ 坐标再次求导。来自系数的项为 $-2(v_k/n_k)\partial_zD^0$。来自两刀相交直线的项为
$$
 \frac{2n'_k}{n_k}\frac{(v\cdot w)^2}{12\|w\|_\infty^3},
 \qquad w=n\times n'.
$$
这是沿 $x=tw$ 在 $|t|\le1/(2\|w\|_\infty)$ 上积分 $(v\cdot x)^2$。剩下两项来自投影边界 $p\cdot X=\pm1/2$。例如正边界一项在方向符号确定后，用沿截线的 $(b\cdot X)^2$ 平均乘以投影密度，再乘 $1/(n_k|n_k|)$ 并取相应符号。线段两端的值若为 $l,r$，平方的平均为 $(l^2+lr+r^2)/3$。符号不确定时包围两种符号；若不能证明叉积远离零，只使用一阶界。[导数和六面裁切的整数实现](affine_angular_interpolation.hpp)

对一个长度为 $d$ 的方向区间，用 $I f$ 表示两端的线性插值。若 $L\le f''\le U$，则
$$
 -\max(0,U)d^2/8\le f-I f\le -\min(0,L)d^2/8 .
$$
同时，一阶导数变化范围为 $R$ 时有 $|f-If|\le Rd/4$。各方向的插值都是凸平均，逐坐标插值只需累加这八个误差；两个刀方向以外的误差为零。于是，七个居中代表体积的任意非负线性组合，均可在同一个八维盒的 256 个角点取下界。应先合并该组合中同一个 $D^0_{ij}$ 的系数，再扣除插值误差。用独立有理函数在分支内求导，另外核对了 15,600 个一阶导数包围、15,596 个二阶导数包围和 12,350 个加权方向下界。区间证明由插值公式和整棵覆盖树给出。

**四个偏移的共同检验。** 记 $\Delta_{ij}=D^0_{ij}-D_{ij}$，$Y_r=V_r+V_{\bar r}$，$O_r=V_r-V_{\bar r}$。单正代表取 $b_r=1$，双正代表取 $b_r=-1$，沿用[局部证明式 (2)](local-candidate-proof.md#angular-radius-35) 的 $a_{r,ij}$，则
$$
 12 V_{r,\pm}=12V_r^0-6b_rt-3\sum_{i<j}a_{r,ij}\Delta_{ij}\ \pm 6O_r.
$$
$V_{r,+}$ 和 $V_{r,-}$ 分别指这对块中被选作代表的一块和其对顶块。七个居中值的另一个下界来自先前的面三角形和棱重叠：把行分母及分子升至双三次 Bernstein 形式，将重叠面积的共同线性项并入同一行。除去零分母控制点后，行值位于相应控制值的凸包。可以在这个下界和上面的角点插值下界中取较强者。

对每一对刀，已有 $-K_-(e_i-e_j)^2\le\Delta_{ij}\le K_+(e_i+e_j)^2$。单正块的 $O_r+t$ 是另外三把刀的 $p_{rj}e_j$ 之和，其中 $p_{rj}$ 是[截面同号系数](affine-one-and-zero-edge-proof.md#单刀截面上的四个互斥事件)。程序求导三刀体积差的分段三次式，沿共同偏移射线积分各导数，并将系数包围在 $[0,1]$ 内。对固定的真实构型，可把这些系数以及 $e_j/H_j$、正核权重、中央四面体的分母固定，在辅助偏移端点上使用同一批系数。这样三刀奇部中的线性项保持原式，$t$ 的上界是线性函数正部的三次方除以正分母。

把偏移盒按四个符号分为十六个正交象限，再按四个归一化绝对值的次序把各象限分为至多 $4!$ 个单形。对每一单形寻找十四个块的非负权重。在加权和中给 $\Delta$ 取不利一侧，给 $t$ 取相应的凸上界或零下界，所得偏移下界为线性项减去若干凸函数。在单形五个顶点核验严格超过 $12q_*$ 乘以权重和，即覆盖该单形内的所有允许偏移。各单形可以有不同权重。权重的提议允许使用浮点数，接受与否只取决于非负整数权重经向外舍入后的严格整数比较。还可直接细分偏移盒，用七个 $Y_r+|O_r|$ 与四个边际关系检验。程序达到深度或节点限制时返回未完成。[完整方向和偏移覆盖](affine_alternating_candidate_cover.cpp)

**等号邻域与覆盖结果。** 尝试把四个标签分别作为候选构型中对称的那一行。若选中标签 $o$，该行先把绝对系数之和归一到三，另外三行把主分量归一到 $\sqrt2$。局部坐标要求的分量比限制都是线性分式不等式，分母为正。只须检查各行的八个角点；含 $\sqrt2$ 的正量比较通过平方化为整数比较。在 $1000$ 中交换标签二、三后还可设 $C_{02}\le C_{03}$。这个交换由立方体对称实现，候选邻域内的等号由局部定理给出。

| 最大列排列 | 方向节点 | 偏移节点 | 最大方向深度 | 未决盒 |
| --- | ---: | ---: | ---: | ---: |
| $1000$ | 1,932,811 | 6,394,866 | 50 | 0 |
| $1001$ | 3,195,195 | 8,938,300 | 64 | 0 |
| $1002$（直接细分） | 6,613,671 | 28,933,407 | 66 | 0 |
| $1002$（截面列和） | 2,033,224 | 6,929,879 | 64 | 0 |
| $1200$（64 个分区之和） | 12,370,888 | 40,147,377 | 68 | 0 |

本轮全量重算从同一根域验证了上述排列。开启 UBSan 后，1001 为 3,194,576 个方向节点和 8,938,300 个偏移节点，1002 的截面列和模式为 2,031,702 和 6,924,946。各根域未决数均为零。表中保留较早运行的数据。[逐排列实测与编译条件](../evidence/README.md)

在 $1000$ 的 966,406 个叶子中，554,290 个由方向比较或已有定理覆盖，412,116 个由联合偏移约束覆盖，后者包括 401,422 个用上述加权单形检验结束的叶子。方向节点数恰为叶子数的两倍减一。

~~~sh
c++ -O3 -std=c++20 affine_alternating_candidate_cover.cpp -o /tmp/cube-candidate-charts
python3 -B affine_centered_moment_certificate.py | /tmp/cube-candidate-charts 1000
~~~

### 排列 $1001$ 的亏量与方向联合比较

该排列表中包括 908,638 次同深度区间收紧。其余方向节点形成有 1,143,279 个叶子的二叉树，490,395 个叶子由方向比较或前置定理结束，652,884 个叶子由联合偏移约束结束。

对给定的偏移盒，将四个符号区域亏量的边际和、七个对顶块总体积及四个三刀奇部同时写在前述 364 个亏量端点上。具体地，令 $d_J=q-V_J$；对四个单正代表 $r$ 有 $O_r=d_{\bar r}-d_r$ 及 $O_r+t=\sum_{j\ne r}p_{rj}e_j$，这里 $j\ne r$ 按单正刀的标签理解。七对总代数式来自上节的 $Y_r$。分别给七个偶部约束、八个 $e_j$ 上下界和八个奇部上下界取非负权重。在每一列上使用 $\sum_{i\ne j}p_{ij}=1+2r_j$ 和 $p_{ij}\ge r_j$ 联合检验各奇部，三维系数多面体只需用整数枚举顶点。这里的 $r_j=t\lambda_j/(B\rho_j)$ 与所选法向及门槛的尺度无关，$r_j\le0.224$ 可一律用最大法向分量归一后的截面三角形证明。[截面分解](affine-one-and-zero-edge-proof.md#单刀截面上的四个互斥事件)

对每个 $\Delta_{ij}$ 的平方上界，在已知的 $e_i\pm e_j$ 区间之外用端点处的切线延伸。这保留真实区间内的值，并使不利误差在整个亏量多面体上仍为凸函数。因而加权约束的下界可由 364 个端点计算。方向项用同一方向盒的角点插值计算；将一个角坐标固定在两端，可在两张面上给出整个加权式的下界 $L,R$。设该坐标在全长中的相对位置为 $x$，则真实加权式至少为 $(1-x)L+xR$，因为插值本身在这一坐标上线性，且误差和偏移项的下界对整盒有效。该下界为正的一段可从方向区间删除。所有收紧端点都向外取整。[亏量与方向证明](affine_deficit_projection.hpp) · [偏移单形比较](affine_offset_simplex.hpp)

~~~sh
python3 -B affine_centered_moment_certificate.py | /tmp/cube-candidate-charts 1001 96 200000000 columns
~~~

### 排列 $1002$ 的两次全域证明

直接细分的覆盖有 3,306,836 个叶子。第二个证明将截面列和也代入每个带符号体积。设 $L_i=\sum_{j\ne i}p_{ij}e_j$，则由 $\sum_j r_j e_j=t$ 及 $\sum_{i\ne j}p_{ij}=1+2r_j$，恰有
$$
 2t=\sum_i L_i-\sum_j e_j,\qquad O_{2^i}=L_i-t.
$$
代入七对块的偶部与奇部后，除六个两刀平移误差外，每一项都是 $e_j$ 的线性式，系数是同一列 $p_{ij}$ 的仿射函数。在真实构型上固定这些系数，沿用偏移单形顶点比较；对每组权重及辅助端点的 $e_j$ 符号，在截面列多面体上取最不利的值。两个合法下界可以逐组权重择强使用。独立的有理体积计算另行核对了 1,440 个射线系数及 2,760 次加权体积比较。第二棵树包括 487,743 次区间收紧及 772,741 个叶子；其中 331,967 个由方向比较或前置定理结束，440,774 个由联合偏移证明结束。

~~~sh
python3 -B affine_centered_moment_certificate.py | /tmp/cube-candidate-charts 1002 96 200000000 plain
python3 -B affine_centered_moment_certificate.py | /tmp/cube-candidate-charts 1002 96 200000000 flow
~~~

表中的整数权重须通过向外舍入的有理不等式；不同的编译选项可以提议不同的权重，进而改变覆盖树。$1000$ 一行采用开启整数未定义行为检查的完整复算，原始优化运行得到 1,932,859 个方向节点，同样没有未决盒。$1001$ 随目录源文件重新编译后又完整复算了一次，得到 3,194,576 个方向节点、8,938,300 个偏移节点和零个未决盒。


<a id="affine-1200"></a>

## 排列 $1200$ 的全方向、全偏移证明

排列 $1200$ 表示四行的最大列分别为 $1,2,0,0$。每行把所选最大分数归一到 1，按行号、列号递增排列另外八个分数 $u_0,\ldots,u_7$。它们的初始闭区间均为 $[0,1]$。严格交替符号要求各分数为正且 $|H_i|<\min_{j\ne i} C_{ij}$。分区 $s=0,\ldots,63$ 的六个二进制位从低到高依次指定坐标 $0,2,4,6,3,7$ 取下半区间 $[0,1/2]$ 或上半区间 $[1/2,1]$。这样便覆盖八维全域，包括最大分数并列和分区交界处。

各区从第六层开始应用前述方向插值、偏移单形及亏量投影。四个偏移始终保留满足顶点符号和体积假设的全部可能取值。每次把方向盒拆成两个闭子盒，或用必要不等式向外舍入地收紧它。叶子或者进入半径 $1/35$ 的闭方向邻域，或者排除该盒中满足 $M\le q_*$ 的所有严格交替构型。达到预算或深度上限而尚有余盒时，程序以失败退出。[逐区覆盖程序](affine_alternating_candidate_cover.cpp) · [局部定理](local-candidate-proof.md#angular-radius-35)

| 核验域 | 完整根域数 | 方向节点之和 | 偏移节点之和 | 最大方向深度 | 未决盒 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 全部 64 个分区 | 64 | 12,370,888 | 40,147,377 | 68 | 0 |
| 第 47 区另行细分 | 64 | 1,104,757 | 4,867,919 | 60 | 0 |
| 第 63 区另行细分 | 64 | 1,616,895 | 5,202,720 | 44 | 0 |

第一行已覆盖整个排列。另一次从未分区的根盒出发的亏量投影核验也闭合，得到 10,335,318 个方向节点、36,939,111 个偏移节点、零个未决盒。分区证明的 4,645,488 个叶子中，2,598,484 个由方向界或前置定理处理，2,047,004 个由联合方向与偏移约束处理。另有 3,079,976 次同深度区间收紧。逐区记方向节点数为 $N$，叶子数为 $L$，同深度收紧次数为 $S$，则每一行记录均满足 $N=2L-1+S$。连接 64 个根需补入 63 个二叉分割节点。由此所得整棵树在略去单后继节点后有 $2\times4,645,488-1=9,290,975$ 个节点。[逐区整数记录](affine-alternating-sector-records.txt) · [记录与分类的核对](check_classification_records.py)

后两行分别在父区间的中点，再按坐标 $1,5,0,2,4,6$ 依次二分。在此复核中利用七对块的体积差收紧单刀偏差。用 $O_r$ 表示下标 $r=0,\ldots,6$ 的代表块与其对顶块之差，前四项是单正代表，后三项为掩码 $R_r=3,5,9$ 的双正代表。令 $L_i=\sum_{j\ne i}p_{ij}e_j$，沿用上一节的 $2t=\sum_i L_i-\sum_j e_j$。则 $O_i=L_i-t$，且 $2O_r=\sum_i\psi_{ri}(e_i+L_i)$，其中 $\psi_{ri}=1-2\mathbf1_{\{i\in R_r\}}$，成员资格按掩码的置位解释。对任意七个整数 $y_r$，令
$$
S=\sum_{i=0}^3y_i,\qquad
\beta_j=S+\sum_{r=4}^6\psi_{rj}y_r,\qquad
\alpha_i=2y_i-S+\sum_{r=4}^6\psi_{ri}y_r .
$$
代入上述等式得到
$$
2\sum_{r=0}^6 y_r O_r
=\sum_{j=0}^3\left(\beta_j+\sum_{i\ne j}\alpha_i p_{ij}\right)e_j .
$$
对七个下标的子集 $A$，有 $\sum_{r\in A}|O_r|\le 2|A|q_*-\sum_{r\in A}Y_r$，也有 $\sum_{r\in A}|O_r|\le1-t-\sum_{r\notin A}Y_r$。对 $Y$ 的和使用已证下界，得到上界 $U_A$。任意超集的上界以及把 $A$ 分成两组后两个上界的和都可收紧 $U_A$。把 $|y_r|$ 的不同正值按 $0=v_0<v_1<\cdots<v_k$ 排列，就有
$$
\left|\sum_r y_r O_r\right|
\le\sum_{a=1}^k(v_a-v_{a-1})U_{\{r:\ |y_r|\ge v_a\}} .
$$
式右只含已证上界。对每个 $j$，在截面列的可行多面体上界定其 $e_j$ 的系数。若该系数区间与零分离，将另三项移到另一侧再作区间除法，便得到 $e_j$ 的新必要区间。浮点求逆用于提议整数 $y$，最后的每一次比较独立使用这些整数验证。[奇部与单刀偏差核验](affine_odd_marginal_bounds.hpp) · [截面列的顶点](affine_section_columns.hpp)

在本目录执行下面的命令即可复算。前一组循环给出全域证明，后一组复算表中的两个子区族。可以按编号分别运行。

~~~sh
c++ -O3 -std=c++20 affine_alternating_candidate_cover.cpp -o /tmp/cube-candidate-charts
python3 -B affine_centered_moment_certificate.py > /tmp/cube-alternating-thresholds.txt || exit 1
for sector in $(seq 0 63); do
    /tmp/cube-candidate-charts 1200 90 120000000 project "$sector" < /tmp/cube-alternating-thresholds.txt || exit 1
done

for sector in 47 63; do
    for patch in $(seq 0 63); do
        /tmp/cube-candidate-charts 1200 94 40000000 inverse "$sector" "$patch" < /tmp/cube-alternating-thresholds.txt || exit 1
    done
done
~~~

由此 $1200$ 中所有允许的方向和四个偏移均有 $M\ge m$，等号只在居中候选构型的对称像处取得。四种包含等号的图与两种已严格排除的图合在一起，给出了整个严格交替顶点类型的结论。一般位置到任意切法的极限及全局等号条件见[原问题的定理](README.md#global-equality)。
