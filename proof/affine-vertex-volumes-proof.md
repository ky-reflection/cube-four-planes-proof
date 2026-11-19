# Exact volumes for independent affine offsets in the alternating vertex chart

Work in $[-1,1]^3$ with volumes normalized by eight. Use the parity vertices
$$
v_0=(1,1,1),\quad v_1=(1,-1,-1),\quad v_2=(-1,1,-1),\quad v_3=(-1,-1,1).
$$
Write the cuts as $n_i\cdot x=H_i$, and let $C_{ij}=-n_i\cdot v_j$ for $j\ne i$. The strict alternating vertex condition is precisely $C_{ij}>|H_i|$ for all $j\ne i$, since the sum of the four parity vertices is zero. Let $p_j$ be the intersection of the three cuts other than $j$. This note first gives exact volume formulas for all directions in this chart when the four triple points are inside the cube. It then proves a candidate-threshold exclusion for arbitrary independent offsets in a six-parameter family and an explicit neighborhood of its directions.

## Boundary areas and an exact affine recurrence

The evaluation matrix $E_{ij}=n_i\cdot v_j$ has positive diagonal, negative off-diagonal, and zero row sums. It has rank three and a positive left kernel vector. To see this, divide each row by its diagonal to obtain $I-P$, where $P$ is stochastic with strictly positive off-diagonal entries. A fixed right vector of $P$ is constant by the maximum principle. A stationary probability row vector exists by taking a limit of Cesàro averages and is strictly positive since $P^2$ is positive. Undo the row scaling to obtain the positive relation $\sum_i\lambda_i n_i=0$. In particular every triple is independent. Reflect through the center if necessary to choose $\sum_i\lambda_iH_i\ge0$, so the all-positive open cell is empty. At $p_j$ the evaluations other than $j$ are zero and evaluation $j$ is nonpositive.

For a cube face $F$, let $A_J(F)$ be the area fraction carrying the positive sign set $J$. The two parity-positive vertices of $F$ give the default pair of labels $B$. On the face each cut toggles its sign on a triangle at one corner. For $i\in B$ it removes a positive sign and for $i\notin B$ it adds a positive sign. For a triangle associated with the corner label $j$ and cut $i$, its two legs are
$$
\frac{C_{ij}\pm H_i}{C_{ij}+C_{ik}},\qquad
\frac{C_{ij}\pm H_i}{C_{ij}+C_{il}},
$$
where $k,l$ are the other labels and the plus sign is used for removal. Its area fraction is half their product. Same-type triangles are at opposite corners and have disjoint interiors. An adjacent removal and addition triangle overlap in a triangle of area fraction $(\ell+r-1)_+^2/(2d)$, with $\ell,r$ their legs along the common edge and $d=\ell/a+r/b$ for their perpendicular legs $a,b$. The ratio of two legs of one triangle is independent of its offset. No three triangles have overlapping interiors. Thus $A_B$ is one minus the four triangle areas plus the overlaps, $A_{B\triangle\{i\}}$ is the triangle area for cut $i$ minus its overlaps, and $A_{B\triangle\{i,j\}}$ is their overlap for opposite types. All other area fractions are zero.

For a function $A$ on the six faces define
$$
\Phi_p(A)=\frac16\sum_F (1-\nu_F\cdot p) A(F),
$$
where $\nu_F$ is the outward unit normal of the face. This is the volume of the pyramids from $p$ to the specified face patches. Fix a nonempty proper sign set $J$, and choose any $j\notin J$. The union of the cells with signs $J$ and $J\cup\{j\}$ leaves the sign of cut $j$ free. Its other three bounding cutting planes all pass through $p_j$. Consequently, its restriction to the cube is exactly the union of the pyramids from $p_j$ to its six face patches. The resulting recurrence is
$$
V_J=\Phi_{p_j}(A_J+A_{J\cup\{j\}})-V_{J\cup\{j\}}. \tag{1}
$$
Begin with $V_{\{0,1,2,3\}}=0$, then compute triples, pairs, and singletons. Subtract these from one to get the bounded tetrahedron volume. The formula treats all offsets independently, including concurrence as a limit. It expresses each volume as a piecewise cubic in the offsets at fixed directions. The only branch tests inside this chart are the affine signs of edge-overlap differences. All directional coefficients and triple vertices are rational in $C$ and $H$. The [exact Python implementation](affine_vertex_volumes.py) was cross-checked against independent three-dimensional rational volume integration on 120 arrangements, with equality for all sixteen volumes.

## All offsets for coordinate-rescaled tetrahedral directions

Let the four normals be obtained by multiplying the coordinates of $v_i$ by common positive factors $a,b,c$. Reorder the axes so that $a\le b\le c$. The alternating vertex chart requires $c<a+b$. Normalize the largest off-diagonal magnitude to one and put
$$
2a=u+v,\quad 2b=1+v,\quad 2c=1+u,\qquad v=uz,\quad 0<u,z\le1.
$$
Then $C_{ij}$ equals $1,u,v$ according as $i\mathbin{\mathrm{xor}}j$ equals $1,2,3$, respectively. Write $H_i=v t_i$; the allowed offsets are exactly $-1<t_i<1$, independently for each $i$.

For these directions the four normal weights are equal. By reflection and the coordinate sign changes preserving parity, take $T=\sum_i t_i\ge0$ and $t_0\le t_i$ for all $i$. For a face with default pair $B$, the projection of $p_j$ on its outward normal is
$$
\nu_F\cdot p_j = r_B\bigg(\sum_{i\in B}t_i-\mathbf 1_{j\in B}T\bigg),
$$
where $r_B$ is $v/(u+v)$, $v/(1+v)$, or $v/(1+u)$ according to the face axis. This follows from the two scores whose sum is the relevant coordinate times its doubled scale. The expression in parentheses is a signed sum of two offsets $t_i$, and $r_B\le1/2$. Thus **every triple intersection is inside the cube** for every strictly allowed offset. In particular the proof in this section needs no three-cut lower bound.

For an overlap between cuts $i\in B$, $j\notin B$, the two longitudinal legs use the same centered numerator $A$ and denominator $A+B'$, with distinct entries $A,B',E'$ drawn from $1,u,v$. The two perpendicular legs have denominator $A+E'$. Hence the overlap area fraction simplifies to
$$
\frac14\left(\frac{A-B'+v(t_i-t_j)}{A+B'}\right)_+^2\frac{A+B'}{A+E'}.
$$
If a denominator contains both $u$ and $v$, cancel the factor $u$ in the longitudinal ratio or in a triangle leg. For the last factor, a uniform bound of four holds on positive support. Only $B'=1$ and $A,E'$ both small needs a check. Positive support then forces $A+2v>1$, hence $A+E'>1/2$. The other cases have denominator at least one.

For extra feasibility tests use the signed positive-cap deficit
$$
e_i=\frac{(u+v+uv)H_i-H_i^3/3}{(u+v)(1+u)(1+v)}.
$$
This is the elementary sum-of-three-uniforms formula on $|H_i|<v$. Inclusion-exclusion in a rectangular box gives a cumulative volume proportional to $(H_i+1+u+v)^3-(H_i+1)^3-(H_i+u)^3-(H_i+v)^3$, because exactly these four truncated-power terms are nonzero on this interval. Expanding and subtracting one half gives the displayed deficit. If every block is at most $M$, one side has at most seven pieces and the other at most eight, so $1/2-7M\le e_i\le8M-1/2$.

The [integer checker](rescaled_tetrahedral_affine_cover.cpp) uses these formulas on boxes in $(u,z,t_0,t_1,t_2,t_3)$. Every scalar interval is rounded outward on a fixed dyadic grid. The offset-order and sum tests contract a box only by necessary linear inequalities. For volume bounds, (1) is expanded before evaluation to cancel common face terms. For instance for a pair $J$ with omitted labels $i,j$, it becomes
$$
V_J=\Phi_{p_i}(A_J)+(\Phi_{p_i}-\Phi_{p_j})(A_{J\cup\{i\}}).
$$
The face multiplier in a difference is $r_BT(\mathbf 1_{i\in B}-\mathbf 1_{j\in B})$. For singletons there is one more term of the same kind. The checker tries all choices of omitted labels. A box is accepted only when the cap constraint fails throughout it or one entire cell has a lower bound above the threshold. Binary subdivision covers every point not already excluded by the offset contractions.

The direct run at $m=(1+2\sqrt2)/48$ has no unresolved or pending boxes. Thus **none of these coordinate-rescaled tetrahedral directions can give a counterexample, for any of the four independent offsets preserving the vertex signs**. This includes nonconcurrent arrangements with 15 positive-volume cells.

## Allowing unequal directions with all four offsets still free

The same checker with a one-thousandth reserve proves a uniform margin above $m+1/1000$ in the rescaled family. Use it to cover directions that are no longer exact coordinate rescalings. In pivot chart $1032$, let $u_*$ be the largest of the four entries $C_{i,i\mathbin{\mathrm{xor}}2}$ and $v_*$ the largest of the four entries $C_{i,i\mathbin{\mathrm{xor}}3}$. If
$$
\sum_i\bigl(u_*-C_{i,i\mathbin{\mathrm{xor}}2}+v_*-C_{i,i\mathbin{\mathrm{xor}}3}\bigr)\le1/1000,
$$
then there is again no counterexample at any offsets preserving the vertex signs. In particular one may allow a spread of $1/8000$ in each of the two kinds of nonpivot entries.

Compare with the rescaled normals having common entries $1,u_*,v_*$ and the same offsets. Each reference entry is at least its actual entry, so the reference stays in the alternating vertex chart. Interchange the two larger axes if necessary to put the reference in the order used by the checker. The sup norm of the change of each linear score over the cube is at most the sum of the changes of its two nonpivot entries. This follows from $n_i=\sum_{j\ne i}C_{ij}(v_i-v_j)/4$ and $\|v_i-v_j\|_1=4$. Each reference score has marginal density at most one, by integrating along a coordinate with normal component of magnitude at least $1/2$. The loss in any chosen sign cell is at most the sum of four one-sided strip volumes, hence at most $1/1000$. The reference margin therefore proves the strict candidate-threshold claim for the actual cuts.

## Reproduction

```sh
c++ -O3 -std=c++20 rescaled_tetrahedral_affine_cover.cpp -o /tmp/affine-cover
/tmp/affine-cover 48 3000000
/tmp/affine-cover 54 5000000 1000
```

The optional third argument is the reserve in millionths, added with outward rounding. The two runs visit 164,759 and 241,741 nodes, respectively. Both end with zero unresolved and zero pending boxes. The interval lower bounds were also cross-checked against exact volume evaluation at 427 rational configurations. The [column-weight proof](column-weight-affine-proof.md) covers another direction family with all offsets free and shows that the strict alternating vertex condition alone puts every triple point inside the cube.

## A box with eight independently varying direction coordinates

Normalize one maximal off-diagonal entry in each row to one. In the following matrix keep the four entries equal to one fixed and vary each of the other eight off-diagonal entries independently by at most $23/1000$.
$$
C_*=
\begin{pmatrix}
0&1&.620&.850\\
1&0&.550&.930\\
.820&.630&0&1\\
.580&.910&1&0
\end{pmatrix}.
$$
Every such direction tuple and every set of independent offsets preserving its strict alternating vertex signs has a block with volume strictly above $m$. At $C_*$ the rows are neither column-constant up to rescaling nor the coordinate-rescaled tetrahedral directions, as seen already from the first two rows.

For general scores $C$, let $E$ be the evaluation matrix, with diagonal $R_i=\sum_{k\ne i}C_{ik}$ and off-diagonal $-C_{ik}$. Write a triple point in affine coordinates $p_j=\sum_i\alpha_i v_i$ with $\sum_i\alpha_i=1$. For $i\ne j$, put $y_i=\alpha_i-\alpha_j$. Then $p_j=\sum_{i\ne j}y_i v_i$, and $E^{(j)}y=H_{-j}$, where $E^{(j)}$ is the principal minor deleting $j$. Its determinant $D_j$ is the sum of the 16 directed tree products rooted at $j$. Indeed, expanding $\prod_{i\ne j}R_i$ sums all choices of an outgoing edge from each nonroot, and the other terms in the determinant remove precisely the choices with a directed cycle. The entries of $(E^{(j)})^{-1}$ are ratios with this positive denominator. For a diagonal entry $i,i$, with the other two nonroots $a,b$, its numerator is $(C_{ai}+C_{aj})R_b+C_{ab}(C_{bi}+C_{bj})$. For an off-diagonal entry $i,k$, with the remaining nonroot $l$, its numerator is $C_{ik}R_l+C_{il}C_{lk}$. Both follow by expanding the corresponding cofactor. Since $E$ has zero row sums and rank three, its adjugate has identical rows $(D_0,D_1,D_2,D_3)$. These give the weights for reflecting to $\sum_j D_j H_j\ge0$.

For each row put $a,b,d$ equal to its off-diagonal scores. The same cap expansion gives the signed positive-cap deficit
$$
e_i=\frac{(ab+ad+bd)H_i-H_i^3/3}{(a+b)(a+d)(b+d)}.
$$
The checker starts with every offset independently bounded in absolute value by the least upper endpoint among its row's three scores. This encloses the entire offset domain, including when it varies with the directions. It checks the reflection and cap inequalities, then uses (1) to enclose the 14 nonempty proper sign-cell volumes jointly. In the shared affine arithmetic the four offsets are variables, while all eight direction intervals contribute to the remainder. Leg fractions are intersected with $[0,1]$, which is valid on the feasible subset of a box. Every direction interval is rounded outward, and an exclusion uses only the integer-verified comparisons described in the [column-weight proof](column-weight-affine-proof.md#joint-lower-bounds-on-the-seven-variable-boxes).

Run the [direction-box checker](affine_direction_box_cover.cpp) from this directory.

~~~sh
c++ -O3 -std=c++20 affine_direction_box_cover.cpp -o /tmp/direction-cover
/tmp/direction-cover 23 40 60000
~~~

The run visits 23,443 nodes with zero unresolved and zero pending boxes. The affine and scalar enclosures were also compared with exact volumes on 241 rational configurations in 54 direction and offset boxes. This verifies one full-dimensional box of the twelve-parameter chart. The [center-pyramid proof](affine-pyramid-offset-proof.md) covers the full direction region with normalized scores at least $7/10$ across all pivot patterns, as well as several other all-offset regions for reciprocal pivots. The [six affine pivot covers](reciprocal-pivot-proof.md#affine-candidate-charts) complete the remaining directions.

## 两个交替角星平移后的精确分歧

沿 $v_i,v_j$ 的唯一同号坐标轴投影，在边长为一的中心立方体中，两刀可写成两张截高图 $F_i,F_j$，两外侧沿竖直线的朝向相同。设 $F_r=z_r+a'_r X+b'_r Y$，令 $a_r=|a'_r|,\ b_r=|b'_r|$，以及 $l_r=(a_r+b_r-1)/2$。两图的横向斜率符号相反。在上面 $Z=1/2$，两图超出的区域是对角顶点处互不重叠的三角形，下面 $Z=-1/2$ 同理。因而两截高之差取绝对值后，上下截断所减少的期望等于四个三角形上的超高体积之和。设 $A=a_i+a_j,\ B=b_i+b_j$，以及 $J(d)=\mathbb E|AX+BY+d|$，则
$$
D_{ij}=J(z_i-z_j)
-\sum_{r=i,j}\frac{(l_r+z_r)^3+(l_r-z_r)^3}{6a_rb_r}.
$$
严格角星保证 $|z_r|<l_r$ 且 $l_r+|z_r|<\min(a_r,b_r)$。分子来自直角三角形底上的线性高度积分。函数 $J$ 只需对正方形上一条直线两侧积分，分段公式已列在[截高图分歧证明](affine-one-and-zero-edge-proof.md#两把截高图与同一坐标半立方体的分歧)。

若 $C_{ij}$、$H_i$ 使用本文边长为二的约定，令 $K_i$ 为第 $i$ 行三个非对角分数的三组两两和之积，投影轴为 $k$。有 $z_i=H_i/(2n_{ik})$。将上式减去零偏移的值，得到
$$
D^0_{ij}-D_{ij}
=\frac{C_{ij}H_i^2}{K_i}+\frac{C_{ji}H_j^2}{K_j}
-\big[J(z_i-z_j)-J(0)\big].
$$
例如令第 $i$ 行其余两项为 $b,c$，则 $l_i=C_{ij}/(b+c)$、$|z_i|=|H_i|/(b+c)$，两横向斜率的绝对值是 $(C_{ij}+b)/(b+c)$ 和 $(C_{ij}+c)/(b+c)$。代入展开式即得该行的二次项。此式允许不同角星各自选择方向和偏移。用完整十六块的有理体积，对 360 组方向及偏移逐一检查六组刀对，共 2,160 组分歧体积和平移差均相等。[精确公式与复核](affine_alternating_pair.py)

## 三刀同号体积的奇部与七对符号块

固定遗漏的刀号 $o$，记其他三把刀的正侧交集体积为 $P_o(H)$，交点为 $p_o$。记 $F_{ij}$ 为含 $v_i,v_j$ 的立方体面，$f_{ij}=\nu_{F_{ij}}\cdot p_o$。改变全部偏移的符号使 $p_o$ 取负。以下恒等式给出 $12[P_o(-H)-P_o(H)]$。

先对每个 $k\ne o$，取余下的两号 $i,j$，置 $g=C_{ko}$，累加
$$
\frac{4gH_k+2f_{ij}(g^2+H_k^2)}{(g+C_{ki})(g+C_{kj})}.
$$
这是 $F_{ij}$ 上补入第 $k$ 个正号的三角形在相反偏移下的体积差。它来自展开 $(1+f_{ij})(g+H_k)^2-(1-f_{ij})(g-H_k)^2$。含 $v_o$ 的三个面需要同时补入两个正号，而对应的两个对角三角形内部不交，故贡献为零。

还要减去移除三角形的重叠。对每个无序对 $i,j\ne o$，取第三号 $k$，置
$$
L=(C_{ik}+C_{io})^{-1},\quad R=(C_{jk}+C_{jo})^{-1},\quad
x=C_{ik}L+C_{jo}R-1,\quad y=H_iL-H_jR,
$$
$$
d_F=(C_{ik}+C_{ij})L+(C_{jo}+C_{ji})R,\qquad
d_K=(C_{jk}+C_{ji})R+(C_{io}+C_{ij})L.
$$
$F=F_{ik},K=F_{jk}$ 上这一对重叠的有向边长差分别是 $u=x+y$ 和 $-u$。置 $v=x-y$，则在相反偏移下，它们变为 $v,-v$。有 $d_F-d_K=2x$。将两个平面等式 $n_i\cdot p_o=H_i$, $n_j\cdot p_o=H_j$ 分别乘以 $2L,2R$ 再相减，沿公共棱的分量抵消，余下 $f_{ik}d_K-f_{jk}d_F=2y$。记
$$
a=(1/d_F+1/d_K)/2,\qquad b=(f_{ik}/d_F+f_{jk}/d_K)/2.
$$
将两个面的正部平方重叠公式合并，用 $2w_+^2=w^2+w|w|$ 以及上述两等式，所得累加项为
$$
4axy-2b(x^2+y^2)-\frac{|x+y|^3-|x-y|^3}{d_Fd_K}.
$$
上述三项三角形和三项重叠给出全部体积差。固定方向后，除一次项外的余项受偏移的三次方控制。例如
$$
\left||x+y|^3-|x-y|^3-6x|x|y\right|\le 2|y|^3
$$
可由 $|x|\ge |y|$、$|x|\le |y|$ 分别展开得到。

记 $O_r=V_r-V_{15\mathbin{\mathrm{xor}}r}$，其中 $r$ 是四位正号掩码。由中心反射得 $O_{2^i}=P_i(-H)-P_i(H)-V_0+V_{15}$。再令 $e_i=1/2-\Pr(i\text{ 为正})$，$\sigma_{ri}$ 在 $r$ 含 $i$ 时为 $1$，否则为 $-1$。另外三对由
$$
O_r=-\tfrac12\sum_{i=0}^3\sigma_{ri}(e_i+O_{2^i}),\qquad r=3,5,9
$$
给出。事实上，单刀符号的期望为 $-2e_i$。把十六块按相反符号配对后，和式取内积 $\sigma_r$ 时常数项消去，三个两正两负的符号向量互相正交且平方范数均为四，每个单号向量贡献 $2\sigma_{ri}$。这一推导也允许两种中央四面体的朝向。[分歧、三刀奇部与逐块有理核对](affine_alternating_pair.py)
