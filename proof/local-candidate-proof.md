# Exact local minimum at the centered candidate

This note supplies the Taylor bound and the minimax step for [local_certificate.py](local_certificate.py). Use the centered side-one cube and the twelve-parameter plane chart in the script. Denote the candidate block volume by $m=(1+2\sqrt2)/48$. All twelve blocks listed by the script as active have volume $m$ at the candidate. We show that at least one of these blocks has volume greater than $m$ whenever the parameter vector $p\ne0$ has $\|p\|_\infty\le 1/2500000$.

## Uniform Taylor bound for the active volumes

Along a unit sup-norm parameter direction, write any triple vertex system as $(A_0+tA_1)x(t)=b_0+t b_1$. Include the six cube-face equations among the ten supporting planes. Opposite cube faces cannot meet. For other triples singular at zero, the script checks that a Cramer numerator has magnitude greater than one. All entries of the determinant matrices, including those for the numerators, and their derivatives have absolute value at most two over the stated radius. A determinant therefore changes by at most $144|t|$. Such a triple cannot have a common point in the cube, since at a common point its Cramer numerators equal its determinant times the coordinates. For nonsingular triples whose zero-parameter point is outside the cube, the script checks that some coordinate has magnitude greater than $7/10$ and that $\|A_0^{-1}\|_\infty<5$. A point inside the cube that solves the perturbed equations would satisfy $\|A_0x-b_0\|_\infty\le3|t|$, placing it within $15|t|$ of the outside point, again impossible.

For every remaining triple, the script bounds $\|A_0^{-1}A_1\|_\infty\le4$ and $\|x'(0)\|_\infty\le5/2$. The exact identity
$$
x(t)=x(0)+t(I+t A_0^{-1}A_1)^{-1}x'(0)
$$
gives derivative bounds $251/100$, $201/10$, and $241$ for the first, second, and third derivatives in sup norm. For example the bound for derivative order $k\ge1$ is $(5/2)k!4^{k-1}/(1-4|t|)^{k+1}$, by the geometric series.

Between changes of incidence, triangulate the cell boundary into at most 28 triangles. Indeed a three-dimensional convex polyhedron with at most ten facets has at most sixteen vertices by Euler's formula and the minimum vertex degree three. A triangulation using these vertices has at most 28 triangles. The oriented sum of the tetrahedron determinants from the origin to the triangles equals the cell volume, even when the origin is outside the cell. Each triangle contributes at most
$$
3(241)(1/2)^2+18(201/10)(251/100)(1/2)+6(251/100)^3
$$
to the absolute third derivative. This follows by the product rule, the coordinate bound $1/2$ at a feasible vertex, and $|\det(u,v,w)|\le6\|u\|_\infty\|v\|_\infty\|w\|_\infty$, with the factor $1/6$ for volume. Multiplying by 28 gives less than $26000$.

The first two derivatives agree across the incidence changes. To see this, differentiate the volume by integrating the signed motion of each moving facet. Eliminate one coordinate on its plane using the component that stays nonzero in the chart. The resulting planar integrand is the signed affine velocity $(n_i'\cdot x-h_i')/|n_{i,k}|$. The second derivative consists of the derivative of this integrand and the boundary terms from the moving edges of the projected facet polygon. A disappearing facet collapses to a point, since collapse to a segment would require its plane and the two planes incident to that segment to form a feasible singular triple. A persistent edge likewise has a unique neighboring supporting plane. The projected normal defining an edge is bounded away from zero at its limit, since the two planes are among a feasible nonsingular vertex triple there. Thus all the projected facet and edge integrals have the same limits on both sides, with vanishing area or length integrals for disappearing features. These are the derivative formulas used by the script at zero, with the polygon integrals evaluated by triangulation and the edge integrals by integrating products of affine functions. Incidence changes occur at finitely many parameter values, because each candidate vertex is a rational function of $t$ and the face tests are polynomial after clearing denominators. Identically vanishing tests can be retained on the whole interval. This proves that Taylor's formula at zero has remainder of absolute value at most $(26000/6)\|p\|_\infty^3$ for each active block.

## Weighted Taylor comparison

Write $\delta=\|p\|_\infty$, and denote the active gradients and Hessians at zero by $J_i,H_i$. The script gives positive weights $w_i>1/21$ with $\sum_iw_i=1$, $\sum_iw_iJ_i=0$. Set $F=\sum_i w_iH_i$ and $Q(p)=p^\top Fp/2$. It verifies $\sum_{jk}|(H_i)_{jk}|<6$ and $\sum_{jk}|F_{jk}|<5/4$, with the matrices symmetrized.

If all the active volumes are at most $m$, the Taylor bound implies $J_i p\le(3+(26000/6)\delta)\delta^2$. Their weighted sum is zero, so $J_i p\ge-20(3+(26000/6)\delta)\delta^2$. The script checks an inverse norm bound of 70 for a nine-row and nine-column submatrix of $J$. Use that submatrix to find $e$ supported on its nine columns with the same selected row values as $p$. The checked arithmetic gives $\|e\|_\infty\le4300\delta^2$. Then $k=p-e$ is in the kernel of all active gradients, since the three independent displayed kernel vectors in the script give rank exactly nine.

For $k=\alpha r+\beta u+\gamma v$ in that kernel, the script gives
$$
Q(k)>\alpha^2/50+(\beta^2+\gamma^2-\beta\gamma)/25
$$
unless $k=0$. Every row of $(u,v)$ has squared dual norm at most $4b^2/3$ for the form $\beta^2+\gamma^2-\beta\gamma$, with $b=(8+9\sqrt2)/14<2$. To check this, the squared dual norm of coefficients $(u_j,v_j)$ is $4(u_j^2+u_jv_j+v_j^2)/3$. For same-sign entries their absolute sum is at most $b$, and for opposite signs the expression in parentheses is at most the greater squared entry. Both zero-offset rows and the first two parameter rows also satisfy the bound. Since the coordinates of $r$ are at most one in absolute value, Cauchy–Schwarz now gives $\|k\|_\infty^2\le (50+25\cdot4b^2/3)Q(k)\le250Q(k)$.

The script checks $4300\delta<1/500$ over the stated radius. Thus $\|k\|_\infty\ge(499/500)\delta$, and $|Q(p)-Q(k)|\le (5/8)4300(1001/500)\delta^3$ using the weighted Hessian norm. The checked final inequality
$$
\bigl[(5/8)4300(1001/500)+26000/6\bigr]\delta<(499/500)^2/250
$$
shows that the weighted sum of the active volumes exceeds $m$ for $\delta>0$. This contradicts the assumed volume bound and proves strict local minimality.

## 十二个活跃块的统一二阶证书

在同一十二参数坐标中，严格局部极小的半径可提高到 **$1/100$**。以下证明对整个闭盒同时成立。唯一的例外点 $p=0$ 正好取得 $m$。

先把[交替顶点体积公式](affine-vertex-volumes-proof.md)中的两个相邻面合并。设面 $F$ 上交叠的两角三角形来自刀 $i,j$，相邻面为 $K$，沿公共棱的腿长超额为 $\Delta$，两个面上的斜率和分别为 $d_F,d_K>0$。相邻面的超额为 $-\Delta$。某个三刀母块省去刀 $b$，顶点 $p_b$ 到两个面的边长二坐标距离为 $h_F,h_K$。若 $b=i$ 或 $b=j$，这对交叠的正负贡献直接抵消。否则 $p_b$ 位于刀 $i,j$ 上，且
$$
h_F d_K-h_K d_F=-2\Delta .
$$
事实上，将棱长归一到一，记到 $K,F$ 的内向距离为 $y,z$。在棱上两条切割线的有向间距为 $\Delta$，在 $F$ 面上随 $y$ 以斜率 $-d_F$ 改变，在 $K$ 面上随 $z$ 以斜率 $d_K$ 改变。两个切割平面各自可沿棱坐标解出一个仿射函数，所以在它们相交时有 $\Delta-d_Fy+d_Kz=0$。代入 $h_F=2z,h_K=2y$ 即得等式。

因此两座棱锥对同一个符号的交叠修正，在乘上该符号的系数前，等于
$$
\frac1{12}\left[
  \frac12\left(\frac{h_F}{d_F}+\frac{h_K}{d_K}\right)\Delta^2
  -\frac{|\Delta|^3}{d_Fd_K}
\right].
$$
这是把 $h_F(\Delta_+)^2/(12d_F)+h_K((-\Delta)_+)^2/(12d_K)$ 展开后代入上式得到的。三个面内角重新排序、或者交叠从一个面换到另一个面时，都可以用同一个二阶可微的表达式。

设 $I$ 为十二个活跃符号，取前文已验证的一阶平衡权重 $w_i^0$。在局部坐标上定义
$$
F(p)=\sum_{i\in I}\left(w_i^0+B_i p\right)\bigl(V_i(p)-m\bigr).
$$
矩阵 $B$ 的每个系数是分母 $8192$ 的有理数，全部列在[区间二阶证书](affine_local_merit.cpp)的常数表中。选择这些数旨在使二阶型正定，证明仅使用表中数值。$F(0)=0$，并且 $\nabla F(0)=0$ 由 $V_i(0)=m$ 和 $\sum_i w_i^0\nabla V_i(0)=0$ 精确得到。程序验证整个盒子上每个权重都严格为正。

计算时用全正块为空的递推公式定义 $F$ 的延拓。给定实际切法，必要时作中心反射，保持方向不变、反转四个偏移，即可使实际体积属于这一分支。原点到目标参数的线段即使穿过共点构型，也可以对同一个光滑延拓积分。几何公式在目标点给出正确的体积。

在每条从原点出发的径向线段上使用三阶 Taylor 公式。记 $A=\nabla^2F(0)$，并定义径向平均
$$
D_{ijk}(p)=3\int_0^1(1-t)^2\partial_{ijk}F(tp)\,dt.
$$
于是 $2F(p)=p^TAp+\tfrac13\sum_{ijk}D_{ijk}(p)p_ip_jp_k$。立方体棱上换序时，$|s|^3$ 的前两阶导数连续，两侧的三阶导数分别为 $\pm6$。沿射线，边换序的条件是一个有理函数为零，分母在该邻域非零；零点只有有限多个，或者该函数恒为零。分段积分即得公式；恒为零的项沿射线导数为零，可以给其外层三阶导数选零。这里的三阶张量只用来包住射线的三阶导数。

程序用 $2^{-40}$ 为单位向外舍入，传播值及前三阶偏导数，将 $[0,1]$ 等分为十六段。在第 $j$ 段用半径 $j/(16\cdot100)$ 的盒子包住 $tp$。各段权重为 $[(17-j)^3-(16-j)^3]/16^3$；加权求得 $D_{ijk}$ 的上下界。零点的矩阵 $A$ 单独计算。给定位移的坐标符号 $s_i$ 后，在 $0\le |p_k|\le 1/100$ 上，对每项 $s_i s_jD_{ijk}p_k$ 取下界。这样得到一张有理对称矩阵 $A_s$，在对应符号域中满足 $2F(p)\ge p^TA_sp$。全部 $2^{12}=4096$ 张矩阵的顺序主子式都由无舍入的整数消元证为正。故每个非零 $p$ 都有 $F(p)>0$，至少一个活跃块严格超过 $m$。[三阶区间运算](affine_local_merit.cpp) · [整数正定性检验](affine_local_merit_certificate.py)

~~~sh
c++ -O2 -std=c++20 affine_local_merit.cpp -o /tmp/cube-local-merit
python3 -B affine_local_merit_certificate.py /tmp/cube-local-merit
~~~

## 上述八维方向邻域允许所有独立平移

仅要求前八个方向坐标满足 $|p_i|\le1/100$，四个偏移可以任取。整个方向邻域仍有 $M\ge m$，等号只在候选点及等价表示处出现。

若 $M\le m$，每个单位法向的实际偏移至多为 $0.067285$。在局部坐标中四个未归一化的法向长度都小于 $2.02$，所以方程右侧的四个偏移 $|d_i|<0.14$。在边长二的立方体中使用 $H_i=2d_i$，全部方向分数至少为 $2-\sqrt2-0.02>0.28$，故整段都保持严格交替顶点符号。中心反射允许选取全正块为空的分支。[偏移界](affine-frame-height-proof.md#用中央四面体的棱长收紧候选阈值)

对这个参数盒，把每个方向分数、三刀交点和六个面上的面积表示为十二个共同变量的仿射式加有界余项。三刀交点使用克拉默公式，分母在使用前检验符号，再按面棱锥递推得到十四个体积的下界。对每个参数盒，先检验全正块为空的定向条件，然后寻找一组非负权重，使十四个体积的加权和超过同权重的 $m$。权重建议可用浮点数求，但最后的不等式用整数向下核验。若四个局部偏移也都落入 $[-1/100,1/100]$，即应用上节在精确 $m$ 处的局部定理。方向端点在数值计算中向外舍入；交给局部定理时使用方向域本身的精确假设。[仿射式的误差和权重核验](affine_box_linear.hpp)

这棵树共有 43,279 个节点。17,544 个叶子由定向或体积排除，4,096 个叶子落入局部定理，没有未决叶子。由于局部定理在非零扰动时为严格不等式，而其余叶子均排除 $M\le m$，等号结论也随之得到。[方向与全偏移覆盖程序](affine_candidate_direction_cover.cpp)

~~~sh
c++ -O3 -std=c++20 affine_candidate_direction_cover.cpp -o /tmp/cube-candidate-directions
/tmp/cube-candidate-directions
~~~

<a id="angular-radius-35"></a>

## 方向半径 $1/35$ 内的全部偏移

沿用上节的八个方向坐标 $x$ 及四个偏移。若 $\|x\|_\infty\le1/35$，则任意偏移均有 $M\ge m$，等号只在候选点及其等价表示处出现。这里的方向半径按法向分量的归一化定义，坐标顺序与上节相同。

**居中部分。** 将十二个活跃块配成六对。七个代表中省略小单正块 $r=1$（即标签零），其余每对使用共同的权重 $w_r(x)$。起始权重沿用 $w_A=1/(6\sqrt 2), w_B=(1-1/\sqrt2)/6$，八个方向的一次系数见[有理权重表](affine_angular_merit_weights.hpp)，权重在该盒内均为正。令
$$
F_0(x)=2\sum_r w_r(x)(V_r^0(x)-m).
$$
$F_0(0)=0$，且 $\nabla F_0(0)=0$ 仍由十二个活跃块的精确一阶平衡给出。计算 $F_0$ 时，仅需居中三角形和相对棱的重叠。若两分支的分母为 $d_+,d_-$、边长超额为 $\delta$，则 $d_+-d_-=2\delta$，重叠项等于
$$
\tfrac12(1/d_++1/d_-)\delta^2-|\delta|^3/(d_+d_-).
$$
这使得同一套三阶区间运算可直接用于八变量的加权和。对每个方向坐标的正负分别限制径向盒，使用上节的十六段 Taylor 积分，全部 $2^8=256$ 张下界矩阵均以整数主子式证为正定。因此 $F_0(x)>0$ 对每个 $x\ne0$ 成立。[八维计算](affine_angular_merit.cpp) · [主子式复算](affine_local_merit_certificate.py)

**平移部分。** 假设 $M\le m$，把全正侧朝向空块。上节的方向归一化下，外法向到中心的距离至多 $0.066765$，所以本盒内实际偏移仍保留严格交替顶点分法。具体地，在边长二的坐标中，最小顶点评分至少 $2-\sqrt2-2/35$，最大的法向长度至多 $\sqrt{2+2(1+1/35)^2}$ 或 $\sqrt{2(1+1/35)^2+(1+2/35)^2}$。将较大值乘以 $2(0.066765)$ 仍小于最小评分。

写 $Y_r=V_r+V_{\bar r}$，$O_r=V_r-V_{\bar r}$，$\Delta_{ij}=D^0_{ij}-D_{ij}$。对单正代表令 $b_r=1$，两正代表令 $b_r=-1$。令 $a_{r,ij}=-1$ 当单正标签不在 $i,j$ 中，令 $a_{r,ij}=1$ 当两正代表分离 $i,j$，其余为零。由七对体积的恒等式，
$$
\sum_r w_r(Y_r+|O_r|-2V_r^0)
=\sum_r w_r |O_r| - t\sum_r b_rw_r
-\tfrac12\sum_{i<j}\Delta_{ij}\sum_r a_{r,ij}w_r. \tag{2}
$$
对 $s=\max_i |e_i|>0$，置 $u_i=e_i/s$，在八张 $u_i=\pm1$ 的三维面上验证 (2) 除以 $s$ 后严格为正。其余 $u_j$ 在 $[-1,1]$ 独立取值。取到正最大值时用 $0\le s\le 7q-1/2+T$，取到负最大绝对值时上界可用 $7q-1/2$。这里 $q=0.079759,T=0.008452$。将每行顶点评分和门槛以相同常数缩放，核验使用每行评分除以二的坐标。门槛 $H_i$ 由[单刀亏量三次式](affine-vertex-volumes-proof.md)得到 $H_i/s=u_i/(e_i/H_i)$。对固定 $H_i$，将三个行分数同时放大只会使中央条带体积减小，因为法向分量的绝对值均为两行分数之和的一半；对单个分数也逐坐标单调。门槛为零时以连续性取导数。若所有行分数不超过一，则由该单刀公式、$|H_i|\le\min_j C_{ij}$ 和单调性得 $e_i/H_i\ge(3-1/3)/8=1/3$，故核验可取 $|H_i|\le3(7q-1/2+T)$。

在 (2) 中，$t/s$ 直接用正核权重的中央四面体三次式计算。$\Delta_{ij}/s$ 用[两刀精确平移公式](affine-vertex-volumes-proof.md#两个交替角星平移后的精确分歧)与中点上、下界取交。$O_r/s$ 用[三刀体积差](affine-vertex-volumes-proof.md#三刀同号体积的奇部与七对符号块)的线性项和带因子 $s^2$ 的三次项计算。这些表达式在径向盒子包含 $s=0$ 时仍有界。另用整数向量 $y$，从四个单刀符号等式求得
$$
L=y\cdot(-2u+(t/s)\mathbf1)-(y\cdot\sigma_1)O_1/s
   =\sum_{r\ne1}(y\cdot\sigma_r)O_r/s .
$$
对 $\mu\ge0$，$\sum_{r\ne1}w_r|O_r|/s$ 至少为 $\mu|L|+\sum_{r\ne1}(w_r-\mu|y\cdot\sigma_r|)|O_r|/s$。程序在若干有理 $\mu$ 处取其中最大的下界，随后检验 (2)。转移核验覆盖 **5,115 个方向节点、2,558 个方向叶子及 2,306,972 个径向节点**。零偏移由居中结论覆盖，非零偏移由 (2) 的严格性排除等号。[偏移转移证书](affine_angular_offset_transfer.cpp) · [奇部和两刀界的区间运算](affine_alternating_moments.hpp)

~~~sh
c++ -O3 -std=c++20 affine_angular_merit.cpp -o /tmp/cube-angular-merit
python3 -B affine_local_merit_certificate.py --angular /tmp/cube-angular-merit 35
c++ -O3 -std=c++20 affine_angular_offset_transfer.cpp -o /tmp/cube-angular-transfer
/tmp/cube-angular-transfer 35 20
~~~
