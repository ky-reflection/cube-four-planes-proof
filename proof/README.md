# 四个平面切割单位立方体的精确最优值

**定理。** 用四个任意仿射平面切割单位立方体，记最大块的体积为 $M$，记所有切法的下确界为 $C$。允许每个平面独立选择方向和偏移，则

$$
\boxed{C=\min M=m=\frac{1+2\sqrt2}{48}=0.0797588984322122937\ldots}
$$

以立方体中心为原点，等号由以下四个平面取得。除立方体对称和四个平面的重排外，等号构型唯一。

$$
\begin{aligned}
x+y+z&=0, & \sqrt2\,x-y-z&=0,\\
-x+\sqrt2\,y-z&=0, & -x-y+\sqrt2\,z&=0 .
\end{aligned}
$$

它们切出 14 块，12 块的体积各为 $m$，另两块的体积各为 $(3-2\sqrt2)/8$。下界由几何分类和覆盖连续参数域的计算机辅助证明给出。叶子上的严格比较使用有理数或向外舍入的整数区间。[构型体积](edge-separated-family-proof.md#the-centered-minimum) · [全局下界与等号证明](#global-equality) · [最后一类的完整覆盖](reciprocal-pivot-proof.md#affine-1200)

图解和一页截面计算见[知乎版证明](../publication/zhihu-proof.md)。从源码验证主定理见[复现入口](#复现入口)。

<a id="global-equality"></a>

## 全局下界和等号的证明

**第一步。** 先取一般位置的四个平面，令 $q=79759/10^6$。若 $M\le q$，任一刀的每一侧体积至少为 $1-8q>1/6$。只有一个顶点的一侧至多是体积 $1/6$ 的角单纯形。翻转坐标使法向分量非负，则截下的顶点集合对坐标逐项向下封闭。选择顶点数较少的一侧，所得四种可能为一条棱的两个端点、同一面上的三个顶点、一整面、一个顶点及其三个邻点。最后一种称为角星。三刀归约保证每组三刀的交点在立方体内部。对顶点分法按立方体对称和刀号重排枚举，四刀共有 3,746 种组合类型。应用同块顶点、棱切刀、公共凸包及重复角星的体积限制后，余下两棱 206 类、一棱 774 类和无棱 757 类。逐类的全方向、全偏移证明给出前两组及无棱中 756 类均有 $M>q$。[三刀归约](affine-gram-bm-proof.md#passing-to-arbitrary-cuts) · [类型计数与几何排除](affine-edge-and-vertex-proof.md#仿射四刀的分类计数) · [两棱的 206 类](affine-two-edge-proof.md) · [一棱与非交替无棱类](affine-one-and-zero-edge-proof.md#all-nonalternating)

**第二步。** 唯一留下的顶点类型是四个不同的角星。定向后，每行有三个正的顶点评分，把其中一个最大项归一到 1。最大项并列时可任取一个。每一行记下所选列，得到一个没有不动点的四标签映射。$3^4=81$ 种选择在同时重排标签下分为 $1000,1001,1002,1032,1200,1230$ 六个轨道，大小依次为 $12,12,24,3,24,6$。正四面体的每个顶点置换均由立方体对称实现。$1032$ 与 $1230$ 上有 $M>q$。在其余四个排列上，覆盖程序将落入局部定理闭邻域的方向交给该定理，并在其余允许的方向与偏移上排除 $M\le q_*$。这里取
$$
q_*=\frac{\lceil 2^{40}m\rceil}{2^{40}}
=\frac{87695836245}{1099511627776},\qquad m<q_*<q.
$$
被区间收紧删去的部分也满足相应的排除条件。局部定理在八个方向坐标的半径 $1/35$ 内允许四个偏移任取，并给出 $M\ge m$，等号只在所列构型取得。故一般位置的切法都有 $M\ge m$，其中满足 $M\le q_*$ 的方向必在这些闭邻域的对称像之并内。[六种图的全域证明](reciprocal-pivot-proof.md#两对互选的完整仿射证明) · [局部定理及等号条件](local-candidate-proof.md#angular-radius-35)

**第三步。** 对所有平面取单位法向。平面参数收敛时，在极限四个平面以外，每个符号块的指示函数逐点收敛。由于平面的三维体积为零，16 个符号块（空块计零）的体积都连续，故其最大值 $M$ 连续。一般位置的切法稠密，所以 $M\ge m$ 也适用于平面经过顶点、平面平行、四刀共点等退化情形。

若一个构型满足 $M=m$，用一般位置构型趋近它。由于 $m<q_*$，序列从某一项起满足 $M<q_*$。第二步迫使这些项的方向落入局部定理的闭方向域及其有限个对称像之并。每组方向是紧坐标盒在非零法向归一化下的连续像，也就紧且闭。坐标置换、刀号和符号只有有限种选择，故其并集仍闭。因此极限方向仍在其中。对极限构型的实际偏移应用同一局部定理，等号就强制四刀过中心，且各方向恰为上述构型的一个对称像。反过来，上述 14 块的精确体积证明等号取得。这同时给出原下确界的值及所有取到它的切法。[闭方向域上的全部偏移](local-candidate-proof.md#angular-radius-35) · [等号构型的精确体积](edge-separated-family-proof.md#the-centered-minimum)

## 构型与下界估计

**构型的核验。** 上述各块的精确体积来自立方体面上的三角形和以中心为顶点的棱锥分解。独立的局部证书还在 $\mathbb Q(\sqrt2)$ 中核对了这 16 个符号区域的体积，以及候选点处的导数和二阶比较。局部结论允许八个方向参数和四个偏移同时变化，十二个方向及偏移坐标同时变化时，证书覆盖半径为 $1/100$。单独按八个方向坐标计，半径已扩大到 $1/35$，四个偏移可任取。[方向邻域与全部偏移](local-candidate-proof.md#angular-radius-35) · [构型和体积公式](edge-separated-family-proof.md#the-centered-minimum) · [局部证明](local-candidate-proof.md) · [局部证书](local_certificate.py)

若假设所有块都不超过 $0.079759$，中央四面体的体积至多为 $0.008452$。各刀的两侧都至少占 $0.433235$，四把刀到中心的距离都至多为 $0.066765$。[四面体棱长与候选阈值的核验](affine-frame-height-proof.md#用中央四面体的棱长收紧候选阈值)

**下界证明使用的几何归约。** 前两步也用于精确最优值的分类证明。第三步给出较早的独立下界。

1. **三刀归约。** 如果某三刀切出至多六块，第四刀之后仍有一块体积至少 $1/12$。如果这三刀切出七块，整数证书给出三刀最大块严格大于 $4/25$，加第四刀后仍严格大于 $2/25$。所以在候选阈值下，每组三刀都实现八个符号区域。由凸性可知任意三刀独立，且交点位于立方体内部。四刀共点时得到 14 块，非共点时得到 15 块。[归约证明](affine-gram-bm-proof.md#passing-to-arbitrary-cuts) · [七块情形证书](seven_scaled_charts_certificate.cpp)
2. **把体积假设转成法向约束。** 以中央四面体为基准定向。用各符号块与阈值之间的体积差控制每刀的半空间体积、成对符号不同的体积，以及单位法向的 Gram 矩阵。成对夹角与立方体中对应体积之间的比较由角度导数界和若干定角证书支撑。[体积与夹角比较](affine-gram-bm-proof.md#pair-disagreement-and-cap-deficits) · [中央四面体及 Gram 界](affine-frame-height-proof.md) · [角度导数证书](angular_moment_certificate.cpp) · [定角证书](wedge_certificate.cpp)
3. **同时使用两刀和三刀的体积。** 将两刀的 $++$ 区域与反射后的 $--$ 区域取中点，得到由两个半空间体积偏差之和的平方控制的误差界。再用立方体八个角超出给定半径的体积，约束每组三刀的同号区域，把立体角转成四个法向内积和的不等式。六对刀和四组三刀共用同一组切块体积，因此可在同一批 28 个端点上检验。用这组限制与“两正两负”总体积比较，127 个权重子域节点覆盖了所有可能，推出 $M>383/5000$。紧性保证 $C$ 也严格大于该值。[新下界证明](affine-triple-radial-proof.md) · [有理数证书](affine_triple_radial_certificate.py) · [四线径向修正](four-line-radial-proof.md)

## 按切法分类的结果

这里的**交替顶点符号**指可以标号并定向，使同一奇偶类的四个顶点 $v_i$ 分别只有第 $i$ 刀为正，对顶点的符号相反。用边长为二的中心立方体计算时，体积除以八。令 $C_{ij}=-n_i\cdot v_j$（$j\ne i$），严格交替符号等价于每个 $C_{ij}>|H_i|$，其中 $H_i$ 是平面偏移。每行把一个最大分数归一到 $1$ 后，按最大项所在列分类，在对称性下有六种排列。这个符号类别内，三刀交点始终在立方体内部，非共点必有 15 块。由六个面的角三角形、重叠面积和棱锥递推，可精确求出各块体积。[一般体积公式](affine-vertex-volumes-proof.md) · [三刀交点的几何证明](column-weight-affine-proof.md#a-useful-fact-for-the-entire-alternating-affine-chart)

| 范围 | 已证结论及证据 |
| --- | --- |
| 至少三刀过立方体中心，所有方向、第四刀任意 | 最小值为 $m$，等号唯一到立方体对称与刀的重排。四刀居中时，无平面经过顶点的非交替类型均有 $M>2/25$，经过顶点的情形有 $M\ge2/25$。[居中分类定理和第四刀推论](centered-classification-proof.md) |
| 一把整面刀、三个不同角星，四把刀各自任意平移 | 严格顶点类型有 $M>0.079759$，闭包上有 $M\ge0.079759>m$。216 张有向图的 57 个轨道已全部证明。[全偏移定理与逐图记录](affine-one-and-zero-edge-proof.md#face-three-stars-all-offsets) |
| 恰有一把棱切刀，四个偏移均任意 | 余下 774 类全部有 $M>0.079759$，各类闭包上有 $M\ge0.079759>m$。[一棱定理与逐类核验](affine-one-and-zero-edge-proof.md) |
| 恰有两把棱切刀，四个偏移均任意 | 余下 206 类的严格顶点分法均有 $M>0.079759$，各类闭包上有 $M\ge0.079759>m$。[证明和核验记录](affine-two-edge-proof.md) |
| 至少三把棱切刀，或两把给出同一角星的刀，四个偏移均任意 | 整个类型都有 $M>2/25$。另有两把沿正方形对角线方向切下棱的切法，若第三刀垂直于这两棱，在任意三刀偏移及任意第四刀下均有 $M\ge1/12$。具体法向见[定理与证明](affine-edge-and-vertex-proof.md)。 |
| 候选点周围的八维方向邻域，四刀任意平移 | 局部坐标中八个方向扰动均至多 $1/35$ 时，全部独立偏移都有 $M\ge m$，等号只在候选构型取得。[局部定理与全偏移覆盖](local-candidate-proof.md#angular-radius-35) |
| 三个方向参数，四个独立偏移 | 每行适当缩放后 $C_{ij}=c_j$，等价于四个法向分别沿 $v_i-b$，其中 $b$ 在四个奇偶顶点围成的四面体内部。全部允许偏移都已覆盖，最小值为 $m$，非共点时严格大于 $m$。[证明和证书](column-weight-affine-proof.md) |
| 坐标缩放的四面体法向，四个独立偏移 | 这一方向族及一个明确的方向邻域全部排除。另有一个不依赖该对称结构的八方向参数盒，所有允许偏移也已覆盖。[两个全偏移结果](affine-vertex-volumes-proof.md#all-offsets-for-coordinate-rescaled-tetrahedral-directions) |
| 六种最大项排列的共同方向区域 | 八个非最大分数均在 $[0.70,1]$ 时，四个偏移可独立取满，且 $M>0.079759>m$。[六个中心棱锥的证明](affine-pyramid-offset-proof.md#all-pivot-patterns-near-the-regular-directions) |
| 两对互选最大项，排列 $1032$ | 整个八维方向域及全部独立偏移均有 $M>0.079759>m$。683 个参数盒节点覆盖全部方向。[完整仿射证明与复算](reciprocal-pivot-proof.md#两对互选的完整仿射证明) |
| 四个最大项首尾相接，排列 $1230$ | 严格交替类型的全部方向和独立偏移都有 $M>0.079759>m$。方向覆盖有 13,625,249 个节点，最大深度 44。[四循环的全偏移证明](reciprocal-pivot-proof.md#affine-four-cycle) |
| 三行最大项同列，排列 $1000$ | 全部方向及四个独立偏移都有 $M\ge m$，仅候选构型及其等价表示取得等号。1,932,811 个方向节点和 6,394,866 个偏移节点覆盖全域。[插值、偏移单形与等号证明](reciprocal-pivot-proof.md#affine-candidate-charts) |
| 最大列排列 $1001$ | 全部方向及四个独立偏移均有 $M\ge m$，等号只在候选构型的对称像取得。3,195,195 个方向节点全部闭合。[亏量与方向联合比较](reciprocal-pivot-proof.md#affine-candidate-charts) |
| 最大列排列 $1002$ | 全部方向和独立偏移均有 $M\ge m$，等号仅在候选构型的等价表示处取得。直接细分与截面列和两种证明分别检查 6,613,671 和 2,033,224 个方向节点。[两种全域证明](reciprocal-pivot-proof.md#affine-candidate-charts) |
| 最大列排列 $1200$ | 全部方向及四个独立偏移均有 $M\ge m$，等号仅在上述构型的对称像取得。64 个分区共检查 12,370,888 个方向节点，未决数为零。[分区证明](reciprocal-pivot-proof.md#affine-1200) · [逐区记录](affine-alternating-sector-records.txt) |
| 存在一对刀满足互选条件 | 若这对中有一刀的两个非最大分数都不超过 $0.23$，就有 $M>0.079759$，其余方向分数仍可自由变化。另有数个两行分数上界的组合。[两刀证据](reciprocal-pivot-proof.md#checking-the-two-cut-bounds-exactly) |

还有三条适用于一般顶点符号的限制。在最大块体积至多 $2/25$ 时，每块的闭包至多含三个立方体顶点，且不能同时含一对对顶点。四个两顶点半空间也不能恰好把八个顶点配成完美匹配。在 $M\le m$ 的条件下，若所有“两正两负”的棱上区间都不超过 $10^{-35}$，则只能是居中的等号构型。[对顶点排除](antipodal-vertex-proof.md) · [顶点匹配排除](vertex-matching-proof.md) · [棱上区间的定量限制](edge-stability-proof.md)

## 分类细节

### 居中切法的分类定理

四刀过中心且均未经过顶点时，每刀严格划分八个顶点为四对四。这样的无向分法只有七种，三种分开一对相对的面，四种将一个顶点及其三个邻点和对面的四个顶点分开，后者简称角星。允许不同平面有相同的顶点分法。把四刀的无序多重集按立方体的 48 个对称变换合并，得到 22 种；应用已有的同块顶点限制后剩 20 种。[顶点限制](antipodal-vertex-proof.md) · [枚举程序](vertex_partition_census.py)

| 20 种居中顶点分法 | 种数 | 在候选阈值处的处理 |
| --- | ---: | --- |
| 四种角星各出现一次 | 1 | 即严格交替顶点类别，已有 $M\ge m$ 的证明。 |
| 同一角星至少出现两次 | 8 | 下述两刀检验证明 $M>2/25$，其余两刀可任取。 |
| 一把面切刀、三种不同角星 | 1 | 已证 $M>2/25>m$，覆盖整个闭方向域。[证明与复现](centered-face-three-stars-proof.md) |
| 两把面切刀、两个不同角星 | 4 | ZZAD、ZZAB、ZYAD、ZXAD 均有 $M>2/25$。 |
| 三把面切刀、一个角星 | 3 | ZZZA、ZZYA、ZYXA 均有 $M>2/25$。 |
| 四把面切刀，面分法至少两种 | 3 | ZZZY、ZZYY、ZZYX 均有 $M>2/25$。 |

八种排除共用同一个证明。若两刀给出同一角星，将两刀的正侧都朝向同一个角。选顶点 $(1,1,1)$，将法向的系数和归一到 2，则两法向都在闭三角形 $n_1+n_2+n_3=2,\ 0\le n_i\le1$ 中。若两刀分歧体积为 $D$，其同号部分被其余两刀分成至多八块，故 $M\ge(1-D)/8$。在该等边三角形上每边 32 等分，得到 561 个格点。坐标置换可将第一法向排序，余下 57,120 对不同格点逐一以有理数裁切，最大分歧恰为 $1/3$。任一点到最近格点的距离至多 $\sqrt{2/3}/32$，整个三角形上的法向长度至少为 $2/\sqrt3$。沿连接线段归一化的角路程至多 $1/(32\sqrt2)$。由每弧度 $0.51$ 的角度导数界，两法向的总误差小于 $23/1000$。因此连续域上 $D<1069/3000<9/25=1-8(2/25)$，从而 $M>2/25$。这也适用于只有这两刀过中心、另外两刀带偏移的情形。[两刀证书](centered_star_pair_certificate.py) · [角度导数证书](angular_moment_certificate.cpp)

表中 $X,Y,Z$ 表示面切刀分开的相对面，$A,B,C,D$ 表示四种无向角星。每次出现一个字母都允许独立选方向。最后 10 类的 88 个根域共分为 261,839,653 个叶子，每个叶子都排除了 $M\le2/25$。四把面刀均同轴时由体积 $1/6$ 的棱锥直接排除，四个角星均相同时由同角星两刀估计排除，补齐最初的 22 类。含顶点的平面可从面切刀一侧逼近，极限仍有 $M\ge2/25>m$。因此只有交替类能取得最小值。对顶块的体积比较再将定理推广到第四刀任意平移。[十类的逐项记录、边界与等号证明](centered-classification-proof.md)

### 允许四个独立偏移的顶点分类

若切法处于一般位置且 $M\le q=79759/10^6$，则每刀两侧各至少含两个顶点，原因是每侧体积至少 $1-8q>1/6$，而单顶点侧至多为一个体积 $1/6$ 的角单纯形。翻转坐标使各法向系数非负，含顶点较少的一侧在立方体的坐标逐项顺序下向下封闭。含两个顶点时是一条棱，含三个时是一个面上的三个顶点；含四个时是整面或角星。故单刀有 $12+24+3+4=43$ 种无向分法。一般四刀有 3,746 种无序顶点分法轨道，应用同块顶点限制和已证的双顶点侧限制后仍有 3,177 种。这些是组合类型数量，每类内部有连续的方向和偏移。[双顶点侧限制](vertex-matching-proof.md) · [枚举程序](vertex_partition_census.py)

初始的 3,746 类还可以用 Burnside 引理独立计算。对立方体的某个对称，若它在 43 种单刀分法上的循环长度为 $d_1,d_2,\ldots$，其固定的四刀无序多重集数是 $\prod_j(1-t^{d_j})^{-1}$ 的 $t^4$ 系数。48 个对称按此系数的取值分组，固定数与相应的对称个数分别为 $(3,8),(13,12),(15,8),(345,3),(545,6),(885,4),(1413,6),(163185,1)$。加权和为 $179808$，除以 48 即得 $3746$。

三把或四把棱切刀合计的 114 种轨道现已全部排除。对其余轨道，每组选取两刀或三刀的同侧顶点，计算它们的凸包是否已超过剩余切刀所能容纳的体积，再排除出现两把同角星刀的轨道，共剩 1,737 种。两棱的 206 种及一棱的 774 种已全部排除。证明使用一个必为空的完整符号块，将相应三刀母块允许的体积从 $2q$ 收紧为 $q$，其中 $q=0.079759$。无棱的 757 类中，756 类由非交替覆盖结束；交替四角星的六个最大列排列见上表。[两棱分类证明](affine-two-edge-proof.md) · [一棱和无棱的分类证明](affine-one-and-zero-edge-proof.md)

采用“三刀加一刀”。前三刀至多有 377 种顶点分法轨道。固定前三刀和第四刀的方向 $n$，对前三刀的每个块 $P$ 定义 $F_P(h)=\operatorname{vol}(P\cap\{n\cdot x\le h\})$。严格反例的充要条件是存在同一个 $h$，对八个块均有
$$
\operatorname{vol}(P)-m < F_P(h)<m.
$$
各 $F_P$ 连续单调，各块的可行偏移形成开区间，也可为空或整条实线。有一块体积至少 $2m$ 就可立即排除。否则可以用两个块对偏移提出的冲突要求排除一整个方向区域；对一个前三刀类型出具证书时，要覆盖其全部参数和第四刀的全部方向。现有第四刀系数比较可以交叉检查这些排除。[第四刀仿射证书](affine_fourth_certificate.py)

## 复现入口

从本目录运行以下命令，可从源码重新验证主定理所依赖的计算步骤。复算使用 macOS 或 Linux 的 POSIX 文件锁。需要 Python 3.9 或更高版本，以及支持 C++20、128 位整数和未定义行为检测的 Clang 或 GCC。输出目录由程序新建。并行数可按本机资源调整。

~~~sh
python3 -B verify_global_proof.py --jobs 4 ../cube-global-replay
python3 -B check_classification_records.py --replay ../cube-global-replay
~~~

中断后可用 `python3 -B verify_global_proof.py --resume --jobs 4 ../cube-global-replay` 续跑。每个沿用的成功项都会重新核对日志与输入。

程序保存源码快照、SHA-256、编译器版本、逐项命令和原始输出。它先从组合枚举核对覆盖域，再重算三刀归约、角度和体积界、所有保留的非交替代表、六个交替轨道及等号邻域。每个代表和分区均须由对应程序正常退出。逐项结果写入 `results.jsonl`，最后的 `summary.json` 给出总数和失败清单。主定理只需每个域有一条完整证明，例如 $1002$ 使用截面列和证明，$1200$ 使用全部 64 个主分区。其余交叉检查和较早下界保留各自的复现命令。

以下命令在本目录执行。核心数值比较使用有理数或向外舍入的整数区间。部分分盒程序用浮点线性规划寻找权重，排除参数盒之前会用整数重新核验。

~~~sh
# 分类清单与逐区记录的一致性核对
python3 -B check_classification_records.py

# 顶点分类及同角星两刀排除
python3 -B vertex_partition_census.py
python3 -B centered_star_pair_certificate.py

# 两棱的 206 个代表与独立偏移的逐类覆盖
python3 -B vertex_partition_census.py --two-edges > /tmp/cube-two-edge-types.txt
c++ -std=c++20 -O3 affine_vertex_type_cover.cpp -o /tmp/cube-two-edges
/tmp/cube-two-edges < /tmp/cube-two-edge-types.txt

# 独立偏移下的三棱定理和重复角星定理
c++ -std=c++20 -O3 affine_edge_triples.cpp -o /tmp/cube-affine-edge-triples
/tmp/cube-affine-edge-triples
c++ -std=c++20 -O3 affine_star_pair.cpp -o /tmp/cube-affine-star-pair
/tmp/cube-affine-star-pair

# 一把面切刀配三个不同角星，以及最后十个多面切刀类
c++ -std=c++20 -O3 centered_face_three_stars.cpp -o /tmp/cube-face-stars
for word in XBCD ZZAD ZZAB ZYAD ZXAD ZZZA ZZYA ZYXA ZZZY ZZYY ZZYX
do
    /tmp/cube-face-stars "$word" || exit 1
done

# 候选构型、一般局部证书与最近两组仿射排除
python3 -B local_certificate.py
python3 -B affine_pyramid_row_certificate.py
python3 -B affine_pyramid_matching_certificate.py
python3 -B reciprocal_pivot_certificate.py

# 全局下界的三刀归约和新的联合体积证书
c++ -std=c++20 -O3 seven_scaled_charts_certificate.cpp -o /tmp/cube-seven
/tmp/cube-seven --horizontal 30
/tmp/cube-seven --vertical 36
python3 -B affine_triple_radial_certificate.py

# 七参数方向及全部独立偏移
python3 -B column_weight_local_certificate.py
c++ -std=c++20 -O3 column_weight_affine_cover.cpp -o /tmp/cube-column
/tmp/cube-column 98 5000000 24
~~~

最后一个排列的逐区复算命令在[排列 $1200$ 的证明](reciprocal-pivot-proof.md#affine-1200)。较早的全局下界还依赖[角度和定角比较](affine-gram-bm-proof.md#reproducing-the-arithmetic)、[扩展定角列表](affine-frame-height-proof.md#reproduction)和[成对分歧插值的定角列表](affine-weight-graph-proof.md#reproduction)。三刀证书须运行上面两个模式，无参数运行只检查较小子域。居中六种排列的命令见[第一组](centered-bernstein-proof.md#reproduction)和[第二组](centered-bernstein-completion.md#reproduction)。有些居中程序会把候选点附近的盒子交给单独的局部证书，须将分盒结果和对应局部证书一起检查。

### 本次从源码全量重算

2,005 个任务均已通过。1,736 个非交替顶点类、交替类的六种方向排列、三刀归约及局部等号证明均已覆盖，所有必需分区未决数为零。每项日志的输出摘要和退出状态均已核对。可查阅[本次验证报告](../evidence/README.md)和[逐项运行账本](../evidence/run/results.jsonl)。

### 既有证书的运行记录

| 证据 | 运行结果 |
| --- | --- |
| 全局最小值及等号 | 按组合枚举重新核对两棱、一棱和无棱的逐类记录，以及 81 种最大列选择的六个轨道。$1200$ 的 64 个分区均已闭合。另以精确体积复算等号构型的 16 个符号块，方向半径 $1/35$ 的 256 个二阶型和 2,306,972 个偏移径向节点复核通过。 |
| 较早的全局下界 $0.0766$ | 三刀归约的两组覆盖分别检查 70,079 和 218,663 个节点，均无未决盒。新证书核验 28 个体积端点，权重单纯形覆盖共 127 个节点、64 个叶子，最大深度 14，最小矛盾余量大于 $0.000018569$。两项角度导数检验和 45 次定角检验均已复跑通过。 |
| 此前的全局界 | $1/14+11/2500=0.075828571428\ldots$ 的证明使用 68,103 个权重盒节点，矛盾余量约为 $0.00032930$。[原证明](affine-weight-graph-proof.md) · [原证书](affine_weight_graph_certificate.py) |
| 候选点及居中交替符号 | 一般局部证书通过。六种居中排列均复现了已记录的结果，交给局部证明的盒子均在其适用半径内。 |
| 至少两把面切刀的十个居中类 | 逐类核对 349,119,508 个节点、261,839,653 个叶子，最大深度 23，全部排除 $M\le2/25$。[分项结果](centered-classification-proof.md) |
| 独立偏移的分盒证明 | 七参数覆盖检查 852,217 个节点，其中 16 个盒子由局部证书处理，未决数为零。坐标缩放方向的两次覆盖，以及非对称方向盒的覆盖，均无未决盒。 |
| 中心棱锥与互选两刀 | 逐行界、合并界和互选证书均通过，包括 $[0.70,1]$ 的共同方向区域，以及 $1032$ 下非最大分数均不超过 $0.694$ 的区域。 |

## 文件索引

| 用途 | 相关文件 |
| --- | --- |
| 全局下界的中间界及替代径向估计 | [较早的 Gram 界](affine-gram-bound-proof.md) · [对应证书](affine_gram_bound.py) · [仿射成对比较证书](affine_gram_bm_certificate.py) · [中央四面体与 Gram 证书](affine_frame_height_certificate.py) · [四线径向证书](four_line_radial_certificate.py) · [耦合径向笔记](coupled-radial-proof.md) · [对应证书](coupled_radial_certificate.py) · [早期全局下界常数](global_lower_bound_constants.py) |
| 三刀归约的基础工具及早期证书 | [整数体积工具](seven_integer_volume.hpp) · [三刀方向盒](seven_box_certificate.py) · [三刀 Gram 估计](seven_gram_certificate.py) · [偶数符号块比较](even_seven_certificate.py) · [对称七块情形](symmetric_seven_certificate.py) |
| 居中交替符号的公式与分盒 | [面面积和重叠公式](vertex-chart-overlap-proof.md) · [整数分盒](centered_chart_cover.cpp) · [Bernstein 核验](centered_bernstein_certificate.py) · [居中候选邻域](centered_candidate_neighborhood.py) |
| 独立偏移的体积计算与分盒 | [精确体积计算](affine_vertex_volumes.py) · [共享几何公式](affine_box_geometry.hpp) · [仿射区间运算](affine_box_linear.hpp) · [坐标缩放方向](rescaled_tetrahedral_affine_cover.cpp) · [非对称方向盒](affine_direction_box_cover.cpp) · [列权方向分盒](column_weight_affine_cover.cpp) · [列权候选邻域](column_weight_local_certificate.py) |
| 中心棱锥、互选两刀与顶点约束 | [逐行棱锥界](affine_pyramid_row_certificate.py) · [合并棱锥界](affine_pyramid_matching_certificate.py) · [互选两刀证书](reciprocal_pivot_certificate.py) · [对顶点证书](antipodal_vertex_certificate.py) · [棱上分离族的仿射证书](edge_separated_affine_certificate.py) |
| 其他特定方向族和第四刀估计 | [固定体对角法向](fixed-tetrahedral-normals-proof.md) · [前三刀循环对称](cyclic-family-proof.md) · [三个独立斜率](independent-slopes-proof.md) · [候选法向附近的偏移约束](candidate_translates_certificate.py) · [第四刀估计](fourth_plane_certificate.py) · [第四刀仿射估计](affine_fourth_certificate.py) |
