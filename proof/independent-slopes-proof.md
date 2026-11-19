# A centered family with three independent slopes and an arbitrary fourth normal

Use the candidate's centered vertex chart and the notation of [the face-overlap formulas](vertex-chart-overlap-proof.md). After rescaling the rows, assume the off-diagonal entries of each high row $i\in\{1,2,3\}$ are $C_{i0}=x_i>0$ and $C_{ij}=1$ for $j\notin\{0,i\}$. Row $0$ has arbitrary positive off-diagonal entries. This gives three independent high-row slopes and two additional parameters for the fourth normal. In spatial coordinates the three high normals are proportional to $(a_1,-1,-1),(-1,a_2,-1),(-1,-1,a_3)$, where $a_i=2/(1+x_i)$, after a choice of cube symmetries. The fourth normal is arbitrary within its vertex sign chart.

In this family the least possible largest cell volume is

$$
m=\frac{1+2\sqrt2}{48}.
$$

Equality forces $x_1=x_2=x_3=\sqrt2-1$ and all three entries of row $0$ equal, giving the candidate construction.

## The volume constraints

All cell variables in this proof are multiplied by $12$. Put

$$
L=12m,\qquad u_* = 1-\frac1{\sqrt2},\qquad
u_i=\frac{x_i}{1+x_i},\qquad g(u)=u-u^2,\qquad g_*=g(u_*).
$$

The identities $2-2L=u_*^2$ and $3(1-L)+3g_*=3/4$ will be used below. In high row $i$, the triangle terms are $t_{i0}=u_i^2$ and $t_{ij}=(1-u_i)/2$ for the other two columns. For $\{i,j,k\}=\{1,2,3\}$, the overlap at the edge with endpoints labeled $0,j$ is

$$
o_{0j}=\frac{(u_i-u_k)^2}{3-2\min(u_i,u_k)}.
\tag{1}
$$

Indeed, when $u_i\ge u_k$, the two legs along the edge are $u_i,1-u_k$. Their perpendicular legs are $u_i,1/2$. The opposite order swaps $i,k$.

Write $O_0=\sum_j o_{0j}$, $O_H=\sum_{1\le i<j\le3}o_{ij}$, and $T_0=\sum_j t_{0j}$. Let $s_j$ be singleton $j$ and $p_j$ the pair cell of the split $\{0,j\}\mid\{i,k\}$. The face formulas give

$$
\begin{aligned}
s_0&=\sum_j u_j^2-O_0,\\
s_j+p_j&=2-u_j^2+O_0-2o_{0j},\\
\sum_j p_j-3L
&=\left[\sum_j g(u_j)+2O_0-3g_*\right]
  -\left[T_0-2O_H-\frac34\right].
\end{aligned}
\tag{2}
$$

The last two identities follow by substituting the triangle terms above. The row $0$ term $t_{0j}$ and the two overlaps involving $j$ and another high label cancel in $s_j+p_j$.

Suppose all cells are at most $L$. Write $a\le b\le c$ for the ordered $u_i$. The largest $o_{0j}$ is the term from the extreme gap $c-a$. For the gap with fixed upper endpoint $c$, $(c-a)^2/(3-2a)$ decreases in $a$, and for fixed $a$ it increases in $c$. At the index whose $u_j=a$, we therefore have $O_0-2o_{0j}\ge0$. The second line of (2) gives $a\ge u_*$.

Let $\mu=(a+b+c)/3$ and $V=\sum_j(u_j-\mu)^2$. For three numbers in an interval of length $R$, $V\le2R^2/3$, since $3V$ is the sum of the three squared gaps and that sum attains its maximum at endpoints. If $a\ge1/2$, then $O_0\le3V$ and $V\le1/6$. If $a<1/2$, the two overlaps involving $a$ have denominator at least $2$. The remaining squared gap is at most half the sum of the three squared gaps, so $O_0\le9V/4$. In this case $V\le1/3$, since $1-a\le1-u_*=1/\sqrt2$. By the first line of (2) we have, in the two cases respectively,

$$
L\ge s_0\ge3\mu^2-\frac13,\qquad
L\ge s_0\ge3\mu^2-\frac5{12}.
$$

Both imply $\mu<17/25$. On the other hand every denominator in (1) is at most $3-2u_*=1+\sqrt2$. Therefore $O_0\ge3(\sqrt2-1)V$ and

$$
\sum_j g(u_j)+2O_0\ge3g(\mu)+\frac{37}{25}V.
\tag{3}
$$

Here $\sum_j g(u_j)=3g(\mu)-V$ and $6(\sqrt2-1)-1>37/25$.

## The fourth row and its overlaps

Write the three entries of row $0$ in increasing order as $A\le B\le C$. Define $r_{AB}=(B-A)/(A+B)$, $r_{BC}=(C-B)/(B+C)$, and $r_{AC}=(C-A)/(A+C)$. From the triangle formula,

$$
T_0=1-\frac{2ABC}{(A+B)(A+C)(B+C)},\qquad
T_0-\frac34=
\frac{r_{AC}^2-(1-r_{AC}^2)r_{AB}r_{BC}}4.
\tag{4}
$$

In particular $0\le T_0-3/4\le1/4$. The first inequality also follows from $(A+B)(A+C)(B+C)\ge8ABC$.

At an edge joining two high labels, take their row $0$ entries $A',B'$ with $A'\ge B'$, and let $C'$ be the remaining entry. The high row supplying the other cut has equal values at the edge's two endpoint labels. Its edge leg is $1/2$ and its perpendicular leg is $1-u_i$, where $i$ is the remaining high index. Thus the overlap equals

$$
\frac{(A'-B')^2}{4(A'+B')^2}
\left[\frac{A'+C'}{A'+B'}+\frac1{2(1-u_i)}\right]^{-1}.
\tag{5}
$$

Let $\widehat O(D)$ be the sum of these three overlap expressions with $1/[2(1-u_i)]$ replaced by a common $D$. We claim that

$$
2\widehat O(2)\ge T_0-\frac34.
\tag{6}
$$

Here is an exact coefficient check. Put $t=A$, $X=B-A$, $Y=C-B$, and set $p=A+B,\ q=B+C,\ r=A+C$, $U=2p+q,\ J=2q+r,\ K=2r+q$. Multiply the difference in (6) by $4pqr^2 UJK$. Substituting (4) and (5) gives

$$
2X^2qr^2JK+2Y^2pr^2UK-q^2(X+Y)^2pUJ+4ACXYUJK.
$$

In the following table, a row for $t^k$ lists the coefficients of $t^k X^{7-k-j}Y^j$ in order $j=0,\ldots,7-k$.

| $k$ | Coefficients |
|---:|---|
| 5 | 288, 288, 288 |
| 4 | 864, 1392, 1392, 480 |
| 3 | 952, 2136, 2536, 1488, 288 |
| 2 | 456, 1324, 1964, 1608, 600, 72 |
| 1 | 80, 288, 570, 658, 380, 90, 6 |
| 0 | 0, 0, 22, 65, 67, 27, 3, 0 |

All are nonnegative. This proves (6) and makes it strict for positive $t$ unless $X=Y=0$.

The true $O_H$ is at least $\widehat O(D)$ for $D=1/[2(1-c)]$. If $c\le3/4$, then $D\le2$, and (6) gives $T_0-2O_H\le3/4$. If $c>3/4$, each first denominator term in (5) is at least $1/2$. Therefore $\widehat O(D)\ge\frac{5/2}{D+1/2}\widehat O(2)$. In combination with (4) and (6) this gives

$$
T_0-2O_H-\frac34
\le\frac{D-2}{4(D+1/2)}=\frac{4c-3}{4(2-c)}.
\tag{7}
$$

## Finishing the comparison

When $c\le3/4$, we have $g(\mu)\ge g_*$ because $u_*\le\mu<17/25<1-u_*$. Equations (2), (3), and (6) give $\sum_jp_j\ge3L$. Equality requires $V=0$, $\mu=u_*$, and equality in (6), so $A=B=C$.

When $c>3/4$, use $V\ge\frac32(c-\mu)^2$, which follows by minimizing the squared deviations of the other two values at fixed sum. The first bracket in the last line of (2) is at least

$$
3[g(\mu)-g_*]+\frac{111}{50}(c-\mu)^2.
\tag{8}
$$

As a function of $\mu$ this is concave, so its minimum on $[u_*,17/25]$ is at an endpoint. At $u_*$, it is greater than $1/4$, using $c>3/4$ and $u_*<3/10$. At $17/25$, use $g_*<20711/100000$. The resulting lower bound for (8) is $3147/100000+(111/50)(c-17/25)^2$, which is strictly larger than the right side of (7) for $c\in[3/4,1]$.

For the last rational comparison, multiply the difference by $4(2-c)$. The resulting cubic has power coefficients

$$
\left(\frac{716499}{62500},\ -\frac{4048199}{125000},
\frac{18648}{625},\ -\frac{222}{25}\right).
$$

In the Bernstein basis of degree three on the four successive intervals of length $1/16$ in $[3/4,1]$, the coefficients multiplied by $96000000$ are

| Interval | Coefficients |
|---|---|
| $[12/16,13/16]$ | 20327040, 15096256, 11097572, 8122863 |
| $[13/16,14/16]$ | 8122863, 5148154, 3197420, 2062536 |
| $[14/16,15/16]$ | 2062536, 927652, 608618, 897309 |
| $[15/16,16/16]$ | 897309, 1186000, 2082416, 3378432 |

All are positive. Thus (8) strictly exceeds the second bracket in (2), giving $\sum_jp_j>3L$. This contradicts the assumed caps. The case $c\le3/4$ gave the unique equality parameters, at which the centered candidate attains $m$.
