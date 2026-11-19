"""Pair disagreements for independent affine cuts in the alternating chart."""
from fractions import Fraction as F
from itertools import combinations
from math import prod

from affine_vertex_volumes import normals_from_scores, triple_points


def absolute_mean(a, b, c):
    """E|a X + b Y + c| for X,Y independent uniform on [-1/2,1/2]."""
    a, b = sorted((abs(a), abs(b)), reverse=True)
    c = abs(c)
    if 2*c >= a+b:
        return c
    if 2*c+b <= a:
        return a/4 + (c*c+b*b/12)/a
    return c + (a+b-2*c)**3/(24*a*b)


def disagreement(c, h, i, j):
    """Return the exact normalized cube volume where cuts i and j disagree."""
    c = tuple(tuple(map(F, row)) for row in c)
    h = tuple(map(F, h))
    normals = normals_from_scores(c)
    k = (i^j)-1
    x, y = (d for d in range(3) if d != k)
    slopes, heights = [], []
    losses = F()
    for r in (i, j):
        n = normals[r]
        a, b = abs(n[x]/n[k]), abs(n[y]/n[k])
        z = h[r]/(2*n[k])
        edge = (a+b-1)/2
        assert abs(z) <= edge and edge+abs(z) <= min(a,b)
        slopes.append((a,b))
        heights.append(z)
        losses += ((edge+z)**3+(edge-z)**3)/(6*a*b)
    ax, ay = map(sum, zip(*slopes))
    return absolute_mean(ax, ay, heights[0]-heights[1])-losses


def offset_loss(c, h, i, j):
    """D(0)-D(h), written as two row quadratics minus an even convex increment."""
    c = tuple(tuple(map(F, row)) for row in c)
    h = tuple(map(F, h))
    normals = normals_from_scores(c)
    k = (i^j)-1
    x, y = (d for d in range(3) if d != k)
    a, b, z = sum(abs(normals[r][x]/normals[r][k]) for r in (i,j)), \
              sum(abs(normals[r][y]/normals[r][k]) for r in (i,j)), \
              h[i]/(2*normals[i][k]) - h[j]/(2*normals[j][k])
    row_sum = F()
    for r, s in ((i, j), (j, i)):
        entries = [c[r][v] for v in range(4) if v != r]
        row_sum += c[r][s]*h[r]**2/prod(p+q for p,q in combinations(entries,2))
    return row_sum - (absolute_mean(a,b,z)-absolute_mean(a,b,F(0)))


def triple_odd(c, h, omitted):
    """Volume of three positive sides at -h minus their volume at h."""
    c = tuple(tuple(map(F, row)) for row in c)
    h = tuple(map(F, h))
    rows = [i for i in range(4) if i != omitted]
    point = triple_points(normals_from_scores(c), h)[omitted]

    def face_height(i, j):
        return point[(i ^ j)-1] * (1 if i == 0 or j == 0 else -1)

    def edge_terms(i, j, k):
        left = 1 / (c[i][k] + c[i][omitted])
        right = 1 / (c[j][k] + c[j][omitted])
        x = c[i][k]*left + c[j][omitted]*right - 1
        y = h[i]*left - h[j]*right
        d = (c[i][k]+c[i][j])*left + (c[j][omitted]+c[j][i])*right
        return x, y, d

    total = F()
    for k in rows:
        i, j = [v for v in rows if v != k]
        g = c[k][omitted]
        denominator = (g+c[k][i])*(g+c[k][j])
        total += (4*g*h[k] + 2*face_height(i,j)*(g*g+h[k]**2)) / denominator
    for i, j in combinations(rows, 2):
        k = next(v for v in rows if v not in (i, j))
        x, y, df = edge_terms(i, j, k)
        _, _, dk = edge_terms(j, i, k)
        total += 2*(1/df+1/dk)*x*y
        total -= (face_height(i,k)/df+face_height(j,k)/dk)*(x*x+y*y)
        total -= (abs(x+y)**3-abs(x-y)**3)/(df*dk)
    return total/12


if __name__ == '__main__':
    from random import Random
    from affine_vertex_volumes import volumes
    rng = Random(532910)
    verified = 0
    for _ in range(360):
        c = [[F(rng.randrange(1,30)) if i != j else F(0) for j in range(4)] for i in range(4)]
        h = [min(s for s in row if s)*rng.randrange(-9,10)/10 for row in c]
        v = volumes(c,h)
        singles = []
        caps = []
        for i in range(4):
            mask = 1 << i
            odd = triple_odd(c,h,i)
            assert odd == v[0]+v[mask]-v[15]-v[15^mask]
            singles.append(odd-v[0]+v[15])
            caps.append(F(1,2)-sum(v[r] for r in range(16) if r & mask))
        for mask in (3,5,9):
            pair_odd = -sum((1 if mask & (1 << i) else -1)*(caps[i]+singles[i])
                            for i in range(4))/2
            assert pair_odd == v[mask]-v[15^mask]
        for i,j in combinations(range(4),2):
            exact = sum(mass for m,mass in enumerate(v) if ((m>>i) ^ (m>>j)) & 1)
            assert disagreement(c,h,i,j) == exact
            assert disagreement(c,(0,)*4,i,j)-exact == offset_loss(c,h,i,j)
            verified += 1
    print(verified, 'pair disagreements and offset differences checked exactly')
    print(360*4, 'triple differences and', 360*3, 'pair-cell differences checked exactly')
