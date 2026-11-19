"""Exact affine cell volumes in the alternating vertex chart.

Rows of C give the three positive numbers -n_i dot v_j, j != i.
The input offsets are in the side-two cube [-1,1]^3. Volumes are
normalized by eight. See affine-vertex-volumes-proof.md.
"""

from fractions import Fraction as F
from itertools import combinations


VERTICES = ((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))
LABELS = range(4)
ALL = 15


def dot(a, b):
    return sum(x*y for x,y in zip(a,b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def normals_from_scores(c):
    normals = []
    for i in LABELS:
        values = tuple(-c[i][j] if j != i else sum(c[i]) for j in LABELS)
        normals.append(tuple(dot(values,column)/4 for column in zip(*VERTICES)))
    return tuple(normals)


def triple_points(normals, offsets):
    points = []
    for omitted in LABELS:
        i,j,k = (a for a in LABELS if a != omitted)
        a,b,c = normals[i], normals[j], normals[k]
        directions = (cross(b,c), cross(c,a), cross(a,b))
        denominator = dot(a,directions[0])
        assert denominator != 0
        points.append(tuple(dot(column,(offsets[i],offsets[j],offsets[k]))/denominator
                            for column in zip(*directions)))
    return tuple(points)


def face_cells(c, offsets, pair):
    """Areas as fractions of one full face, indexed by positive sign set."""
    face_mask = sum(1 << i for i in pair)
    shift = pair[0] ^ pair[1]
    triangles = [F(0)] * 4
    overlap_sum = [F(0)] * 4
    for i in LABELS:
        j = i ^ shift
        k,l = (a for a in LABELS if a not in (i,j))
        h = offsets[i] if i in pair else -offsets[i]
        triangles[i] = (c[i][j]+h)**2/(2*(c[i][j]+c[i][k])*(c[i][j]+c[i][l]))
    areas = [F(0)] * 16
    for i in pair:
        for j in LABELS:
            if j in pair:
                continue
            a,b = i ^ shift, j ^ shift
            u = (c[i][a]+offsets[i])/(c[i][a]+c[i][b])
            v = (c[j][b]-offsets[j])/(c[j][a]+c[j][b])
            d = (c[i][a]+c[i][j])/(c[i][a]+c[i][b]) + (c[j][b]+c[j][i])/(c[j][a]+c[j][b])
            overlap = max(F(0),u+v-1)**2/(2*d)
            areas[face_mask ^ (1<<i) ^ (1<<j)] = overlap
            overlap_sum[i] += overlap
            overlap_sum[j] += overlap
    areas[face_mask] = 1-sum(triangles)+sum(overlap_sum)/2
    for i in LABELS:
        areas[face_mask ^ (1<<i)] = triangles[i]-overlap_sum[i]
    assert min(areas) >= 0 and sum(areas) == 1
    return tuple(areas)


def volumes(c, offsets):
    """Return all 16 exact volumes when all four triple points are interior."""
    c = tuple(tuple(map(F,row)) for row in c)
    offsets = tuple(map(F,offsets))
    assert all(c[i][i] == 0 and abs(offsets[i]) < min(c[i][j] for j in LABELS if i != j) for i in LABELS)
    normals = normals_from_scores(c)
    points = triple_points(normals,offsets)
    assert all(abs(v) < 1 for point in points for v in point)
    exterior = all(dot(normals[i],points[i]) <= offsets[i] for i in LABELS)
    interior = all(dot(normals[i],points[i]) >= offsets[i] for i in LABELS)
    assert exterior or interior
    if not exterior:
        # Center reflection changes offsets and reverses all cell labels.
        return tuple(reversed(volumes(c,tuple(-h for h in offsets))))
    faces = []
    for pair in combinations(LABELS,2):
        axis = (pair[0] ^ pair[1])-1
        sign = 1 if pair[0] == 0 else -1
        faces.append((axis,sign,face_cells(c,offsets,pair)))
    result = [F(0)] * 16
    for cardinality in (3,2,1):
        for mask in range(1,ALL):
            if bin(mask).count('1') != cardinality:
                continue
            j = next(i for i in LABELS if not (mask & (1<<i)))
            larger = mask | (1<<j)
            total = sum((1-sign*points[j][axis])*(areas[mask]+areas[larger])
                        for axis,sign,areas in faces)/6
            result[mask] = total-result[larger]
    result[0] = 1-sum(result)
    assert min(result) >= 0
    return tuple(result)
