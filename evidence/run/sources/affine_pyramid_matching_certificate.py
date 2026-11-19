"""All-offset exclusions from the cube-center pyramid bounds."""

from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import isqrt

from affine_pyramid_row_certificate import CAP, add, controls, multiply
from centered_bernstein_certificate import plus, scale, times


M = F(79759,1000000)
G = 7*M - F(1,2)
FACE_LIMIT = F(1,6) - M
ROUND = 10**10
LABELS = range(4)
MASKS = range(1,15)


def round_up(a):
    return F(-((-a*ROUND)//1),ROUND)


def sqrt_lower(a):
    return F(isqrt((a*ROUND**2)//1),ROUND)


def sqrt_upper(a):
    return sqrt_lower(a) + F(1,ROUND)


@cache
def deficit_vertices(tetrahedron):
    total = 2*G + tetrahedron
    count = total//M
    residual = total - count*M
    result = set()
    for full in combinations(MASKS,count):
        for last in MASKS:
            if last in full:
                continue
            e = tuple(M*sum(bool(j & (1<<i)) for j in full)
                      + residual*bool(last & (1<<i)) - G for i in LABELS)
            result.add(tuple(sorted(e)))
    return result


def max_square_deficits(tetrahedron):
    assert M < 2*G < 2*M and 2*G + tetrahedron < 3*M
    endpoints = {F(0),tetrahedron}
    transition = 2*M - 2*G
    if transition < tetrahedron:
        endpoints.add(transition)
    return max(sum(x*x for x in e) for t in endpoints for e in deficit_vertices(t))


def slope_bounds(a,b):
    lower = (1+b/3)/(1+b)**2
    for _ in range(3):
        h = min(b,CAP/lower)
        lower = (1+b/2-h*h/(6*b))/(1+b)**2
    upper = (1+a/2)/(1+a)**2
    return lower,upper


def max_tetrahedron(d_low,d_high,lower,upper):
    bound = M
    for _ in range(20):
        improved = F(0)
        for e in deficit_vertices(bound):
            h = [x/(lower if x >= 0 else upper) for x in e]
            for count in range(5):
                small = 4-count
                value = max(F(0),d_low*sum(h[:small]) + d_high*sum(h[small:]))
                volume = value**3/(12*d_low**small*d_high**count)
                improved = max(improved,volume)
        improved = round_up(improved)
        if improved >= bound:
            break
        bound = improved
    return bound


def verify_monotone_loss_coefficient():
    one,x,y = ({(0,0):F(1)}, {(1,0):F(1)}, {(0,1):F(1)})
    a,b = plus(one,x),plus(one,y)
    ab,s = times(a,b),plus(a,b)
    numerator = add(multiply(plus(ab,one),s,plus(x,y)),
                    scale(multiply(a,add(ab,scale(one,-1)),
                                   add(s,scale(one,-1))),-2))
    assert min(controls(numerator,(3,3)).values()) >= 0


def verify_middle(a,b):
    lower,upper = slope_bounds(a,b)
    d_low,d_high = (4*v*(1+v)**2 for v in (a,b))
    tetrahedron = max_tetrahedron(d_low,d_high,lower,upper)
    loss_coefficient = 2*(1+a)**3/(a*(2+a)**2)
    row_losses = 4*b/(6*(1+b)) + loss_coefficient*max_square_deficits(tetrahedron)
    assert row_losses < 4*FACE_LIMIT
    return tetrahedron,row_losses


def minimum_tree_weight(a):
    choices = tuple(product(*(tuple(j for j in LABELS if j != i) for i in (1,2,3))))
    trees = []
    for choice in choices:
        parent = (0,) + choice
        if all(parent[parent[parent[i]]] == 0 for i in (1,2,3)):
            trees.append(choice)
    assert len(trees) == 16
    return min(sum(a**sum(j != k for j,k in zip(pivots,tree)) for tree in trees)
               for pivots in choices)


def verify_all_pivots(a):
    lower,upper = slope_bounds(a,F(1))
    tetrahedron = max_tetrahedron(minimum_tree_weight(a),F(16),lower,upper)
    square_offsets = round_up(max_square_deficits(tetrahedron)/lower**2)
    assert square_offsets < 1
    k_low = a/(4*(1+a)**2)
    k_nonpivot = 1/(4*(1+sqrt_lower(a))**2)
    k_pivot = 1/(4*(1+a)**2)
    k_reference = max((k_low+k_nonpivot)/2,(2*k_low+k_pivot)/3)
    error_nonpivot = max(k_reference-k_low,k_nonpivot-k_reference)
    error_pivot = max(k_reference-F(1,16),k_pivot-k_reference)
    error_norm = sqrt_upper(2*error_nonpivot**2+error_pivot**2)
    nonpivot_quadratic = 1/(8*a*(1+a))
    total_quadratic = k_pivot + 2*nonpivot_quadratic
    total_cubic = k_pivot + 2*nonpivot_quadratic/a
    assert 2-1/(2+a) >= 1/a
    assert 3*CAP < a*(3*a-3*CAP)
    assert total_cubic < 4*total_quadratic
    centered = (1-a/(1+a)**2)/3
    for step in range(20):
        even_loss = square_offsets*(total_quadratic-total_cubic*sqrt_lower(square_offsets)/6)
        three = centered + even_loss - 3*FACE_LIMIT
        if three < 0:
            return tetrahedron,square_offsets,step
        one_norm = 2*sqrt_upper(square_offsets)
        contrasts = (three+error_norm*one_norm)/k_reference
        slope_lower = max(lower,F(3,8)-square_offsets/24)
        mean = (2*G+3*tetrahedron + (upper-slope_lower)*one_norm/2)/((upper+slope_lower)/2)
        improved = round_up((mean**2+contrasts**2)/4)
        assert improved < square_offsets, (a,step,square_offsets,improved)
        square_offsets = improved
    raise AssertionError("offset bounds did not exclude the box")


def main():
    assert 48*M > 1 and (48*M-1)**2 > 8
    assert CAP == 8*M - F(1,2) and 3*CAP < F(1,2)
    verify_monotone_loss_coefficient()
    middle = ((F(9,50),F(1,2)),(F(1,4),F(13,20)),(F(3,10),F(37,50)),
              (F(2,5),F(89,100)),(F(9,20),F(93,100)))
    for a,b in middle:
        tetrahedron,loss = verify_middle(a,b)
        print(f"[{a}, {b}]: tetrahedron <= {float(tetrahedron):.10f}, "
              f"total loss <= {float(round_up(loss)):.10f}",flush=True)
    a = F(7,10)
    tetrahedron,squares,steps = verify_all_pivots(a)
    print(f"All pivot patterns, [{a}, 1]: tetrahedron <= {float(tetrahedron):.10f}, "
          f"offset squares <= {float(squares):.10f} after {steps} refinements",flush=True)


if __name__ == "__main__":
    main()
