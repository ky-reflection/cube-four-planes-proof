"""An exact deficit-simplex certificate for the three-plane seven-cell case.

Orient the three positive halfspaces so that +++ is missing. Let K be the
hypothetical upper bound for each cell, e_S=K-volume(S), and E=7K-1.
The deficits form a simplex e>=0, sum(e)=E. A positive cap has deficit
from 1/2 equal to z_i=1/2-3K+sum(e_S for i in S). Its unit offset
h_i lies in [0,z_i] by the thin-slab theorem.

Let p be the centered overlap of the i,j caps, and let d=K-e_{ij}
be their affine overlap. The Brunn--Minkowski tangent at centered cuts
gives d <= p*(1-(z_i+z_j)/(6*p))_+^3. Indeed, the centered derivative
with respect to each offset is half the central section area, and the
one-sided cap deficit is at most that section area times the offset.
Inverting this inequality yields p>=P(d,z_i+z_j). P is jointly concave:
its hypograph is the epigraph in d of the convex perspective
p*(1-s/(6*p))_+^3. It also satisfies P<=d+s/2<=3K/2.

At any certified angle anchor, a piecewise affine concave nondecreasing
function of P bounds the dot product of the two unit normals from below.
Thus, fixing one anchor for each pair, the lower quadratic form
lambda^T G lambda - 4*(lambda dot z)^2 is concave in the deficits for
each lambda>=0. Positive definiteness at all vertices of a deficit
simplex therefore certifies positivity throughout it for lambda>=0.
But separation of the missing +++ region gives a nonzero lambda>=0
with ||sum(lambda_i*n_i)|| <= 2*sum(lambda_i*h_i), a contradiction.

For the stronger check, choose a separator whose resultant is on a cube
axis. The function ||sum(lambda_i*n_i)||_1/2 - lambda dot h is affine
on each polygon cut out of the weight simplex by the three coordinate
zero lines. Its minimum is attained at one of the polygon vertices.
It is positive on the boundary of the simplex because every pair of
positive halfspaces meets the cube interior. Hence a minimum occurs at
an interior vertex with at least two zero resultant coordinates.

For a resulting unit-axis component r_a, let r0*e be the point of the
normal triangle nearest the origin. A lower dot bound rho_ij excludes
an edge minimizer when r_a^2 < (1+min rho_ij)/2. In the interior case
all normals have the same e-component r0, so r0 = r_a*abs(e dot axis).
The projections into e-perpendicular have lengths sqrt(1-r0^2).
Their three angular gaps phi_i sum to 2*pi. Positive relation weights
are proportional to sin(phi_i), by the planar cofactor identity.
Consequently the determinant D of two weighted projected normals is
(1-r0^2)*prod(sin(phi_i))/(sum(sin(phi_i)))^2. At axial height s>0
the centered all-positive cone has triangular area r0^2*s^2/(2*D).
Integrate up to the cube support ||e||_1/2, or split off a corner cap.

Write w for the centered all-positive volume and Delta=sum(z_i).
The sum of the three centered pair overlaps is 1/2+2w by symmetry.
Their affine sum is 1/2+t-Delta, where t is the all-negative volume.
The all-negative cell and the three single-positive cells give
t >= .01 + Delta/2 when K<=.16. For every pair, write s=z_i+z_j
and x=s/(6p). Since d>=.04 and s<=.28, the overlap inequality implies
x<.34. Then p-d >= s*(3-3x+x^2)/6 >= .349*s. Thus
w >= .005 + .099*Delta, as well as the bound from the three P roots.

Verify each non-default angular anchor (A,B,U) with:
    wedge_certificate --pair A B U
and verify the angular derivative bounds using angular_moment_certificate
with no arguments and with --upper. The default wedge_certificate run
certifies the anchor whose negative dot product is 1/3.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import factorial

from seven_box_certificate import sqrt_lower, sqrt_upper

PAIRS = ((0, 1), (0, 2), (1, 2))
PATTERNS = ((), (0,), (1,), (2,), *PAIRS)
DERIVATIVE_LOWER = F(23, 100)
DERIVATIVE_UPPER = F(51, 100)
GRID = 1 << 23


def anchor_from_half_angle(cosine, sine, upper_thousandths):
    squared = cosine**2 + sine**2
    return (F(sine**2-cosine**2, squared),
            F(2*cosine*sine, squared), F(1)-F(upper_thousandths, 1000))


ANCHOR_61 = anchor_from_half_angle(100, 61, 384)
ANCHOR_605 = anchor_from_half_angle(200, 121, 381)
ANCHORS = (
    (-F(1, 3), None, F(143, 250)),
    ANCHOR_61,
    ANCHOR_605,
    anchor_from_half_angle(12, 7, 371),
    anchor_from_half_angle(7, 4, 366),
    anchor_from_half_angle(200, 113, 363),
)
WIDE_ANCHORS = ANCHORS + (
    anchor_from_half_angle(5, 4, 469),
    anchor_from_half_angle(9, 7, 459),
    anchor_from_half_angle(5, 2, 287),
    anchor_from_half_angle(2, 1, 335),
)


def dot_lower(overlap, anchor):
    rho, sine, disagreement = anchor
    gap = 1 - 2*overlap - disagreement
    if gap >= 0:
        sine = sine if sine is not None else sqrt_upper(F(8, 9))
        return rho - sine*gap/DERIVATIVE_LOWER
    sine = sine if sine is not None else sqrt_lower(F(8, 9))
    return rho - sine*gap/DERIVATIVE_UPPER


def floor_grid(value):
    return (value.numerator * GRID) // value.denominator


def ceil_grid(value):
    return F(-floor_grid(-value), GRID)


def overlap_lower(overlap, cap_sum):
    """Round down the root of the Brunn--Minkowski tangent bound."""
    lo, hi = 0, GRID//2
    sixth = cap_sum/6
    while hi-lo > 1:
        mid = (lo+hi)//2
        p = F(mid, GRID)
        if max(F(0), p-sixth)**3 < overlap*p*p:
            lo = mid
        else:
            hi = mid
    return F(lo, GRID)


def positive_definite(diagonal, offdiag):
    a, d, f = diagonal
    b, c, e = offdiag
    return a > 0 and a*d > b*b and a*d*f+2*b*c*e > a*e*e+d*c*c+f*b*b


def bisect_simplex(simplex):
    i, j = max(combinations(range(7), 2), key=lambda edge:
               sum((a-b)**2 for a, b in zip(simplex[edge[0]], simplex[edge[1]])))
    midpoint = tuple((a+b)/2 for a, b in zip(simplex[i], simplex[j]))
    for endpoint in (i, j):
        child = list(simplex)
        child[endpoint] = midpoint
        yield tuple(child)


def initial_simplex(target):
    excess = 7*target-1
    return tuple(tuple(excess if i == j else F(0) for j in range(7))
                 for i in range(7))


def cap_deficits(target, deficit):
    return tuple(F(1, 2)-3*target
                 + sum(e for e, s in zip(deficit, PATTERNS) if i in s)
                 for i in range(3))


def certify(target, anchors=ANCHORS, max_nodes=20000):
    """Cover the full deficit simplex, or report an unresolved region."""
    assert F(1, 7) < target < F(1, 6)
    assert all(dot_lower(3*target/2, anchor) < 0 for anchor in anchors)
    # For acute angles the lower bound is trivial from rho>=0 and the
    # preceding check. For obtuse angles it is the tangent to cosine at
    # the anchor, using the appropriate one-sided bound on D'(theta).
    anchor_count = len(anchors)
    @lru_cache(maxsize=50000)
    def allowed(deficit):
        cap = cap_deficits(target, deficit)
        diagonal = tuple(floor_grid(1-4*z*z) for z in cap)
        off = []
        for k, (i, j) in enumerate(PAIRS):
            p = overlap_lower(target-deficit[4+k], cap[i]+cap[j])
            off.append(tuple(floor_grid(dot_lower(p, anchor)-4*cap[i]*cap[j])
                             for anchor in anchors))
        if not positive_definite(diagonal, tuple(map(max, off))):
            return 0
        mask = 0
        for index, values in enumerate(product(*off)):
            if positive_definite(diagonal, values):
                mask |= 1 << index
        return mask

    pending = [(initial_simplex(target), 0)]
    closed = deepest = visited = 0
    while pending:
        simplex, depth = pending.pop()
        visited += 1
        common = (1 << anchor_count**3)-1
        for vertex in simplex:
            mask = allowed(vertex)
            if not mask:
                raise RuntimeError(f"unresolved vertex at depth {depth}: {vertex}")
            common &= mask
        if common:
            closed += 1
            deepest = max(deepest, depth)
            continue
        if visited >= max_nodes:
            raise RuntimeError(f"unresolved deficit simplex after {visited} nodes")
        pending.extend((child, depth+1) for child in bisect_simplex(simplex))
    return closed, deepest


def cosine_bound(overlap, anchor):
    """Tighter constant dot lower bound using the cosine addition formula."""
    rho, sine, disagreement = anchor
    gap = 1-2*overlap-disagreement
    angle = gap/(DERIVATIVE_LOWER if gap >= 0 else DERIVATIVE_UPPER)
    tangent_bound = dot_lower(overlap, anchor)
    if abs(angle) > 1:
        return max(-F(1), tangent_bound)
    if angle >= 0:
        # pi-theta = 2 atan(t), with t=sin(theta)/(1-cos(theta)).
        tangent = (sine if sine is not None else sqrt_lower(F(8, 9)))/(1-rho)
        if angle > 2*(tangent-tangent**3/3):
            return max(-F(1), tangent_bound)
    cosine_upper = sum((-1)**i*angle**(2*i)/factorial(2*i) for i in range(5))
    sine_upper = sum((-1)**i*angle**(2*i+1)/factorial(2*i+1)
                     for i in range(5 if angle >= 0 else 4))
    if sine is None:
        sine = sqrt_upper(F(8, 9)) if angle >= 0 else sqrt_lower(F(8, 9))
    return max(tangent_bound, rho*cosine_upper-sine*sine_upper)


def maximum_weighted_offset(rho, caps):
    """Bound lambda dot caps subject to the Gram constraint on lambda."""
    diagonal = tuple(1-4*q*q for q in caps)
    offdiag = tuple(r-4*caps[i]*caps[j] for r, (i,j) in zip(rho, PAIRS))
    if positive_definite(diagonal, offdiag):
        return None
    a, d, f = diagonal
    b, c, e = offdiag
    # On lambda_1+lambda_2+lambda_3=1 the quadratic has 2x2 Hessian form.
    u, v, w = a-2*c+f, b-c-e+f, d-2*e+f
    determinant = u*w-v*v
    if u <= 0 or determinant <= 0:
        return max(caps)
    g, h = c-f, e-f
    center = ((v*h-w*g)/determinant, (v*g-u*h)/determinant)
    minimum = f+g*center[0]+h*center[1]
    if minimum > 0:
        return None
    x, y = caps[0]-caps[2], caps[1]-caps[2]
    radius_squared = -minimum*(w*x*x-2*v*x*y+u*y*y)/determinant
    return min(max(caps), caps[2]+x*center[0]+y*center[1]+sqrt_upper(radius_squared, GRID))


def projected_area_lower(rho, square):
    """Lower bound on the weighted area determinant in the nearest plane."""
    dots = tuple(min(F(0), (r-square)/(1-square)) for r in rho)
    sine = tuple(sqrt_lower(1-r*r, GRID) for r in dots)
    lower = tuple(min(sine[i], -sine[(i+1)%3]*dots[(i+2)%3]
                      -sine[(i+2)%3]*dots[(i+1)%3]) for i in range(3))
    if min(lower) <= 0:
        return F(0)
    # For each coordinate, xyz/(x+y+z)^2 is unimodal with an interior
    # maximum, so the rectangle's minimum is at one of its eight corners.
    sine_area = min(x*y*z/(x+y+z)**2
                    for x,y,z in product(*((s,F(1)) for s in lower)))
    # Put theta_i=pi-phi_i. The same determinant equals
    # (1-square)*product(tan(theta_i/2))/2. Holding the third angle
    # fixed, the log of the product is concave under an opposite
    # variation of the first two. Its minimum is at a vertex of the
    # simplex theta_i >= pi-arccos(dots[i]), sum(theta_i)=pi.
    squares = tuple((1+c)/(1-c) for c in dots)
    tangent_lo = tuple(sqrt_lower(q, GRID) for q in squares)
    tangent_hi = tuple(sqrt_upper(q, GRID) for q in squares)
    tangent_area = min(tangent_lo[i]*tangent_lo[j]*(1-tangent_hi[i]*tangent_hi[j])
                       / (2*(tangent_hi[i]+tangent_hi[j])) for i,j in PAIRS)
    return (1-square)*max(sine_area, tangent_area)


def cone_upper(rho, offset, refine=False):
    """Upper bound on the centered +++ volume from an axis Farkas witness."""
    square = 4*offset**2
    if square >= (1+min(rho))/2:
        return F(1, 2)
    area = projected_area_lower(rho, square)
    crude = square*F(147,50)/(48*area) if area else F(1, 2)
    if not refine:
        return crude
    bounds = []
    for step in range(16):
        lo, hi = F(step,16), F(step+1,16)
        if hi**2*(1+F(627,1000)**2) <= 1:
            endpoint = hi
        elif lo**2*(1+F(626,1000)**2) >= 1:
            endpoint = lo
        else:
            endpoint = None
        factor = (F(147,50) if endpoint is None else min(
            F(147,50), endpoint**2*(endpoint+sqrt_upper(2*(1-endpoint**2), GRID))**3
        ))
        area = projected_area_lower(rho, hi**2*square)
        if not area:
            bounds.append(F(1, 2))
            continue
        full_cone = square*factor/(48*area)
        if not lo or hi == 1:
            bounds.append(full_cone)
            continue
        # For e=(a,b,c) in absolute coordinates, where a is its axial
        # component, split according to b+c <= q*sqrt(2*(1-a^2)).
        # In that case the support length is smaller. Otherwise b*c
        # is bounded below and the cube cap clips off the cone tip.
        q = F(9,10)
        transverse = sqrt_lower(2*(1-hi**2), GRID)
        narrow = full_cone*((hi+q*transverse)/(hi+transverse))**3
        support = (hi+sqrt_upper(2*(1-lo**2), GRID))/2
        coefficient_product = lo*(1-hi**2)*(q*q-F(1,2))
        tip = (support/8)**3/(6*coefficient_product)
        wide = min(full_cone, full_cone*F(7,8)**3+tip)
        bounds.append(max(narrow, wide))
    return min(crude, max(bounds))


def certify_with_cone(target, anchors=WIDE_ANCHORS, max_nodes=20000):
    """Check the deficit simplex with both the Gram and cone constraints."""
    assert F(1, 7) < target <= F(4, 25)
    # The maximum of (1+sqrt(2)*t)^3/(1+t*t)^(5/2) is below 2.94.
    lo, hi = F(626,1000), F(627,1000)
    s_lo, s_hi = F(1414,1000), F(1414214,1000000)
    assert s_lo*s_lo < 2 < s_hi*s_hi
    assert s_lo*(3-2*lo*lo) > 5*lo and s_hi*(3-2*hi*hi) < 5*hi
    assert (1+s_hi*hi)**6 < F(147,50)**2*(1+lo*lo)**5

    @lru_cache(maxsize=50000)
    def data(deficit):
        caps = cap_deficits(target, deficit)
        overlap = tuple(overlap_lower(target-deficit[4+k], caps[i]+caps[j])
                        for k,(i,j) in enumerate(PAIRS))
        return caps, overlap

    @lru_cache(maxsize=50000)
    def correlation(overlap):
        return F(floor_grid(max(cosine_bound(overlap,a) for a in anchors)), GRID)

    def closed(simplex):
        points = tuple(data(vertex) for vertex in simplex)
        rho = tuple(correlation(min(p[k] for _,p in points)) for k in range(3))
        offsets = tuple(h for caps,_ in points
                        if (h:=maximum_weighted_offset(rho,caps)) is not None)
        if not offsets:
            return True
        w_lower = max(min((sum(p)-F(1,2))/2 for _,p in points),
                      F(1,200)+F(99,1000)*min(sum(z) for z,_ in points))
        offset = ceil_grid(max(offsets))
        return cone_upper(rho, offset) < w_lower or cone_upper(rho, offset, True) < w_lower

    pending = [(initial_simplex(target), 0)]
    accepted = deepest = visited = 0
    while pending:
        simplex, depth = pending.pop()
        visited += 1
        if closed(simplex):
            accepted += 1
            deepest = max(deepest, depth)
            continue
        if depth >= 36 or visited >= max_nodes:
            raise RuntimeError(f"unresolved simplex after {visited} nodes at depth {depth}: {simplex[0]}")
        pending.extend((child,depth+1) for child in bisect_simplex(simplex))
    return accepted, deepest


if __name__ == "__main__":
    target = F(4, 25)
    assert 24*target - 1 > 0 and (24*target - 1)**2 > 8
    leaves, depth = certify_with_cone(target)
    print(f"seven cells: maximum > {target}, leaves={leaves}, depth={depth}")
