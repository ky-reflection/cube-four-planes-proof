"""Certify orientation families for three cuts with one missing sign pattern.

Each normalized normal has one coefficient fixed at 1 or -1. The cuts
are oriented so that +++ is the missing pattern. The command-line demos
prove that every seven-cell arrangement in either family has a cell
larger than 4/25, for arbitrary offsets.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import factorial, isqrt, prod

ZERO = F(0)
ONE = F(1)
HALF = F(1, 2)
TARGET = F(4, 25)
SIGNS = tuple(product((0, 1), repeat=3))


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def cross(x, y):
    return tuple(x[i-2]*y[i-1]-x[i-1]*y[i-2] for i in range(3))


def interval_add(x, y):
    return x[0] + y[0], x[1] + y[1]


def interval_negate(x):
    return -x[1], -x[0]


def interval_multiply(x, y):
    values = tuple(a*b for a, b in product(x, y))
    return min(values), max(values)


def interval_abs(x):
    return max(ZERO, x[0], -x[1]), max(abs(x[0]), abs(x[1]))


def sqrt_upper(value, scale=65536):
    ceiling = (value.numerator * scale**2 + value.denominator - 1) // value.denominator
    lower = isqrt(ceiling)
    return F(lower + (lower**2 < ceiling), scale)

def sqrt_lower(value, scale=65536):
    numerator = value.numerator * scale**2 // value.denominator
    return F(isqrt(numerator), scale)


def cap(weights, h):
    """Volume of n dot x > h for the centered unit cube."""
    w = tuple(a for a in weights if a)
    if not w:
        return ZERO if h >= 0 else ONE
    width = sum(w)/2 - h
    numerator = ZERO
    for pattern in product((0, 1), repeat=len(w)):
        numerator += (-1)**sum(pattern)*max(ZERO, width-dot(pattern,w))**len(w)
    return numerator/(factorial(len(w))*prod(w))


@lru_cache(maxsize=16384)
def cap_quantile(weights, q, steps=15):
    """An enclosing interval for the positive offset yielding volume q."""
    a, b = ZERO, sum(weights)/2
    for _ in range(steps):
        midpoint=(a+b)/2
        if cap(weights,midpoint) > q:
            a=midpoint
        else:
            b=midpoint
    return a, b


def clip(polygon, inequality):
    """Clip an oriented planar polygon by a*x+b*y+c >= 0."""
    if not polygon:
        return []
    a,b,c=inequality
    def value(point):
        return a*point[0]+b*point[1]+c
    previous=polygon[-1]
    prev_value=value(previous)
    answer=[]
    for point in polygon:
        current=value(point)
        if (prev_value<0)!=(current<0):
            scale=prev_value/(prev_value-current)
            answer.append(tuple(p+scale*(q-p) for p,q in zip(previous,point)))
        if current>=0:
            answer.append(point)
        previous,prev_value=point,current
    return answer


def affine_integral(polygon, linear):
    """Integrate a*x+b*y+c over the polygon using signed triangles."""
    if len(polygon)<3:
        return ZERO
    a,b,c=linear
    def value(point):
        return a*point[0]+b*point[1]+c
    x=polygon[0]
    answer=ZERO
    for y,z in zip(polygon[1:],polygon[2:]):
        det=(y[0]-x[0])*(z[1]-x[1])-(y[1]-x[1])*(z[0]-x[0])
        answer+=det*(value(x)+value(y)+value(z))/6
    return answer


def subtract(x,y):
    return tuple(a-b for a,b in zip(x,y))


def volume(bounds, halfspaces):
    """Exact volume in an axis-aligned box intersected with halfspaces.

    Each extra halfspace is a*x+b*y+c*z >= h. Integrate the difference
    between the active upper and lower height functions.
    """
    (xa,xb),(ya,yb),(za,zb)=bounds
    square=[(xa,ya),(xb,ya),(xb,yb),(xa,yb)]
    lower=[(ZERO,ZERO,za)]
    upper=[(ZERO,ZERO,zb)]
    for a,b,c,h in halfspaces:
        if c:
            pair=lower if c>0 else upper
            plane=(-a/c,-b/c,h/c)
            if plane not in pair:
                pair.append(plane)
        else:
            square=clip(square,(a,b,-h))
    result=ZERO
    for lo in lower:
        for hi in upper:
            diff=subtract(hi,lo)
            polygon=clip(square,diff)
            for candidate in lower:
                if candidate!=lo:
                    polygon=clip(polygon,subtract(lo,candidate))
            for candidate in upper:
                if candidate!=hi:
                    polygon=clip(polygon,subtract(candidate,hi))
            result+=affine_integral(polygon,diff)
    return result


CUBE=((-HALF,HALF),)*3
OCTANTS=tuple(tuple((ZERO,HALF) if j else (-HALF,ZERO) for j in p) for p in SIGNS)


def robust_cell_volume(normal_box, offsets, mask, precise):
    """Lower bound on a cell for every normal triple in the parameter box.

    Use either one uniform erosion or the exact common part by octants.
    The common part is convex in the combined space of x and offsets.
    """
    signs=[1 if mask & (1 << i) else -1 for i in range(3)]
    if not precise:
        spaces=[]
        for row,h,s in zip(normal_box,offsets,signs):
            middle=[s*(lo+hi)/2 for lo,hi in row]
            error=sum(hi-lo for lo,hi in row)/4
            spaces.append((*middle,s*h+error))
        return volume(CUBE,spaces)
    answer=ZERO
    for octant in OCTANTS:
        spaces=[]
        for row,h,s in zip(normal_box,offsets,signs):
            coeffs=[]
            for (lo,hi),(a,b) in zip(row,octant):
                coeffs.append(s*(lo if (a==0)==(s==1) else hi))
            spaces.append((*coeffs,s*h))
        answer+=volume(octant,spaces)
    return answer


def vertices(inequalities):
    """Vertices of a bounded polytope given by a dot h <= b."""
    result=set()
    for i,j,k in combinations(range(len(inequalities)),3):
        a,u=inequalities[i]
        b,v=inequalities[j]
        c,w=inequalities[k]
        ab=cross(a,b)
        det=dot(ab,c)
        if not det:
            continue
        bc,ca=cross(b,c),cross(c,a)
        x=tuple((u*p+v*q+w*r)/det for p,q,r in zip(bc,ca,ab))
        if all(dot(row,x)<=rhs for row,rhs in inequalities):
            result.add(x)
    return sorted(result)


def axis_data(normal_box, axis):
    x,y=(axis+1)%3,(axis+2)%3
    cofactors=[]
    for j,k in ((1,2),(2,0),(0,1)):
        cofactors.append(interval_add(
            interval_multiply(normal_box[j][x],normal_box[k][y]),
            interval_negate(interval_multiply(normal_box[j][y],normal_box[k][x])),
        ))
    det=(ZERO,ZERO)
    for cof,row in zip(cofactors,normal_box):
        det=interval_add(det,interval_multiply(cof,row[axis]))
    return cofactors,interval_abs(det)


def positive_part_mean(weights, threshold):
    """Mean positive part of a centered sum of scaled uniform variables."""
    nonzero = tuple(w for w in weights if w)
    degree = len(nonzero) + 1
    distance = sum(nonzero, ZERO) / 2 - threshold
    return sum((-1)**sum(bits) * max(ZERO, distance - dot(bits, nonzero))**degree
               for bits in product((0, 1), repeat=len(nonzero))) / (
                   factorial(degree) * prod(nonzero))


def opposite_coordinate_pair(normal_box, i, j):
    """One convex upper bound on a pair overlap, preserving cancellation."""
    choices = []
    for axis in range(3):
        for first, second in ((i, j), (j, i)):
            a = normal_box[first][axis]
            b = interval_negate(normal_box[second][axis])
            if min(a[0], b[0]) <= 0:
                continue
            combined = tuple(interval_abs(interval_add(
                interval_multiply(b, normal_box[first][k]),
                interval_multiply(a, normal_box[second][k]),
            ))[1] for k in range(3) if k != axis)
            denominator = a[0] * b[0]
            choices.append((positive_part_mean(combined, ZERO) / denominator,
                            first, second, a[0], b[0], combined))
    if not choices:
        return None
    _, first, second, a, b, combined = min(choices)
    return lambda h: positive_part_mean(combined, b*h[first] + a*h[second]) / (a*b)


def possible_pair_masses(normal_box, offsets, target):
    pairs = ((0, 1), (0, 2), (1, 2))
    bounds = [opposite_coordinate_pair(normal_box, i, j) for i, j in pairs]
    if not any(bounds):
        return True
    # The upper bounds and their sums are convex in the offsets. Their
    # maxima on the relaxed offset polytope occur at its vertices.
    values = tuple(tuple(bound(h) if bound else 3*target for bound in bounds)
                   for h in offsets)
    if any(max(v[i] for v in values) < 1 - 6*target for i in range(3)):
        return False
    if max(map(sum, values)) < 1 - 4*target:
        return False
    return all(max(v[i]+v[j] for v in values) >= 1 - 5*target
               for i, j in combinations(range(3), 2))


def offset_polytopes(normal_box, target, axes=(0,1,2)):
    # Monotonicity in coefficient magnitudes follows by conditioning on one
    # uniform variable and using that convolution with a symmetric unimodal
    # density is even and decreasing away from zero.
    magnitudes=tuple(tuple(interval_abs(component) for component in row) for row in normal_box)
    norms=[sqrt_upper(sum(hi*hi for _,hi in row)) for row in magnitudes]
    support=[sum(hi for _,hi in row)/2 for row in magnitudes]
    max_weights=tuple(tuple(hi for _,hi in row) for row in magnitudes)
    core=[]
    upper_offsets=[]
    for i,row in enumerate(magnitudes):
        minimum,maximum=zip(*row)
        lo=cap_quantile(minimum,3*target)[0]
        hi=cap_quantile(maximum,1-4*target)[1]
        # Barthe--Koldobsky, with the normal rescaled to unit length.
        hi=min(hi,norms[i]*(4*target-HALF))
        upper_offsets.append(hi)
        row_vector=tuple(F(i==j) for j in range(3))
        core.append((row_vector,hi))
        core.append((tuple(-q for q in row_vector),-lo))
    # The three cap masses sum to at least 2-5M. Their deficits from 1/2
    # dominate the corresponding Euclidean offsets.
    core.append((tuple(1/q for q in norms),5*target-HALF))
    beta=1-6*target
    for axis in axes:
        cof,det=axis_data(normal_box,axis)
        for sign in (1,-1):
            signed=cof if sign==1 else tuple(interval_negate(v) for v in cof)
            if any(hi<=0 for _,hi in signed):
                continue
            nu_lo=[max(ZERO,lo) for lo,_ in signed]
            nu_hi=[hi for _,hi in signed]
            # Every double-positive block has mass at least beta. Given the
            # missing-sign inequality, it also forces every witness weight
            # to be comparable with the largest weight.
            largest_lower=max(nu_lo)
            nu_lo=[max(lo,beta*largest_lower/(length+hi))
                   for lo,length,hi in zip(nu_lo,support,upper_offsets)]
            if any(lo>hi for lo,hi in zip(nu_lo,nu_hi)):
                continue
            # With all seven cells at most 4/25 the central triple-positive
            # cone has volume w >= 1/200 + (99/1000)*sum(cap deficits).
            # An axis witness gives w <= det^2/(48*nu_1*nu_2*nu_3).
            product_lower=prod(nu_lo)
            min_det=det[0]
            if product_lower:
                threshold=48*product_lower*F(547,50000)
                if det[1]**2<threshold:
                    continue
                min_det=max(min_det,sqrt_lower(threshold))
            # Missing +++ requires this weighted combination to separate.
            necessary=[(tuple(-q for q in nu_hi),-min_det/2)]
            if product_lower:
                max_w=det[1]**2/(48*product_lower)
                deficit_limit=(max_w-F(1,200))*F(1000,99)
                necessary.append((tuple(1/q for q in norms),deficit_limit))
            # Every double-positive cell has volume at least beta. The
            # witness forces it into a slab of bounded width.
            for k in range(3):
                i,j=(k+1)%3,(k+2)%3
                row=tuple(ZERO if p==k else nu_lo[p] for p in range(3))
                bound=det[1]/2+nu_hi[k]*support[k]-beta*max(nu_lo[i],nu_lo[j])
                necessary.append((row,bound))
            offsets = vertices(core+necessary)
            # For positive offsets the cap as a function of its offset is
            # convex. The upper-envelope cap sum takes a maximum at a vertex.
            if offsets and not any(sum(cap(w, h) for w, h in zip(max_weights, v))
                                   >= 2 - 5*target for v in offsets):
                offsets = []
            if offsets and not possible_pair_masses(normal_box, offsets, target):
                offsets = []
            yield axis,sign,offsets


def check_box(normal_box, target=TARGET, axes=(0,1,2), precise=False):
    """Return exact certificates for branches, or None if any is unresolved.

    Restrict axes only when a separate argument guarantees that one of the
    retained axes supplies the missing-sign witness.
    """
    if not isinstance(target, (int, F)):
        raise TypeError("target must be rational")
    target = F(target)
    if not F(1,7)<=target<=TARGET:
        raise ValueError("target must lie in [1/7, 4/25]")
    if len(normal_box)!=3 or any(len(row)!=3 for row in normal_box):
        raise ValueError("expected three normal intervals in each of three rows")
    if any(not isinstance(v, (int, F)) for row in normal_box for pair in row for v in pair):
        raise TypeError("interval endpoints must be rational")
    normal_box = tuple(tuple((F(lo), F(hi)) for lo, hi in row) for row in normal_box)
    if not axes or any(axis not in (0, 1, 2) for axis in axes):
        raise ValueError("expected a nonempty collection of coordinate axes")
    for row in normal_box:
        if any(lo>hi or lo < -1 or hi>1 for lo,hi in row) or not any(lo==hi and abs(lo)==1 for lo,hi in row):
            raise ValueError("normal chart needs a fixed component of magnitude one and other components in [-1,1]")
    proof=[]
    for axis,sign,polytope in offset_polytopes(normal_box,target,axes):
        if not polytope:
            proof.append((axis,sign,None,0,None))
            continue
        for mask in (1,2,4,0,3,5,6):
            certified=target+ONE
            for h in polytope:
                v=robust_cell_volume(normal_box,h,mask,precise)
                certified=min(certified,v)
                if certified<=target:
                    break
            if certified>target:
                proof.append((axis,sign,mask,len(polytope),certified))
                break
        else:
            return None
    return proof


def diagonal_box(lower, upper):
    return tuple(tuple((ONE, ONE) if i == j else (lower, upper)
                       for j in range(3)) for i in range(3))


def verify_diagonal_box(radius, max_depth=16):
    """Bisect a box with all six off-diagonal coefficients independent."""
    axes = tuple((i, j) for i in range(3) for j in range(3) if i != j)
    stack = [(diagonal_box(-HALF - radius, -HALF + radius), 0)]
    closed = deepest = 0
    while stack:
        box, depth = stack.pop()
        proof = check_box(box, precise=True)
        if proof is not None:
            closed += 1
            deepest = max(deepest, depth)
            continue
        if depth == max_depth:
            raise RuntimeError(f"unresolved diagonal box {box}")
        i, j = max(axes, key=lambda index: box[index[0]][index[1]][1]
                   - box[index[0]][index[1]][0])
        lo, hi = box[i][j]
        middle = (lo + hi)/2
        for subinterval in ((lo, middle), (middle, hi)):
            subbox = [list(row) for row in box]
            subbox[i][j] = subinterval
            stack.append((subbox, depth + 1))
    return closed, deepest


def verify_band(noise, max_depth=14):
    """Bisect a common slope and retain six independent perturbations."""
    stack = [(F(-5, 7), F(-1, 3), 0)]
    closed = deepest = 0
    while stack:
        lo, hi, depth = stack.pop()
        proof = check_box(diagonal_box(lo - noise, hi + noise), precise=depth >= 3)
        if proof is not None:
            closed += 1
            deepest = max(deepest, depth)
            continue
        if depth == max_depth:
            raise RuntimeError(f"unresolved band over [{lo}, {hi}]")
        middle = (lo + hi)/2
        stack.extend(((lo, middle, depth + 1), (middle, hi, depth + 1)))
    return closed, deepest


if __name__ == "__main__":
    leaves, depth = verify_diagonal_box(F(1, 7))
    print(f"certified independent radius 1/7: leaves={leaves}, depth={depth}")
    leaves, depth = verify_band(F(1, 35))
    print(f"certified common-slope band radius 1/35: leaves={leaves}, depth={depth}")
