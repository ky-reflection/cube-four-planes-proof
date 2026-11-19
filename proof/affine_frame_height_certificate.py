"""Exact height and Gram bounds for the central tetrahedron of four cuts.

This program bounds parameters under a hypothetical largest-block upper bound.
It uses the geometry in affine-frame-height-proof.md. Run with Python 3.
"""

import argparse
from fractions import Fraction as F
from itertools import combinations
from math import prod

from affine_gram_bm_certificate import (
    ANCHORS as STANDARD_ANCHORS,
    GRID,
    LOWER_ANGULAR_SLOPE,
    ONE_THIRD,
    UPPER_ANGULAR_SLOPE,
    mean_norm_squared_bound,
    round_up,
    sqrt_down,
    triangle_area_ratio,
    weight_bounds,
    weight_product_lower_bound,
)

DEFAULT_DELTA = F(62, 10000)
CANDIDATE_UPPER = F(79759, 1000000)
WIDE_SHARP_PAIRS = (
    (12, 7, 371), (7, 4, 366), (200, 113, 363), (20, 11, 357), (15, 8, 350),
    (2, 1, 335), (19, 9, 323), (20, 9, 312), (5, 2, 287), (20, 7, 261),
)
ANCHORS = STANDARD_ANCHORS + tuple(
    (F(a*a - b*b, a*a + b*b), 1 - F(upper, 1000))
    for a, b, upper in WIDE_SHARP_PAIRS
)


def negative_dot_bound(disagreement):
    """Upper bound on the negative dot product from a disagreement upper bound."""
    bounds = []
    for initial_dot, initial_mass in ANCHORS:
        if disagreement <= initial_mass:
            backwards = (initial_mass - disagreement) * sqrt_down(1 - initial_dot**2) / UPPER_ANGULAR_SLOPE
            bounds.append(max(F(0), round_up(initial_dot - backwards)))
            continue

        def excludes(grid_index):
            candidate = F(grid_index, GRID)
            midpoint = (candidate + initial_dot) / 2
            return (LOWER_ANGULAR_SLOPE * (candidate - initial_dot))**2 > (
                (disagreement - initial_mass)**2 * (1 - midpoint**2)
            )

        lower = int(GRID * initial_dot) + 1
        upper = GRID
        if not excludes(upper):
            bounds.append(F(1))
            continue
        while lower < upper:
            middle = (lower + upper) // 2
            if excludes(middle):
                upper = middle
            else:
                lower = middle + 1
        bounds.append(F(lower, GRID))
    return min(bounds)


def fixed_sum_vertices(size, lower, upper, total):
    for count in range(size):
        remainder = total - count * lower - (size - 1 - count) * upper
        if lower <= remainder <= upper:
            yield (lower,)*count + (upper,)*(size - 1 - count) + (remainder,)


def projected_frame_eigen_lower(c,u):
    # The plane containing the projections must support four positive weights.
    # Cyclic gap condition gives 3c_projected-4c_projected^3<=u_projected.
    u=max(F(0),u)
    def excludes(eigenvalue):
        v=c+eigenvalue/2
        return (3*v-u-eigenvalue/2)*(1-eigenvalue)>4*v**3
    if not excludes(F(0)):
        return F(0)
    lower,upper=0,GRID
    while upper-lower>1:
        middle=(lower+upper)//2
        if excludes(F(middle,GRID)):
            lower=middle
        else:
            upper=middle
    return F(lower,GRID)


def frame_det_from_eigen(q,eigen_lower):
    def bound(e):return e*(6-q-4*e+e*e)
    result=eigen_lower**2*(4-2*eigen_lower)
    for cut in (F(1,5),F(1,4),F(3,10),F(2,5),F(1,2),F(3,4),F(1)):
        if cut>=eigen_lower:
            result=max(result,min(bound(eigen_lower),bound(cut),cut**2*(4-2*cut)))
    return result


def frame_scale_lower(c,u,mean_squared,weight_hi):
    s_lo,s_hi = F(-2), -2+8*mean_squared
    vertices = [(-c,)*k + (u,)*(6-k) for k in range(7) if s_lo<=-k*c+(6-k)*u<=s_hi]
    for total in (s_lo,s_hi):
        vertices.extend(fixed_sum_vertices(6,-c,u,total))
    if not vertices:
        return F(0)
    q = max(sum(x*x for x in w) for w in vertices)
    eig_mean = F(4,3)
    var = max(F(0),2*q-F(4,3))
    det = eig_mean**3-eig_mean*var/2-(var/3)*(sqrt_down(var/6)+F(1,GRID))
    det=max(det,frame_det_from_eigen(q,projected_frame_eigen_lower(c,u)))
    full=1//weight_hi
    weight_square=full*weight_hi**2+(1-full*weight_hi)**2
    return sqrt_down(det/weight_square) if det>0 else F(0)


def frame_denominator_lower(c,scale,alpha_floor):
    lo, hi = weight_bounds(c)
    hi -= alpha_floor
    lo = max(lo, 1-3*hi)
    if lo <= 0:
        return F(0)
    weights = tuple(fixed_sum_vertices(4, lo, hi, 1))
    if not weights:
        return F(0)
    return 6*min(prod(w) for w in weights)*scale


def capped_deficit_vertices(maximum,total):
    full = total // maximum
    rest = total-full*maximum
    for masks in combinations(range(1,15),full):
        heavy = tuple((mask,maximum) for mask in masks)
        if rest:
            for last in range(1,15):
                if last not in masks:
                    yield heavy+((last,rest),)
        else:
            yield heavy


def vertex_caps(support,G):
    return tuple(sum(weight for mask,weight in support if mask & (1<<i))-G for i in range(4))


def capped_B_upper(maximum,delta,T,weight_lo,weight_hi):
    G = 7*delta
    offsets = {tuple(sorted(vertex_caps(support,G),reverse=True)) for support in capped_deficit_vertices(maximum,14*delta+T)}
    remaining = 1-4*weight_lo
    weights = []
    for i in range(4):
        extra = min(remaining,weight_hi-weight_lo)
        weights.append(weight_lo+extra)
        remaining -= extra
    assert remaining==0
    return max(sum(w*(e if e>=0 else F(4,7)*e) for w,e in zip(weights,v)) for v in offsets)


def tetra_height_upper(c,scale,H,G,alpha_floor,cap_B):
    lo, hi = weight_bounds(c)
    hi -= alpha_floor
    lo=max(lo,1-3*hi)
    C = H+F(7,8)
    if not (scale>0 and hi>=F(1,4) and H>0):
        return None
    best=F(0)
    bins=64
    for i in range(bins):
        al=lo+(F(1,4)-lo)*i/bins
        ah=lo+(F(1,4)-lo)*(i+1)/bins
        ws=tuple(fixed_sum_vertices(3,al,hi,max(3*al,1-ah)))
        if not ws:
            continue
        prod_lower=min(prod(w) for w in ws)
        if prod_lower<=0:
            return None
        bcap=min(cap_B,H-(H+F(4,7)*G)*al)
        tnum=C**3*ah**2
        if al>0:
            tnum=min(tnum,bcap**3/al)
        best=max(best,tnum/(6*prod_lower*scale))
    return best


def pair_error(sum_deficit, non_disagreement):
    assert sum_deficit>=0 and non_disagreement>=0
    if sum_deficit>=3*non_disagreement:
        return sum_deficit-non_disagreement
    return sum_deficit**2 / (3 * non_disagreement) - sum_deficit**3 / (27 * non_disagreement**2)


def centered_disagreement_upper(affine,H,G,reverse):
    bound=min(F(1),affine+2*H)
    lo,hi=0,GRID
    if 1-2*H>affine:
        while hi-lo>1:
            mid=(lo+hi)//2
            d=F(mid,GRID)
            if d-pair_error(2*H,1-d)>affine:
                hi=mid
            else:
                lo=mid
        bound=min(bound,F(hi,GRID))
    if reverse:
        for _ in range(20):
            q=1-bound
            err=max(pair_error(H+G,q),2*pair_error(H,q))
            nxt=min(bound,round_up(affine+err))
            if nxt==bound:
                break
            bound=nxt
    return bound


STAR = ((0, 1), (0, 2), (0, 3))
FIVE_EDGES = tuple(combinations(range(4), 2))[1:]


def graph_pair_bound(maximum,delta,t_max,disagree_upper,edges,tetra_credit,use_cut):
    G=7*delta
    uniform_error=pair_error(2*(t_max+G),1-disagree_upper)
    same_q=1-8*maximum-uniform_error
    use_affine_q=same_q>0
    if not use_affine_q:
        same_q=1-disagree_upper
    endpoints={F(0),t_max}
    for i in range(1,(14*delta+t_max)//maximum+1):
        t=i*maximum-14*delta
        if 0<t<t_max:
            endpoints.add(t)
    best=None
    for t in sorted(endpoints):
        summaries=set()
        for support in capped_deficit_vertices(maximum,14*delta+t):
            e=vertex_caps(support,G)
            terms=[]
            for i,j in edges:
                a=sum(w for mask,w in support if bool(mask & (1<<i))!=bool(mask & (1<<j)))
                p=max(e[i],F(0))+max(e[j],F(0));n=max(-e[i],F(0))+max(-e[j],F(0))
                terms.append((a,min(p,n),max(p,n)))
            summaries.add(tuple(sorted(terms)))
        for terms in summaries:
            total=-tetra_credit*t
            for a,p,n in terms:
                q=same_q+a if use_affine_q else same_q
                total+=pair_error(p,q)+pair_error(n,q)
                if use_cut:
                    total+=8*maximum-a
            best=total if best is None else max(best,total)
    return best


def graph_deficit_lower(edge_count, centered_sum_upper, negative_dot, positive_dot):
    bounds = [F(0)]
    for c0, D0 in ANCHORS:
        max_abs = max(c0, positive_dot)
        slope = UPPER_ANGULAR_SLOPE / sqrt_down(1 - max_abs**2)
        gain = (edge_count * D0 - centered_sum_upper) / slope - edge_count * max(F(0), c0 - negative_dot)
        bounds.append(gain)
    return max(bounds)


def positive_dot_and_weight_floor(maximum, delta, t_max, disagree_upper, negative_dot, positive_dot):
    if positive_dot >= 1 or disagree_upper >= 1:
        return positive_dot, F(0)
    uniform_error = pair_error(2 * (t_max + 7 * delta), 1 - disagree_upper)
    if 1 - 8 * maximum - uniform_error <= 0:
        return positive_dot, F(0)
    star_centered = graph_pair_bound(maximum, delta, t_max, disagree_upper, STAR, 0, True)
    star_lower = graph_deficit_lower(3, star_centered, negative_dot, positive_dot)
    weight_lower, _ = weight_bounds(negative_dot)
    alpha_floor = weight_lower * star_lower / (1 + negative_dot)
    mean_squared = mean_norm_squared_bound(negative_dot, positive_dot, alpha_floor)
    five_centered = graph_pair_bound(maximum, delta, t_max, disagree_upper, FIVE_EDGES, 0, True)
    five_lower = graph_deficit_lower(5, five_centered, negative_dot, positive_dot)
    positive_dot = min(positive_dot, round_up(-2 + 5 * negative_dot + 8 * mean_squared - five_lower))
    return positive_dot, alpha_floor


def tetra_edge_upper(maximum, tetrahedron, negative_dot, positive_dot, scale):
    """Bound t using the length of the edge shared by the two largest facets."""
    subdivisions = 128
    G = 7*maximum-F(1,2)
    total = 2*G+tetrahedron
    cap_vectors = {tuple(sorted(vertex_caps(support,G)))
                   for support in capped_deficit_vertices(maximum,total)}
    cap_vectors = [v for v in cap_vectors if not any(
        v != w and all(x<=y for x,y in zip(v,w)) for w in cap_vectors)]
    # Outward rounding of every sorted signed-height bound.
    heights = [tuple(int(round_up(max(e,4*e/7))*GRID) for e in vector) for vector in cap_vectors]
    sine = sqrt_down(1-max(negative_dot,positive_dot)**2)
    if not sine:
        return tetrahedron
    sine_units = int(sine*GRID)
    scale_units = int(scale*GRID)
    upper = weight_bounds(negative_dot)[1]
    limit = int(upper*subdivisions)
    best = 0
    def ceil_div(a,b):
        return -(-a//b)
    n = subdivisions
    height_ceiling = int(round_up(G+tetrahedron)*GRID)+875000
    current_ceiling = int(round_up(tetrahedron)*GRID)
    for a in range(n//4+1):
        for b in range(a,min(limit,(n-a)//3)+1):
            for c in range(b,min(limit,(n-a-b)//2)+1):
                # The last weight is uniquely determined; retain its lower bound.
                d = max(c,n-a-b-c-3)
                if d > limit:
                    continue
                # B is measured in units 1/(GRID*n). Each consecutive height
                # difference is nonnegative, so lower cumulative weights suffice.
                B = max(g[3]*n-(g[3]-g[2])*(a+b+c)-(g[2]-g[1])*(a+b)-(g[1]-g[0])*a
                        for g in heights)
                B = min(B,height_ceiling*(a+1))
                if B<0:
                    continue
                bound = current_ceiling
                if c*d*sine_units:
                    # sqrt(3) < 1732051 / GRID.
                    bound = min(bound,ceil_div(1732051*B**2,6*c*d*sine_units*GRID))
                if a*b*c*d*scale_units:
                    bound = min(bound,ceil_div(n*B**3,6*a*b*c*d*scale_units*GRID))
                best = max(best,bound)
    return F(best,GRID)


def verify(delta = DEFAULT_DELTA, use_edges=False):
    M = F(1, 14) + delta
    tetrahedron = min(M,28 * delta)
    assert tetrahedron + 7 * delta < F(2, 7)
    for _ in range(80):
        offset = tetrahedron + 7 * delta
        tetrahedron = min(tetrahedron, round_up(triangle_area_ratio(offset) * offset))
    offset = tetrahedron + 7 * delta
    negative_dot = F(1)
    positive_dot = F(1)
    disagree_upper = F(1)
    for _ in range(20):
        disagree_upper = min(disagree_upper, centered_disagreement_upper(8 * M, offset, 7 * delta, False))
        agreement_cap = min(7*M,6*M+tetrahedron) if use_edges else 7*M
        agree_upper = centered_disagreement_upper(agreement_cap, offset, 7 * delta, True)
        negative_dot = min(negative_dot, negative_dot_bound(disagree_upper))
        positive_dot = min(positive_dot, negative_dot_bound(agree_upper))
        mean_squared = mean_norm_squared_bound(negative_dot, positive_dot)
        positive_dot = min(
            positive_dot,
            round_up(-ONE_THIRD + 5 * (negative_dot - ONE_THIRD) + 8 * mean_squared),
        )
        positive_dot, alpha_floor = positive_dot_and_weight_floor(
            M, delta, tetrahedron, disagree_upper, negative_dot, positive_dot
        )
        lower_weight, upper_weight = weight_bounds(negative_dot)
        absolute_dot = max(negative_dot, positive_dot)
        determinant_squared = 1 - 3 * absolute_dot**2 - 2 * absolute_dot**3
        mean_squared = mean_norm_squared_bound(negative_dot,positive_dot,alpha_floor)
        scale=frame_scale_lower(negative_dot,positive_dot,mean_squared,upper_weight-alpha_floor)
        if use_edges:
            tetrahedron=min(tetrahedron,tetra_edge_upper(M,tetrahedron,negative_dot,positive_dot,scale))
            offset=tetrahedron+7*delta
        cap_B=capped_B_upper(M,delta,tetrahedron,lower_weight,upper_weight-alpha_floor)
        height_upper=tetra_height_upper(negative_dot,scale,offset,7*delta,alpha_floor,cap_B)
        if height_upper is not None:
            tetrahedron=min(tetrahedron,round_up(height_upper))
            offset=tetrahedron+7*delta
        denominator = frame_denominator_lower(negative_dot,scale,alpha_floor)
        if determinant_squared > 0 and lower_weight > 0:
            denominator = max(denominator,6*weight_product_lower_bound(negative_dot)*sqrt_down(determinant_squared))
        if denominator > 0:
            b_upper = min(cap_B,upper_weight*(3*tetrahedron+21*delta), (1-lower_weight)*offset-lower_weight*F(4,7)*7*delta)
            tetrahedron = min(tetrahedron,round_up(b_upper**3/denominator))
            offset=tetrahedron+7*delta

    return tetrahedron, negative_dot, positive_dot, offset, disagree_upper


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("delta", nargs="?", type=F, default=DEFAULT_DELTA)
    parser.add_argument("--candidate", action="store_true", help="bound the proposed optimum")
    parser.add_argument("--edges", action="store_true", help="also use central tetrahedron edge lengths")
    arguments = parser.parse_args()
    if arguments.candidate:
        assert 48 * CANDIDATE_UPPER - 1 > 0
        assert (48 * CANDIDATE_UPPER - 1)**2 > 8
        delta = CANDIDATE_UPPER - F(1, 14)
    else:
        delta = arguments.delta
    tetrahedron, negative_dot, positive_dot, offset, disagreement = verify(delta, arguments.edges)
    if delta == DEFAULT_DELTA:
        assert tetrahedron <= F(1818, GRID)
        assert negative_dot <= F(490800, GRID) and positive_dot <= F(290961, GRID)
        assert offset <= F(45218, GRID) and disagreement <= F(628163, GRID)
    eigenvalue = projected_frame_eigen_lower(negative_dot, positive_dot)
    if arguments.candidate:
        assert negative_dot <= F(740529, GRID) and positive_dot <= F(385202, GRID)
        assert offset <= F(138072, GRID) and disagreement <= F(718536, GRID)
        assert eigenvalue >= F(49451, GRID)
    if arguments.candidate and arguments.edges:
        assert tetrahedron <= F(8452, GRID) and offset <= F(66765, GRID)
        assert negative_dot <= F(573370, GRID) and positive_dot <= F(143814, GRID)
        assert disagreement <= F(654538, GRID) and eigenvalue >= F(251817, GRID)
    print(f"hypothetical maximum <= {float(F(1, 14) + delta):.12f}")
    print(f"tetrahedron volume <= {float(tetrahedron):.6f}")
    print(f"Gram dots in [{-float(negative_dot):.6f}, {float(positive_dot):.6f}]")
    print(f"absolute signed offsets <= {float(offset):.6f}")
    print(f"centered pair disagreement <= {float(disagreement):.6f}")
    print(f"smallest frame eigenvalue > {float(eigenvalue):.6f}")
