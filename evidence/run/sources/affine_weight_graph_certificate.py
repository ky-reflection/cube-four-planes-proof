"""An exact Gram and radial certificate for four affine cube cuts.

The geometry and C++ pair checks are described in affine-weight-graph-proof.md.
The weight cover uses integer polygon clipping on every continuous weight box.
"""

import argparse
from fractions import Fraction as F
from itertools import combinations

from affine_frame_height_certificate import (
    ANCHORS,
    graph_deficit_lower,
    graph_pair_bound,
    verify as height_verify,
)
from affine_gram_bm_certificate import (
    ANCHOR_ANGLE_LOWER,
    PI_LOWER,
    PI_UPPER,
    SQRT_TWO_UPPER,
    angular_quadratic,
    weight_bounds,
)
from four_line_radial_certificate import excess_lower_bound, verify as verify_radial

PAIR_LINE_ANCHORS=((62,388),(65,403),(68,417),(71,430),(74,443),(77,456),(80,468),(83,480),(86,492),(89,503),(92,514),(95,525),(98,536),(101,547),(104,557),(107,568),(110,578),(113,588))
X_LO,X_HI=F('-.12'),F('.445')
LO,HI=F('.23'),F('.572')

def intercept(b, anchors=PAIR_LINE_ANCHORS, *, x_lo=X_LO, x_hi=X_HI, upper_slope=HI):
    """Certify D(x) >= a+b*x by minimizing the angular envelope minus b*x."""
    pair_bounds = [(F(10000-k*k, 10000+k*k), 1-F(p, 1000)) for k, p in anchors]
    pair_bounds.extend((x, d) for x, d in ANCHORS if x_lo <= x <= x_hi)
    by_x = {}
    for x, d in pair_bounds:
        by_x[x] = max(d, by_x.get(x, F(0)))
    cuts = sorted(by_x.items())
    radius = max(abs(x_lo), abs(x_hi), *(abs(x) for x in by_x))
    assert cuts and x_lo <= x_hi and upper_slope > LO
    assert upper_slope**2 * (1-radius**2) > F('.51')**2
    values = []
    if len(cuts) == 1 and x_lo <= cuts[0][0] <= x_hi:
        values.append(cuts[0][1]-b*cuts[0][0])
    for (x0, d0), (x1, d1) in zip(cuts, cuts[1:]):
        left, right = max(x_lo, x0), min(x_hi, x1)
        if left > right:
            continue
        crossing = (d0-d1-LO*x0+upper_slope*x1)/(upper_slope-LO)
        points = [left, right]
        if left < crossing < right:
            points.append(crossing)
        for x in points:
            values.append(max(d0+LO*(x-x0), d1+upper_slope*(x-x1))-b*x)
    x0, d0 = cuts[0]
    if x_lo < x0:
        values.extend(d0+upper_slope*(x-x0)-b*x for x in (x_lo, min(x_hi, x0)))
    x1, d1 = cuts[-1]
    if x_hi > x1:
        values.extend(d1+LO*(x-x1)-b*x for x in (max(x_lo, x1), x_hi))
    return min(values)


S=1<<22
EDGES=tuple(combinations(range(4),2))
SX=(1,0,-1,-1,0,1)
SY=(0,1,-1,-1,1,0)

def ceil_div(a,b):return -(-a//b)

def floor_q(x):return x.numerator//x.denominator

def ceil_q(x):return -floor_q(-x)

def clip(poly, lines, row):
    if not row[0] and not row[1]:
        return poly if row[2]>=0 else []
    l=len(lines)
    lines.append(row)
    a,b,c=row
    inside=[a*x+b*y<=c*d for x,y,d,_ in poly]
    if all(inside):return poly
    if not any(inside):return []
    new=[]
    for i,(x,y,d,e) in enumerate(poly):
        k=(i+1)%len(poly)
        if inside[i]:new.append((x,y,d,e))
        if inside[i]!=inside[k]:
            aa,bb,cc=lines[e]
            den=a*bb-aa*b
            assert den
            xx=c*bb-cc*b
            yy=a*cc-aa*c
            if den<0:xx,yy,den=-xx,-yy,-den
            new.append((xx,yy,den,l if inside[i] else e))
    return new

def tighten(box):
    lo=[a for a,b in box];hi=[b for a,b in box]
    for _ in range(3):
        for i in range(4):
            lo[i]=max(lo[i],S-sum(hi)+hi[i],lo[i-1] if i else 0)
            hi[i]=min(hi[i],S-sum(lo)+lo[i],hi[i+1] if i<3 else S)
            if lo[i]>hi[i]:return None
    return tuple(zip(lo,hi))

def verify_weight_cover(delta, target, node_limit=400000):
    M=F(1,14)+delta
    t,c,u,_,U=height_verify(delta)
    assert -u >= X_LO and c <= X_HI
    gl=[]
    slopes=tuple((b,intercept(b)) for b in map(F,('.30','.33','.34','.345','.348','.35','.355','.36','.38')))
    star_lower=F(0)
    for mask in range(1,64):
        idx=tuple(i for i in range(6) if mask>>i&1)
        graph=tuple(EDGES[i] for i in idx)
        up=graph_pair_bound(M,delta,t,U,graph,0,True)
        g=max(graph_deficit_lower(len(graph),up,c,u),*((len(graph)*(a+b*c)-up)/b for b,a in slopes))
        if mask==7: star_lower=g
        if g>0: gl.append((idx,floor_q(g*S*S)))
    gl.sort(key=lambda kv:F(kv[1],len(kv[0])),reverse=True)
    gl.insert(0,(tuple(range(6)),floor_q((8*target-2+6*c)*S*S)))
    lower,upper=weight_bounds(c)
    for _ in range(8):
        upper=min(upper,c/(1+c)-lower*star_lower/(1+c))
        lower=max(lower,1-3*upper)
    assert lower>0
    lo=floor_q(lower*S);hi=ceil_q(upper*S)
    cnum,cden=c.numerator,c.denominator
    top_c=ceil_q((c+u)*S)
    bvertex=ceil_q(c*c/(4*(1+c))*S)
    vnum=cnum*S;vden=2*(cnum+cden)
    def bval(w):return w*(cnum*S-(cden+cnum)*w)
    stack=[(((lo,hi),)*4,0)]
    visited=0;max_depth=0
    while stack:
        box,depth=stack.pop();visited+=1;max_depth=max(depth,max_depth)
        if visited>node_limit:return False,visited,depth,box
        box=tighten(box)
        if box is None:continue
        weight_b=[]
        for a,b in box:
            first,second=bval(a),bval(b)
            lower_b=max(0,min(first,second)//(cden*S))
            upper_b=bvertex if a*vden<=vnum<=b*vden else ceil_div(max(first,second),cden*S)
            weight_b.append((lower_b,upper_b))
        (b0,b1,b2,b3)=weight_b
        dl=(b0[0]+b1[0]+b2[0]-b3[1])//2
        dh=ceil_div(b0[1]+b1[1]+b2[1]-b3[0],2)
        kl=(0,0,b0[0],dl,b1[0]-dh,b2[0]-dh)
        kh=(0,0,b0[1],dh,b1[1]-dl,b2[1]-dl)
        inverse=[ceil_div(S**3,box[i][0]*box[j][0]) for i,j in EDGES]
        zhi=[ceil_div(top_c*box[i][1]*box[j][1],S*S) for i,j in EDGES]
        lines=[(0,-1,0),(1,0,zhi[0]),(0,1,zhi[1]),(-1,0,0)]
        poly=[(0,0,1,0),(zhi[0],0,1,1),(zhi[0],zhi[1],1,2),(0,zhi[1],1,3)]
        for k in range(2,6):
            poly=clip(poly,lines,(-SX[k],-SY[k],kh[k]))
            if not poly:break
            poly=clip(poly,lines,(SX[k],SY[k],zhi[k]-kl[k]))
            if not poly:break
        for idx,threshold in gl:
            if not poly:break
            sx=sum(inverse[i]*SX[i] for i in idx)
            sy=sum(inverse[i]*SY[i] for i in idx)
            v=sum(inverse[i]*kh[i] for i in idx)
            poly=clip(poly,lines,(-sx,-sy,v-threshold))
        if not poly:continue
        i=max(range(4),key=lambda k:box[k][1]-box[k][0])
        a,b=box[i]
        if b-a<=1:return False,visited,depth,box
        mid=(a+b)//2
        stack.append((box[:i]+((mid,b),)+box[i+1:],depth+1))
        stack.append((box[:i]+((a,mid),)+box[i+1:],depth+1))
    return True,visited,max_depth,None

def verify(delta=F(".0044"), mean_squared=F(".011"), node_limit=400000,
           coupled_radial=False):
    complete, visited, depth, unresolved = verify_weight_cover(delta, mean_squared, node_limit)
    assert complete, f"Unresolved weight box after {visited} nodes: {unresolved}"
    maximum = F(1, 14) + delta
    t, _, upper_dot, _, disagreement = height_verify(delta)
    loss = graph_pair_bound(maximum, delta, t, disagreement, EDGES, 3, False)
    quadratic = angular_quadratic(upper_dot)
    spherical_fraction = (
        6 * ANCHOR_ANGLE_LOWER
        - (6 * SQRT_TWO_UPPER - F(16, 3) * quadratic) * mean_squared
    ) / PI_UPPER - 3
    assert 0 < spherical_fraction < 1
    radius_squared = (1 - spherical_fraction / 3)**-2 - 1
    assert 0 < radius_squared < 1
    extra = F(0)
    if F(9, 16) <= radius_squared <= F(16, 25):
        verify_radial(verbose=False)
        extra = excess_lower_bound(radius_squared)
    if coupled_radial:
        from coupled_radial_certificate import excess_lower_bound as coupled_excess

        if F(8, 17) <= radius_squared <= F(16, 25):
            extra = max(extra, coupled_excess(radius_squared))
    centered = PI_LOWER * radius_squared / 4 + extra
    margin = centered - loss - 6 * maximum
    assert margin > 0, f"The comparison did not close for delta = {delta} and mean bound {mean_squared}."
    print(f"continuous weight cover: {visited} nodes, depth {depth}, no unresolved boxes")
    print(f"squared mean normal < {float(mean_squared):.9f}")
    print(f"hypothetical maximum <= {float(maximum):.12f}")
    print(f"centered 2+2 volume >= {float(centered):.12f}")
    print(f"affine 2+2 volume >= {float(centered - loss):.12f}")
    print(f"contradiction margin >= {float(margin):.12f}")
    return margin


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--delta", type=F, default=F(".0044"))
    parser.add_argument("--mean-squared", type=F, default=F(".011"))
    parser.add_argument("--node-limit", type=int, default=400000)
    parser.add_argument("--coupled-radial", action="store_true")
    args = parser.parse_args()
    verify(args.delta, args.mean_squared, args.node_limit, args.coupled_radial)
