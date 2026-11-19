"""Exact rational certificate for the far-from-2 parameter boxes.

The accompanying proof reduces the relevant parameters to:
  7/5 <= k <= 3, |k-2| >= 1/10, 0 <= t1 <= t2 <= t3 <= 1/2,
  w(k).t >= (k+1)|k-2|/2, with w=(k-1,1,1) for k<=2
  and w=(1,1,k-1) for k>=2.
The program verifies that the single-positive region for plane 1 has
volume at least 1/6 throughout this set, using rational arithmetic.
"""
from fractions import Fraction as F
from itertools import combinations

HALF=F(1,2)
ZERO=F(0)
ONE=F(1)
SQUARE=[(-HALF,-HALF),(HALF,-HALF),(HALF,HALF),(-HALF,HALF)]
OFFSETS=[(ZERO,ZERO,ZERO),(ZERO,ZERO,HALF),(ZERO,HALF,HALF),(HALF,HALF,HALF)]

def clip(polygon,line):
    """Keep a*x+b*y+c >= 0, preserving counterclockwise order."""
    if not polygon: return []
    a,b,c=line
    def f(pt): return a*pt[0]+b*pt[1]+c
    prev=polygon[-1]
    fp=f(prev)
    out=[]
    for current in polygon:
        fc=f(current)
        if (fp<0)!=(fc<0):
            weight=fp/(fp-fc)
            out.append(tuple(x+weight*(y-x) for x,y in zip(prev,current)))
        if fc>=0:
            out.append(current)
        prev,fp=current,fc
    return out

def integrate(poly,line):
    if len(poly)<3: return ZERO
    a,b,c=line
    def f(p): return a*p[0]+b*p[1]+c
    x=poly[0]
    v=ZERO
    for y,z in zip(poly[1:],poly[2:]):
        det=(y[0]-x[0])*(z[1]-x[1])-(y[1]-x[1])*(z[0]-x[0])
        v+=det*(f(x)+f(y)+f(z))/6
    return v

def min_positive_integral(lines):
    """Integral on SQUARE of max(0,min of the affine functions)."""
    answer=ZERO
    for i,fn in enumerate(lines):
        poly=clip(SQUARE,fn)
        for j,other in enumerate(lines):
            if i==j: continue
            diff=tuple(b-a for a,b in zip(fn,other))
            if all(x==0 for x in diff):
                if j<i: poly=[]
            else:
                poly=clip(poly,diff)
        answer+=integrate(poly,fn)
    return answer

def singleton_bound(k,t,radius):
    p,q,s=t
    q1=(k,-ONE,HALF-p-radius/2)
    q2=(-ONE,k,HALF-q+radius/2)
    q3=(-ONE,k,HALF-s+radius/2)
    cap=(ZERO,ZERO,ONE)
    return min_positive_integral([cap,q1])-min_positive_integral([cap,q1,q2])-min_positive_integral([cap,q1,q3])

def feasible_vertices(k,radius):
    k0=(k-1,ONE,ONE) if k<=2 else (ONE,ONE,k-1)
    target=(k+1)*abs(k-2)/2-3*radius
    values=[sum(a*b for a,b in zip(k0,p))-target for p in OFFSETS]
    out=[p for p,v in zip(OFFSETS,values) if v>=0]
    for i,j in combinations(range(4),2):
        if values[i]*values[j]<0:
            weight=values[i]/(values[i]-values[j])
            out.append(tuple(a+weight*(b-a) for a,b in zip(OFFSETS[i],OFFSETS[j])))
    return out

if __name__=='__main__':
    count=0
    minimum=None
    radius=F(1,200)
    for box in range(160):
        k=F(7,5)+(2*box+1)*radius
        if F(19,10) < k < F(21,10):
            continue
        for t in feasible_vertices(k,radius):
            bound=singleton_bound(k,t,radius)
            count+=1
            if minimum is None or bound<minimum[0]:
                minimum=(bound,k,t)
            if bound<F(1,6):
                raise RuntimeError(('failed',minimum,float(minimum[0])))
    print("certified offset vertices:", count)
    print("minimum verified bound and its parameters:", minimum)
