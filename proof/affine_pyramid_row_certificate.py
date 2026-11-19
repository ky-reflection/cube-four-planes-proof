"""Exact row bounds for the cube-center pyramid exclusion."""

from fractions import Fraction as F
from itertools import product

from centered_bernstein_certificate import bernstein, plus, scale, times


CAP = F(17259, 125000)


def add(*terms):
    result = {}
    for term in terms:
        result = plus(result, term)
    return result


def multiply(*terms):
    result = terms[0]
    for term in terms[1:]:
        result = times(result, term)
    return result


def controls(poly, degrees):
    converted = bernstein(poly, degrees)
    return {index: converted.get(index, F(0))
            for index in product(*(range(degree + 1) for degree in degrees))}


def split(table, axis, degrees):
    """De Casteljau bisection of the tensor Bernstein controls."""
    lower, upper = {}, {}
    degree = degrees[axis]
    axes = [range(d + 1) if i != axis else (0,) for i, d in enumerate(degrees)]
    for base in product(*axes):
        line = [base[:axis] + (j,) + base[axis+1:] for j in range(degree+1)]
        work = [table[index] for index in line]
        for j in range(degree+1):
            lower[line[j]], upper[line[degree-j]] = work[0], work[-1]
            work = [(a+b)/2 for a, b in zip(work, work[1:])]
    return lower, upper


def verify(bound, *, upper_y=F(1), fixed_t=None, x_limit=None, use_cap=False):
    """Bound the loss when x=u*y and the offset magnitude is t*x."""
    dimension = 3 if fixed_t is None else 2
    zero = (0,)*dimension
    one = {zero: F(1)}
    u = {(1,0) + (0,)*(dimension-2): F(1)}
    y = {(0,1) + (0,)*(dimension-2): upper_y}
    t = {(0,0,1): F(1)} if fixed_t is None else {zero: fixed_t}
    x = times(u,y)
    sum_xy = plus(x,y)
    r = plus(one,sum_xy)
    d = multiply(plus(one,u),plus(one,x),plus(one,y))
    offset_factor = add(scale(r,3),
                        scale(multiply(t,r,add(one,u,x)),F(-1,2)),
                        scale(times(t,x),F(1,2)))
    loss_numerator = times(y,add(multiply(r,u,u,plus(one,y)),
                                 times(r,plus(one,x)),
                                 multiply(add(r,one),u,u,t,t,offset_factor)))
    loss_slack = add(scale(times(d,r),6*bound),scale(loss_numerator,-1))
    infeasible = []
    if x_limit is not None:
        infeasible.append(add(x, {zero: -x_limit}))
    if use_cap:
        cap_numerator = multiply(t,x,add(one,u,x,scale(multiply(t,t,u,x),F(-1,3))))
        infeasible.append(add(cap_numerator,scale(d,-CAP)))
    polynomials = [loss_slack, *infeasible]
    degrees = tuple(max(index[i] for poly in polynomials for index in poly)
                    for i in range(dimension))
    stack = [(tuple(controls(poly,degrees) for poly in polynomials),0)]
    visited = 0
    while stack:
        box, depth = stack.pop()
        visited += 1
        if min(box[0].values()) >= 0:
            continue
        if any(min(table.values()) > 0 for table in box[1:]):
            continue
        assert depth < 36, (bound, upper_y, fixed_t, x_limit, depth)
        children = [split(table,depth % dimension,degrees) for table in box]
        stack.extend((child,depth+1) for child in zip(*children))
    return visited


def main():
    upper_m = F(79759, 1000000)
    assert 48*upper_m > 1 and (48*upper_m-1)**2 > 8
    assert CAP == 8*upper_m - F(1,2) > 7*upper_m - F(1,2)
    target = F(1,6) - upper_m
    small_offset = F(869,10000)
    general = F(113,1000)
    one_tiny_row = F(17,2000)
    two_small_rows = F(3,50)
    assert small_offset < target
    assert 3*general + one_tiny_row < 4*target
    assert 2*general + 2*two_small_rows < 4*target
    assert general + two_small_rows + 2*small_offset < 4*target
    cases = (
        ("relative offset at most 307/2500", small_offset,
         {"fixed_t":F(307,2500)}),
        ("smaller nonpivot at most 79/500", small_offset,
         {"fixed_t":F(1),"x_limit":F(79,500)}),
        ("both nonpivots at most 1/60", one_tiny_row,
         {"fixed_t":F(1),"upper_y":F(1,60)}),
        ("both nonpivots at most 27/200", two_small_rows,
         {"fixed_t":F(1),"upper_y":F(27,200)}),
        ("unrestricted row under cap constraint", general,
         {"use_cap":True}),
    )
    for description, bound, parameters in cases:
        print(f"{description}: {verify(bound,**parameters)} boxes, zero unresolved",flush=True)


if __name__ == "__main__":
    main()
