"""Offset confinement near the candidate's first three directions.

Use n_i = sqrt(2) e_i minus the other two coordinate vectors, with each
off-diagonal coefficient allowed to vary independently by at most r.
Write the first three planes n_i dot x = h_i without rescaling.
If the largest final block is at most m = (1 + 2 sqrt(2)) / 48, then
max |h_i| <= 58 r, for 0 <= r <= 1/200. In particular, r=0 forces h=0.
The fourth plane is arbitrary.

Proof. Every cap has volume >= 1-8m, so the computed cap at 1/3 forces
H = max |h_i| < 1/3. Each mixed sign cone at offset zero has volume
at least b = 2m - (17/8)r: starting with the exact candidate, each
inequality loses at most a one-sided slab of width r/sqrt(2).

If the offsets have both signs, choose a mixed sign pattern in which
every constraint relaxes. Compare with translating just a plane of
maximal |h_i| in the same direction. If all offsets have one weak sign,
reflect to make them nonnegative and choose the cone with its lone
positive sign at the smallest offset l. If the maximal offset H is at
another index, this cone contains the one with offsets (l,H,l) in that
order. The latter lies on the offset segment from (0,H,0) to (H,H,H).
When all offsets are equal, just use (H,H,H).

The three endpoint checks at H=1/3 are the two possible relaxations
of one constraint and the common positive offset. Each volume exceeds
L=7/40. Brunn--Minkowski makes the cube root of the volume concave
under these affine offset interpolations. Thus some three-plane cone
has volume at least

    ((1-3H) b^(1/3) + 3H L^(1/3))^3 >= b + H/27.

The last inequality follows by the tangent lower bound for the cube
function and b > 29/200, L-b > 3/200. Since the fourth plane can split
this cone into at most two, H <= 27*(17/8)*r < 58*r for r>0.
"""
from fractions import Fraction as F

from fourth_plane_certificate import common_volume, signed_row
from seven_box_certificate import ONE, ZERO, cap


def check_endpoints():
    lower_a, upper_a = F(1414213, 1000000), F(1414214, 1000000)
    assert lower_a**2 < 2 < upper_a**2
    lower_m, upper_m = (1 + 2 * lower_a) / 48, (1 + 2 * upper_a) / 48
    radius = F(1, 200)
    limit = F(1, 3)
    normals = tuple(tuple((lower_a, upper_a) if i == j else (-ONE-radius, -ONE+radius)
                          for j in range(3)) for i in range(3))
    signs = (1, -1, -1)
    rows = tuple(signed_row(row, sign) for row, sign in zip(normals, signs))
    cap_upper = cap((upper_a, ONE + radius, ONE + radius), limit)
    assert cap_upper < 1 - 8 * upper_m
    endpoint_bounds = tuple(
        common_volume(rows, offsets, precise=True)
        for offsets in ((limit, -limit, -limit), (-limit, ZERO, ZERO),
                        (ZERO, -limit, ZERO))
    )
    target = F(7, 40)
    assert min(endpoint_bounds) > target
    assert 3 / lower_a < F(17, 8)
    assert 2 * lower_m - F(17, 8) * radius > F(29, 200)
    assert 2 * upper_m < F(4, 25)
    assert 3 * F(29, 200) / target * (target - F(4, 25)) > F(1, 27)
    assert 27 * F(17, 8) < 58
    return endpoint_bounds


def check_local_cylinder_constants():
    """Arithmetic for the direction-only neighborhood from the local proof.

    Comparing back to the centered candidate triple loses at most epsilon
    in any one sign cell. If the reference fourth plane is centered, its
    coefficients lambda_i in the first three forms are all positive, or
    some mixed cone of volume 2m is unsplit. Scale max lambda_i to one.

    The coordinates y_i = n_i dot x contain the box [-1/4,1/4]^3, with
    Jacobian |det(n_1,n_2,n_3)| < 7/2. In the cone with lone positive
    index j of largest weight, write t_k=-y_k and t_l=-y_l. On the
    triangle t_k+t_l <= 1/4, lowering the least weight from one gains
    volume at least (1-min lambda)/1344 in the positive piece.

    With all lambda >= 1/2, restrict that triangle to t_k+t_l >= 1/8.
    The appropriate antipodal piece then gains at least
    (3/448)*min(|h_4|,1/16) from shifting the fourth plane in the
    equation sum lambda_i y_i = h_4. Thus |h_4| <= (448/3)*epsilon.
    Normalizing the fourth normal to have coefficient sum three
    multiplies its offset by less than four and bounds the maximum
    coefficient deviation by 10*(1-min lambda). These bounds fit the
    exact twelve-parameter local certificate for the stated radius.
    """
    radius = F(1, 5_000_000_000_000)
    epsilon = 126 * radius
    assert F(17, 8) * 59 < 126
    assert epsilon < F(1, 4000)
    assert 1344 * epsilon < F(1, 2)
    assert epsilon < F(3, 448) / 16
    local_radius = F(1, 2_500_000)
    assert max(radius, 58 * radius, 13440 * epsilon,
               4 * F(448, 3) * epsilon) < local_radius
    return radius


if __name__ == "__main__":
    bounds = check_endpoints()
    print("certified endpoint bounds", *(f"{float(v):.9f}" for v in bounds))
    print("direction-only radius compatible with local certificate",
          check_local_cylinder_constants())
