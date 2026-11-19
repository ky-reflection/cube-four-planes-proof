"""Check the rational comparisons for the lower bound 1/14 + 1/1250.

The two angular anchors are checked by wedge_certificate.cpp, and the
angular derivative bound is checked by angular_moment_certificate.cpp.
"""
from fractions import Fraction as F


def check_fourteen(delta):
    offset_upper = F(1, 14) + 8 * delta
    cap_deficit_lower = F(1, 14) - 6 * delta
    forced_negative_dot = F(1, 2) - 6 * offset_upper**2
    half_angle_tangent = F(61, 100)
    anchor_negative_dot = (1 - half_angle_tangent**2) / (1 + half_angle_tangent**2)
    assert forced_negative_dot > anchor_negative_dot
    centered_overlap_upper = F(384, 2000)
    overlap_upper = centered_overlap_upper * (
        1 - cap_deficit_lower / (3 * centered_overlap_upper)
    ) ** 3
    overlap_lower = F(1, 7) - 12 * delta
    assert overlap_upper < overlap_lower
    return overlap_upper, overlap_lower


def check_fifteen(delta):
    maximum = F(1, 14) + delta
    anchor = F(143, 250)
    slope = F(23, 100)
    gap = F(27, 100)
    assert 8 * maximum + 70 * delta < anchor + slope * gap

    # Quadratic comparison of affine and centered pair disagreement.
    max_abs_dot = F(1, 3) + gap
    sin_squared = 1 - max_abs_dot**2
    quadratic = F(7, 2)
    assert 3 * max_abs_dot**2 < (quadratic - 2)**2 * sin_squared
    assert 3 < quadratic**2 * sin_squared
    assert F(8, 9) < F(943, 1000)**2
    epsilon = F(7, 100)
    assert F(943, 1000) * (
        8 * maximum + quadratic * (70 * delta)**2 - anchor
    ) < slope * epsilon

    # The positive relation bounds the weights and the other dot products.
    r = (F(1, 3) + epsilon) / (F(4, 3) + epsilon)
    mean_normal_norm = 1 - 1 / (4 * r)
    largest_abs_dot = F(41, 100)
    assert F(1, 3) + epsilon < largest_abs_dot
    assert -F(1, 3) + 5 * epsilon + 8 * mean_normal_norm**2 < largest_abs_dot
    det_lower = F(59, 100)
    assert 1 - 3 * largest_abs_dot**2 - 2 * largest_abs_dot**3 > det_lower**2
    weight_product_lower = F(11, 1000)
    assert (1 - 3 * r) * r**2 > weight_product_lower

    # Two substitutions into the tetrahedron volume upper bound.
    central_volume = 28 * delta
    for upper in (F(11, 25000), F(1, 100000)):
        assert (r * (28 * delta + 3 * central_volume))**3 < (
            6 * weight_product_lower * det_lower * upper
        )
        central_volume = upper

    epsilon = F(28, 1000)
    assert 8 * maximum + quadratic * (14 * delta + 2 * central_volume)**2 < (
        anchor + slope * epsilon
    )
    r = (F(1, 3) + epsilon) / (F(4, 3) + epsilon)
    mean_normal_norm = 1 - 1 / (4 * r)
    assert -F(1, 3) + 5 * epsilon + 8 * mean_normal_norm**2 < 0
    assert -F(1, 3) + F(4, 3) * mean_normal_norm**2 < -F(5, 16)
    assert 5 < F(9, 4)**2  # cos(3*pi/5) > -5/16

    # The spherical 2+2 directions have fraction at least 3/5. The cube
    # radial integral is at least 9*pi/64, using pi > 3.14 here.
    centered_mass_lower = 9 * F(314, 100) / 64
    affine_mass_lower = centered_mass_lower - 3 * quadratic * (
        28 * delta + 3 * central_volume
    )**2
    affine_mass_upper = 6 * maximum
    assert affine_mass_lower > affine_mass_upper
    return affine_mass_lower, affine_mass_upper


if __name__ == "__main__":
    delta = F(1, 1250)
    upper, lower = check_fourteen(delta)
    print(f"14 cells: pair overlap at most {float(upper):.9f}, requires {float(lower):.9f}")
    lower, upper = check_fifteen(delta)
    print(f"15 cells: 2+2 mass at least {float(lower):.9f}, allows {float(upper):.9f}")
