"""Verify the angular bounds in antipodal-vertex-proof.md using rationals."""

from fractions import Fraction as F

from affine_gram_bm_certificate import (
    PI_LOWER,
    cosine_taylor,
    negative_dot_bound,
    pair_error,
    round_up,
    sqrt_down,
)


MAXIMUM = F(2, 25)
G = 7 * MAXIMUM - F(1, 2)
H = 8 * MAXIMUM - F(1, 2)


def centered_upper(affine, first, second, opposite):
    magnitudes = (first, second) if opposite else (first + second,)
    assert affine < 1 - sum(magnitudes)
    upper = min(affine + sum(magnitudes), 1 - max(magnitudes) / 3)
    for _ in range(50):
        q = 1 - upper
        improved = min(upper, round_up(affine + sum(pair_error(s, q) for s in magnitudes)))
        if improved == upper:
            break
        upper = improved
    return upper


def projected_dot_upper(affine, first, second, opposite):
    centered = centered_upper(affine, first, second, opposite)
    dot = negative_dot_bound(centered)
    factor = sqrt_down((1 - F(4, 3) * first**2) * (1 - F(4, 3) * second**2))
    return (dot + F(4, 3) * first * second) / factor


def verify():
    # MAXIMUM is strictly greater than (1+2 sqrt(2))/48.
    assert 48 * MAXIMUM > 1 and (48 * MAXIMUM - 1)**2 > 8
    assert 0 < G < H < F(3, 8)
    nn_min, nn_max = F(1235, 1000), F(2220, 1000)
    np_min, np_max = F(1150, 1000), F(2275, 1000)
    pp_min, pp_max = F(1128, 1000), F(2625, 1000)
    for first, second, opposite, lower, upper in (
        (G, G, False, nn_min, nn_max),
        (G, H, True, np_min, np_max),
        (H, H, False, pp_min, pp_max),
    ):
        lower_dot = -projected_dot_upper(8 * MAXIMUM, first, second, opposite)
        upper_dot = projected_dot_upper(7 * MAXIMUM, first, second, not opposite)
        assert cosine_taylor(lower, 9) > upper_dot
        assert cosine_taylor(upper, 8) < lower_dot
        print(f"Projected angles in [{float(lower):.3f}, {float(upper):.3f}]")
    assert pp_min <= np_min <= nn_min and nn_max <= np_max <= pp_max
    assert 2 * np_max + pp_max < 2 * PI_LOWER + pp_min
    assert np_max - np_min < pp_min


if __name__ == "__main__":
    verify()
