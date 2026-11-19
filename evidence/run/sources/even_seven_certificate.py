#!/usr/bin/env python3
"""Check a lower matrix for the seven-cell case with all four even cells active.

Parameters:
    1/7 <= M <= 4/25, p_i >= 0, sum(p_i) = 1.
    s_i = M - (7*M - 1)*p_i
    delta_i = 1/2 - 2*M - s_i
    z_i = max(0, 321/125 - 12*M - 5*s_i)
    g(z) = -1/3 - (943/1000)*z + z*z/7

B_ii = 1-4*delta_i**2, and B_ij = g(z_k)-4*delta_i*delta_j.
The code checks that every parameter box has an entrywise lower bound
given by a symmetric positive definite matrix. All arithmetic is integer.
"""
from itertools import product


def verify():
    n = 20
    scale = 10**8
    base_denominator = 175*n*n
    delta_denominator = 2*base_denominator
    z_denominator = 5*base_denominator
    minimum_det = None
    boxes = 0

    for m_index in range(n):
        m_left = 25*n + 3*m_index
        m_right = m_left + 3
        bounds = {}
        for p_left in range(n+1):
            for p_right in range(p_left, min(p_left+2, n)+1):
                corners = list(product((m_left, m_right), (p_left, p_right)))
                # Delta and z are bilinear in (M,p). Their maxima on a
                # rectangle occur at corners, even when z is truncated at 0.
                delta_num = max(
                    base_denominator - 6*m*n + 2*(7*m-175*n)*p
                    for m, p in corners
                )
                z_num = max(
                    0, max(2247*n*n - 85*m*n + 25*(7*m-175*n)*p for m, p in corners)
                )
                assert 0 <= z_num <= z_denominator
                # g is decreasing on [0,1]. Evaluate it at the upper bound
                # and round each term down, in units of 1/scale.
                g_lower = (
                    (-scale)//3
                    + (-943*scale*z_num)//(1000*z_denominator)
                    + (scale*z_num*z_num)//(7*z_denominator*z_denominator)
                )
                diagonal = scale + (-4*scale*delta_num*delta_num)//(delta_denominator**2)
                bounds[p_left, p_right] = (delta_num, g_lower, diagonal)

        for i in range(n):
            for j in range(n-i):
                d1, g1, a = bounds[i, i+1]
                d2, g2, b = bounds[j, j+1]
                d3, g3, c = bounds[max(0, n-i-j-2), n-i-j]
                x = g3 + (-4*scale*d1*d2)//(delta_denominator**2)
                y = g2 + (-4*scale*d1*d3)//(delta_denominator**2)
                z = g1 + (-4*scale*d2*d3)//(delta_denominator**2)
                det = a*b*c + 2*x*y*z - a*z*z - b*y*y - c*x*x
                assert a > 0 and a*b > x*x and det > 0, (m_index, i, j)
                minimum_det = det if minimum_det is None else min(minimum_det, det)
                boxes += 1

    print(f"verified {boxes} boxes, determinant >= {minimum_det}/{scale**3}")


if __name__ == "__main__":
    verify()
