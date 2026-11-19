#ifndef CUBE_FACE_PAIR_MOMENTS_HPP
#define CUBE_FACE_PAIR_MOMENTS_HPP
#include "affine_box_geometry.hpp"
namespace cube_face_pairs {
using namespace cube_affine_box;
using R = Range;
using V = std::array<R, 3>;
using Caps = std::array<R, 4>;
using Box = std::array<R, 8>;
constexpr Int U = scale, Q = (79759 * U + 999999) / 1000000, T = (8452 * U + 999999) / 1000000,
              G = 7 * Q - U / 2;
constexpr R oneU{U, U}, zero{};
R pt(Int a) { return {a, a}; }
R times(R a, int c) { return c < 0 ? -times(a, -c) : R{a.lo * c, a.hi * c}; }
Int mid(R a) { return (a.lo + a.hi) / 2; }
Int minabs(R a) { return std::max({Int(0), a.lo, -a.hi}); }
Int maxabs(R a) { return std::max(-a.lo, a.hi); }
R square(R a) { return {(pt(minabs(a)) * pt(minabs(a))).lo, (pt(maxabs(a)) * pt(maxabs(a))).hi}; }
R unit(R a) { return {std::max(a.lo, -U), std::min(a.hi, U)}; }
R ratio1(R a, R b) {
    if (b.lo <= 0)
        return {-U, U};
    return unit(a / b);
}
R normgrad(Int a, Int b) {
    if (a == 0)
        return {};
    if (a < 0)
        return -normgrad(-a, b);
    if (a < b)
        return pt(a) / pt(6 * b);
    return pt(U / 4) - square(pt(b) / pt(a)) / 12;
}
R grad(R a, R b) {
    R v{U, -U};
    for (auto x : {a.lo, a.hi})
        for (auto y : {minabs(b), maxabs(b)}) {
            auto z = normgrad(x, y);
            v = {std::min(v.lo, z.lo), std::max(v.hi, z.hi)};
        }
    return v;
}
R jvalue(R a, R b) {
    Int x = maxabs(a), y = maxabs(b);
    if (x < y)
        std::swap(x, y);
    if (!x)
        return {};
    return pt(x) / 4 + square(pt(y)) / pt(12 * x);
}
R det(V a, V b, V c) {
    R v{};
    for (int k = 0; k < 3; ++k) {
        int y = (k + 1) % 3, z = (k + 2) % 3;
        v = v + a[k] * (b[y] * c[z] - b[z] * c[y]);
    }
    return v;
}
} // namespace cube_face_pairs
#endif
