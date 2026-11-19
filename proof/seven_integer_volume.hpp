#ifndef CUBE_SEVEN_INTEGER_VOLUME_HPP
#define CUBE_SEVEN_INTEGER_VOLUME_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <compare>
#include <cstdint>
#include <cstdlib>

namespace seven_integer {
using Int = std::int64_t;
using Wide = __int128;
using Vector = std::array<Int, 3>;
using Bounds = std::array<std::array<Int, 2>, 3>;

// Normals and offsets are expressed in units 1/parameter_grid.
// Spatial bounds are expressed in units 1/space_grid.
constexpr Int parameter_grid = Int{1} << 16;
constexpr Int space_grid = Int{1} << 20;
constexpr Int volume_grid = Int{1} << 24;

enum class Bound { Lower, Upper };

struct Point {
    Int x, y;
    auto operator<=>(const Point &) const = default;
};

struct Line {
    Int a, b, c; // a*X+b*Y <= c
};

struct Height {
    Vector numerator; // z=(a*x+b*y+c)/denominator in physical coordinates
    Int denominator;
};

inline Int turn(Point p, Point q, Point r) { return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x); }

inline Int floor_division(Wide a, Wide b) { return Int(a / b - (a % b < 0)); }

inline Vector difference_numerator(const Height &p, const Height &q) {
    Vector r;
    for (int i = 0; i < 3; ++i)
        r[i] = p.numerator[i] * q.denominator - q.numerator[i] * p.denominator;
    return r;
}

inline bool equal_height(const Height &p, const Height &q) { return difference_numerator(p, q) == Vector{}; }

// Integrate a nonnegative affine function on a rational polygon.
// For a lower bound, the constraints have already been shifted inward;
// for an upper bound, use an outer hull and lift the integrand by a bound
// on its variation within one grid square.
template <Bound bound>
Int polygon_integral(const Line *lines, int count, const Vector &difference, Int denominator) {
    // At most nine lines for four halfspaces. Each pair contributes at most
    // four grid points in the upper-bound calculation.
    Point points[144];
    int length = 0;
    for (int i = 0; i < count; ++i) {
        for (int j = i + 1; j < count; ++j) {
            const Line p = lines[i], q = lines[j];
            Wide det = Wide(p.a) * q.b - Wide(p.b) * q.a;
            if (!det)
                continue;
            Wide x = Wide(p.c) * q.b - Wide(p.b) * q.c;
            Wide y = Wide(p.a) * q.c - Wide(p.c) * q.a;
            if (det < 0) {
                det = -det;
                x = -x;
                y = -y;
            }
            bool feasible = true;
            for (int k = 0; k < count; ++k) {
                if (lines[k].a * x + lines[k].b * y > lines[k].c * det) {
                    feasible = false;
                    break;
                }
            }
            if (feasible) {
                const Int gx = floor_division(x, det), gy = floor_division(y, det);
                points[length++] = {gx, gy};
                if constexpr (bound == Bound::Upper) {
                    points[length++] = {gx + 1, gy};
                    points[length++] = {gx, gy + 1};
                    points[length++] = {gx + 1, gy + 1};
                }
            }
        }
    }
    if (length < 3)
        return 0;
    std::sort(points, points + length);
    length = int(std::unique(points, points + length) - points);
    if (length < 3)
        return 0;
    Point hull[288];
    int m = 0;
    for (int j = 0; j < length; ++j) {
        while (m >= 2 && turn(hull[m - 2], hull[m - 1], points[j]) <= 0)
            --m;
        hull[m++] = points[j];
    }
    const int lower = m;
    for (int j = length - 2; j >= 0; --j) {
        while (m > lower && turn(hull[m - 2], hull[m - 1], points[j]) <= 0)
            --m;
        hull[m++] = points[j];
    }
    const Point p = hull[0];
    auto value = [&](Point a) {
        Wide result =
            Wide(difference[0]) * a.x + Wide(difference[1]) * a.y + Wide(difference[2]) * space_grid;
        if constexpr (bound == Bound::Upper)
            result += std::abs(difference[0]) + std::abs(difference[1]);
        return result;
    };
    const Wide p_value = value(p);
    Wide numerator = 0;
    for (int i = 1; i + 1 < m - 1; ++i) {
        const Point q = hull[i], r = hull[i + 1];
        numerator += Wide(turn(p, q, r)) * (p_value + value(q) + value(r));
    }
    const Wide divisor = Wide(6) * space_grid * space_grid * space_grid * denominator;
    const Wide scaled = numerator * volume_grid;
    if constexpr (bound == Bound::Upper)
        return -floor_division(-scaled, divisor);
    return floor_division(scaled, divisor);
}

// Add p-q >= 0 in the xy plane. Shift the linear inequality before rounding
// candidate vertices down to the spatial grid. Constant inequalities have no
// rounding error and need no buffer.
template <Bound bound> bool append(Line *lines, int &size, const Vector &difference) {
    if (!difference[0] && !difference[1])
        return difference[2] >= 0;
    Int rhs = difference[2] * space_grid;
    if constexpr (bound == Bound::Lower)
        rhs -= std::abs(difference[0]) + std::abs(difference[1]);
    lines[size++] = {-difference[0], -difference[1], rhs};
    return true;
}

// Exact integer enclosure for the intersection of a cube or an octant with
// one to four halfspaces n_i dot x >= h_i. Each normal component and offset must
// lie in [-parameter_grid,parameter_grid]. Coordinates in box must be the
// original cube's boundaries or zero. Returns volume in units 1/volume_grid.
// With these fixed grid sizes, planar line coefficients are <2^34, and all
// intersection predicates, affine integrals, and final divisions fit Wide.
template <Bound bound, std::size_t N = 3>
Int volume(const std::array<Vector, N> &normals, const std::array<Int, N> &offsets, const Bounds &box) {
    static_assert(N >= 1 && N <= 4);
    for (const auto &row : normals)
        for (Int coefficient : row)
            assert(-parameter_grid <= coefficient && coefficient <= parameter_grid);
    for (Int h : offsets)
        assert(-parameter_grid <= h && h <= parameter_grid);
    for (const auto &limits : box) {
        assert(limits[0] <= limits[1]);
        for (Int x : limits)
            assert(x == -space_grid / 2 || x == 0 || x == space_grid / 2);
    }
    const auto [xl, xu] = box[0];
    const auto [yl, yu] = box[1];
    const auto [zl, zu] = box[2];
    Line base[4 + N] = {{1, 0, xu}, {-1, 0, -xl}, {0, 1, yu}, {0, -1, -yl}};
    int nbase = 4;
    Height lower[1 + N] = {{{0, 0, 2 * zl / space_grid}, 2}};
    Height upper[1 + N] = {{{0, 0, 2 * zu / space_grid}, 2}};
    int nl = 1, nu = 1;
    for (std::size_t i = 0; i < N; ++i) {
        const auto [a, b, c] = normals[i];
        if (c == 0) {
            if (!append<bound>(base, nbase, {a, b, -offsets[i]}))
                return 0;
            continue;
        }
        const Int sign = c > 0 ? 1 : -1;
        Height h = {{-sign * a, -sign * b, sign * offsets[i]}, std::abs(c)};
        auto list = sign > 0 ? lower : upper;
        int &length = sign > 0 ? nl : nu;
        bool exists = false;
        for (int j = 0; j < length; ++j)
            if (equal_height(list[j], h))
                exists = true;
        if (!exists)
            list[length++] = h;
    }
    Int result = 0;
    for (int i = 0; i < nl; ++i) {
        for (int j = 0; j < nu; ++j) {
            Line lines[5 + N];
            std::copy(base, base + nbase, lines);
            int size = nbase;
            const Vector difference = difference_numerator(upper[j], lower[i]);
            if (!append<bound>(lines, size, difference))
                continue;
            bool valid = true;
            for (int k = 0; k < nl; ++k)
                if (k != i)
                    valid = valid && append<bound>(lines, size, difference_numerator(lower[i], lower[k]));
            for (int k = 0; k < nu; ++k)
                if (k != j)
                    valid = valid && append<bound>(lines, size, difference_numerator(upper[k], upper[j]));
            if (valid)
                result += polygon_integral<bound>(lines, size, difference,
                                                  upper[j].denominator * lower[i].denominator);
        }
    }
    return result;
}
} // namespace seven_integer
#endif
