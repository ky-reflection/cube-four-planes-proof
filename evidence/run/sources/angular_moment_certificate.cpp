// Certify the bounds 23/100 and 51/100 on the integral of |v dot x| over
// the central section perpendicular to u, for every orthonormal pair (u,v).
// Build with clang++ -std=c++20 -O3 angular_moment_certificate.cpp.
#include <algorithm>
#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string_view>

using Integer = std::int64_t;
using WideInteger = __int128;

constexpr Integer grid = Integer{1} << 20;
constexpr int depth_limit = 13;

struct Point {
    Integer x;
    Integer y;
    auto operator<=>(const Point &) const = default;
};

struct Line {
    Integer a;
    Integer b;
    Integer c;
};

Integer turn(Point a, Point b, Point c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

WideInteger moment_bound(Integer x, Integer y, Integer z, Integer d, bool upper) {
    // These are columns 2 and 3 of the rotation, with common denominator n.
    // The first coordinate in the section is the direction v.
    const std::array<Integer, 3> b = {2 * (x * y - d * z),
                                      d * d - x * x + y * y - z * z,
                                      2 * (y * z + d * x)};
    const std::array<Integer, 3> c = {2 * (x * z + d * y), 2 * (y * z - d * x),
                                      d * d - x * x - y * y + z * z};
    const Integer n = d * d + x * x + y * y + z * z;
    // Shrink every face by two grid units. Rounding both section coordinates
    // down then gives points inside the original section and its positive half.
    const Integer h = (grid / 2 - (upper ? 0 : 2)) * n;
    const Line lines[] = {
        {b[0], c[0], h}, {-b[0], -c[0], h}, {b[1], c[1], h}, {-b[1], -c[1], h},
        {b[2], c[2], h}, {-b[2], -c[2], h}, {-1, 0, 0}};
    Point points[84];
    int count = 0;
    for (int i = 0; i < 7; ++i)
        for (int j = i + 1; j < 7; ++j) {
            auto p = lines[i], q = lines[j];
            WideInteger det = WideInteger(p.a) * q.b - WideInteger(p.b) * q.a;
            if (det == 0)
                continue;
            WideInteger nx = WideInteger(p.c) * q.b - WideInteger(p.b) * q.c;
            WideInteger ny = WideInteger(p.a) * q.c - WideInteger(p.c) * q.a;
            if (det < 0) {
                det = -det;
                nx = -nx;
                ny = -ny;
            }
            bool feasible = true;
            for (auto w : lines)
                if (w.a * nx + w.b * ny > w.c * det) {
                    feasible = false;
                    break;
                }
            if (feasible) {
                const Integer px = Integer(nx / det - (nx % det < 0));
                const Integer py = Integer(ny / det - (ny % det < 0));
                points[count++] = {px, py};
                // The outer hull contains every original vertex and stays
                // in the half-plane of nonnegative first coordinate.
                if (upper) {
                    points[count++] = {px + 1, py};
                    points[count++] = {px, py + 1};
                    points[count++] = {px + 1, py + 1};
                }
            }
        }
    std::sort(points, points + count);
    count = int(std::unique(points, points + count) - points);
    if (count < 3)
        return 0;
    Point hull[168];
    int length = 0;
    for (int i = 0; i < count; ++i) {
        while (length >= 2 &&
               turn(hull[length - 2], hull[length - 1], points[i]) <= 0)
            --length;
        hull[length++] = points[i];
    }
    int lower = length;
    for (int i = count - 2; i >= 0; --i) {
        while (length > lower &&
               turn(hull[length - 2], hull[length - 1], points[i]) <= 0)
            --length;
        hull[length++] = points[i];
    }
    // The polygon first moment is this sum divided by 6 * grid^3.
    // Central symmetry doubles it to give the integral of the absolute value.
    WideInteger moment = 0;
    for (int i = 1; i < length; ++i) {
        auto p = hull[i - 1], q = hull[i];
        moment += WideInteger(p.x + q.x) * (p.x * q.y - p.y * q.x);
    }
    return moment;
}
Integer excluded = 0, accepted = 0, subdivisions = 0;
int deepest = 0;
bool verify(Integer x, Integer y, Integer z, Integer d, int depth, bool upper) {
    Integer ax = std::max<Integer>(0, std::abs(x) - 1);
    Integer ay = std::max<Integer>(0, std::abs(y) - 1);
    Integer az = std::max<Integer>(0, std::abs(z) - 1);
    if (std::max({ax, ay, az}) * 1000 > 415 * d || ax + ay + az > d) {
        ++excluded;
        return true;
    }
    if (d > 6) {
        auto n = moment_bound(x, y, z, d, upper);
        // The frame moves by at most 2*sqrt(3)/d radians in the box.
        // A (1-6/d) homothety of the central section remains in every section.
        // The reciprocal of that inner scaling also gives an outer
        // scaling. More directly, every moved section sits in the
        // (1+6/d) scaling of the center section.
        auto s = upper ? d + 6 : d - 6;
        auto moment = 100 * n * s * s * s;
        auto denominator = WideInteger(grid) * grid * grid * d * d * d;
        if (upper ? moment <= 153 * denominator : moment >= 69 * denominator) {
            ++accepted;
            deepest = std::max(deepest, depth);
            return true;
        }
    }
    if (depth >= depth_limit) {
        std::cout << "unresolved " << x << " " << y << " " << z << " " << d
                  << "\n";
        return false;
    }
    ++subdivisions;
    for (int dx : {-1, 1})
        for (int dy : {-1, 1})
            for (int dz : {-1, 1})
                if (!verify(2 * x + dx, 2 * y + dy, 2 * z + dz, 2 * d,
                            depth + 1, upper))
                    return false;
    return true;
}
int main(int argc, char** argv) {
    const bool upper = argc == 2 && std::string_view(argv[1]) == "--upper";
    if (argc != 1 && !upper) {
        std::cerr << "Usage: angular_moment_certificate [--upper]\n";
        return 2;
    }
    // The same cube symmetries as in wedge_certificate.cpp cover the chart
    // with these two boxes.
    bool ok = verify(1, 1, -1, 4, 1, upper) && verify(1, 1, 1, 4, 1, upper);
    std::cout << (ok ? "certified" : "failed") << " leaves=" << accepted
              << " excluded=" << excluded << " subdivided=" << subdivisions
              << " depth=" << deepest << "\n";
    return ok ? 0 : 1;
}
