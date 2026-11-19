// Independent offsets over an eight-dimensional angular neighborhood of the candidate.
#include "affine_box_linear.hpp"
#include <array>
#include <bit>
#include <iostream>
#include <vector>
using namespace cube_affine_box;
namespace {
constexpr int N = 12;
using A = Linear<N>;
using Box = std::array<Range, N>;
using Vec = std::array<A, 3>;
using Forms = std::array<A, 14>;
A cn(Range x) {
    auto a = A::constant(x.lo);
    a.error = {0, x.hi - x.lo};
    a.range = x;
    return a;
}
A mul(const A &a, Int k) {
    auto x = a * A::constant(k * scale);
    return x;
}
Vec cross(const Vec &a, const Vec &b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
A dot(const Vec &a, const Vec &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
int sign(int row, int k) { return row == 0 || row == k + 1 ? 1 : -1; }
bool bit(int m, int j) { return (m >> j) & 1; }
struct Model {
    std::array<std::array<A, 4>, 4> c{}, inverse_den[4]{};

    std::array<std::array<A, 4>, 4> point_factor[3]{};
    std::array<A, 4> cof{};
    Model(const Box &box) {
        constexpr Int root = 1554944255987;
        static_assert(Wide(root) * root < 2 * Wide(scale) * scale &&
                      Wide(root + 1) * (root + 1) > 2 * Wide(scale) * scale);
        const A I = A::constant(scale), R = cn({root, root + 1});
        std::array<A, 8> p;
        for (int k = 0; k < 8; ++k)
            p[k] = A::variable(box[k], k);
        std::array<Vec, 4> n{Vec{I + p[0], I + p[1], I - p[0] - p[1]}, Vec{R, -I - p[6], -I - p[7]},
                             Vec{-I - p[4], R, -I - p[5]}, Vec{-I - p[2], -I - p[3], R}};
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (i != j) {
                    for (int k = 0; k < 3; ++k)
                        c[i][j] = c[i][j] - (sign(j, k) > 0 ? n[i][k] : -n[i][k]);
                    c[i][j].tighten();
                }
        for (int i = 0; i < 4; ++i)
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    if (a != i && b != i && a != b)
                        inverse_den[i][a][b] = I / (c[i][a] + c[i][b]);
        for (int omit = 0; omit < 4; ++omit) {
            int i = (omit + 1) % 4, j = (omit + 2) % 4, k = (omit + 3) % 4;
            auto x = cross(n[j], n[k]), y = cross(n[k], n[i]), z = cross(n[i], n[j]);
            auto d = dot(n[i], x);
            int s = d.range.lo > 0 ? 1 : -1;
            d = s == 1 ? d : -d;
            if (d.range.lo <= 0)
                throw std::logic_error("singular direction interval");
            cof[omit] = d;
            auto reciprocal = I / d;
            if (s < 0)
                reciprocal = -reciprocal;
            for (int t = 0; t < 3; ++t) {
                point_factor[t][omit][i] = x[t] * reciprocal;
                point_factor[t][omit][j] = y[t] * reciprocal;
                point_factor[t][omit][k] = z[t] * reciprocal;
            }
        }
    }
    bool forms(const Box &box, Forms &output) {
        const A unit = A::constant(scale);
        std::array<A, 4> h{};
        int axes[4] = {8, 11, 10, 9};
        A total;
        for (int i = 0; i < 4; ++i) {
            h[i] = mul(A::variable(box[axes[i]], axes[i]), 2);
            total = total + cof[i] * h[i];
        }
        if (total.range.hi < 0)
            return false;
        auto areas = face_expressions(
            unit,
            [&](int f, int i) {
                int a = i ^ face_shift(f), k = 0;
                while (k == i || k == a)
                    ++k;
                int l = 6 - i - a - k;
                auto num = c[i][a] + (bit(f, i) ? h[i] : -h[i]);
                return fraction(num * inverse_den[i][a][k]) * fraction(num * inverse_den[i][a][l]) /
                       2;
            },
            [&](int f, int i, int j) {
                int a = i ^ face_shift(f), b = j ^ face_shift(f);
                auto first = fraction((c[i][a] + h[i]) * inverse_den[i][a][b]);
                auto second = fraction((c[j][b] - h[j]) * inverse_den[j][a][b]);
                auto d = (c[i][a] + c[i][j]) * inverse_den[i][a][b] +
                         (c[j][b] + c[j][i]) * inverse_den[j][a][b];
                return fraction(positive_square(first + second - unit) / d / 2);
            });
        std::array<std::array<A, 4>, 6> point{};
        for (unsigned f = 0; f < 6; ++f) {
            int face = faces[f];
            int a = face_shift(face) - 1;
            for (int b = 0; b < 4; ++b) {
                for (int i = 0; i < 4; ++i)
                    if (i != b)
                        point[f][b] = point[f][b] + point_factor[a][b][i] * h[i];
                if (!bit(face, 0))
                    point[f][b] = -point[f][b];
            }
        }
        for (int mask = 1; mask < 15; ++mask) {
            std::array<int, 3> omit{};
            int n = 0;
            for (int i = 0; i < 4; ++i)
                if (!bit(mask, i))
                    omit[n++] = i;
            auto [i, j, k] = omit;
            A v;
            for (unsigned f = 0; f < 6; ++f) {
                v = v + (unit - point[f][i]) * areas[f][mask];
                if (n > 1)
                    v = v + (point[f][j] - point[f][i]) * areas[f][mask | (1 << i)];
                if (n == 3)
                    v = v + (point[f][j] - point[f][k]) * areas[f][mask | (1 << i) | (1 << j)];
            }
            output[mask - 1] = v / 6;
            output[mask - 1].tighten();
        }
        return true;
    }
};
int direction_den;
bool local(const Box &b) {
    for (int k = 0; k < 12; ++k) {
        if (k < 8 && direction_den >= 100)
            continue;
        auto r = b[k];
        if (std::max(-r.lo, r.hi) * 100 > scale)
            return false;
    }
    return true;
}
struct Solver {
    int max_depth;
    long long budget, nodes = 0, closed = 0, local_boxes = 0;
    bool cover(Box b, Model &mod, int depth) {
        ++nodes;
        if (nodes > budget)
            return false;
        if (local(b)) {
            ++local_boxes;
            return true;
        }
        Forms v;
        if (!mod.forms(b, v) || exclude_linear<N, 14>(v, candidate_upper())) {
            ++closed;
            return true;
        }
        if (depth >= max_depth)
            return false;
        int k = 0;
        for (int i = 1; i < N; ++i)
            if (b[k].hi - b[k].lo < b[i].hi - b[i].lo)
                k = i;
        auto mid = floor_ratio(Wide(b[k].lo) + b[k].hi, 2);
        Box right = b;
        b[k].hi = right[k].lo = mid;
        if (k >= 8)
            return cover(b, mod, depth + 1) && cover(right, mod, depth + 1);
        Model lo(b);
        if (!cover(b, lo, depth + 1))
            return false;
        Model hi(right);
        return cover(right, hi, depth + 1);
    }
};
} // namespace
int main(int argc, char **argv) {
    int radius = argc > 1 ? std::atoi(argv[1]) : 100;
    if (radius < 100 || argc > 4)
        return 2;
    Solver s{argc > 2 ? std::atoi(argv[2]) : 80, argc > 3 ? std::atoll(argv[3]) : 1000000};
    Box box;
    for (int k = 0; k < 12; ++k) {
        Int h = k < 8 ? ceil_ratio(scale, radius) : ceil_ratio(Wide(14) * scale, 100);
        box[k] = {-h, h};
    }
    direction_den = radius;
    Model m(box);
    bool ok = s.cover(box, m, 0);
    std::cout << radius << " " << (ok ? "proved" : "OPEN") << " nodes " << s.nodes << " volumes "
              << s.closed << " local " << s.local_boxes << std::endl;
    return ok ? 0 : 1;
}
