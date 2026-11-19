// Second derivative at zero and radially averaged third-derivative bounds for the local merit.
// Run affine_local_merit_certificate.py on the compiled executable.
#include "affine_box_geometry.hpp"
#include <bit>
#include <cassert>
#include <cstdio>
#include <cstdlib>
using namespace cube_affine_box;
using std::array;
namespace {
constexpr int D = 12;
#include "affine_merit_jet.hpp"
// Integer slopes of the twelve weights, in units of 1/8192.
constexpr int weight_slopes[12][12] = {
    {231, 231, 94, 94, 129, 206, 129, 206, 308, 296, 58, 58},
    {-102, -129, 128, 206, 90, 91, 209, 132, 308, 58, 296, 58},
    {-129, -102, 206, 128, 209, 132, 90, 91, 308, 58, 58, 296},
    {231, 231, 94, 94, 129, 206, 129, 206, -308, -296, -58, -58},
    {-102, -129, 128, 206, 90, 91, 209, 132, -308, -58, -296, -58},
    {-129, -102, 206, 128, 209, 132, 90, 91, -308, -58, -58, -296},
    {-238, 126, -169, -594, -174, -599, -264, -265, -251, -399, -399, -549},
    {126, -238, -594, -169, -264, -265, -174, -599, -251, -399, -549, -399},
    {113, 113, -269, -269, -595, -169, -595, -169, -251, -549, -399, -399},
    {-238, 126, -169, -594, -174, -599, -264, -265, 251, 399, 399, 549},
    {126, -238, -594, -169, -264, -265, -174, -599, 251, 399, 549, 399},
    {113, 113, -269, -269, -595, -169, -595, -169, 251, 549, 399, 399}};
constexpr array<int, 12> masks{3, 5, 9, 12, 10, 6, 7, 11, 13, 8, 4, 2};
// The local chart orders the normals as 0,3,2,1. Geometry uses 0,1,2,3.
int reorder(int m) { return (m & 5) | ((m & 2) << 2) | ((m & 8) >> 2); }
bool bit(int m, int k) { return (m >> k) & 1; }
int sign(int row, int coord) { return (row == 0 || row == 1 + coord) ? 1 : -1; }
struct EdgeTerm {
    Jet delta, slope;
};
Jet merit(const array<Range, D> &box) {
    const Wide two = Wide(2) * scale * scale;
    Int root = 1554944255987;
    while (Wide(root) * root > two)
        --root;
    while (Wide(root + 1) * (root + 1) <= two)
        ++root;
    Range rt{root, root + 1}, wA = reciprocal(rt) / 6, wB = (one - reciprocal(rt)) / 6,
                              m = (one + exact(2) * rt) / 48;
    const Jet I(one), R(rt);
    array<Jet, D> p;
    for (int k = 0; k < D; ++k)
        p[k] = Jet::variable(box[k], k);
    array<V3, 4> n{V3{I + p[0], I + p[1], I - p[0] - p[1]}, V3{R, -I - p[6], -I - p[7]},
                   V3{-I - p[4], R, -I - p[5]}, V3{-I - p[2], -I - p[3], R}};
    array<Jet, 4> offset{p[8] * exact(2), p[11] * exact(2), p[10] * exact(2), p[9] * exact(2)};
    array<array<Jet, 4>, 4> c{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (i != j) {
                for (int k = 0; k < 3; ++k)
                    c[i][j] = c[i][j] - (sign(j, k) > 0 ? n[i][k] : -n[i][k]);
                assert(c[i][j].v.lo > std::max(-offset[i].v.lo, offset[i].v.hi));
            }
    array<V3, 4> point;
    for (int b = 0; b < 4; ++b) {
        int i = (b + 1) % 4, j = (b + 2) % 4, k = (b + 3) % 4;
        auto x = cross(n[j], n[k]), y = cross(n[k], n[i]), z = cross(n[i], n[j]);
        auto inv = inverse(dot(n[i], x));
        for (int t = 0; t < 3; ++t)
            point[b][t] = (x[t] * offset[i] + y[t] * offset[j] + z[t] * offset[k]) * inv;
    }
    array<Jet, 16> weight;
    Jet result;
    for (int z = 0; z < 12; ++z) {
        Jet w(z < 6 ? wA : wB);
        for (int k = 0; k < 12; ++k)
            w = w + p[k] * ratio(weight_slopes[z][k], 8192);
        assert(w.v.lo > 0);
        weight[reorder(masks[z])] = w;
        result = result - w * m;
    }
    // Expand each cell into three-cut parents. Average all available omitted
    // labels to reduce enclosure width. The all-positive cell is empty.
    array<array<Jet, 16>, 4> parent{};
    for (int mask = 1; mask < 15; ++mask) {
        int count = std::popcount(unsigned(15 ^ mask));
        Jet w = weight[mask] / count;
        for (int b = 0; b < 4; ++b)
            if (!bit(mask, b)) {
                parent[b][mask] = w;
                weight[mask | (1 << b)] = weight[mask | (1 << b)] - w;
            }
    }
    array<array<Jet, 4>, 16> heights;
    for (int face : faces) {
        int sh = face_shift(face);
        for (int b = 0; b < 4; ++b)
            heights[face][b] = I - (bit(face, 0) ? point[b][sh - 1] : -point[b][sh - 1]);
        auto factor = [&](int mask) {
            Jet w;
            for (int b = 0; b < 4; ++b)
                w = w + heights[face][b] * parent[b][mask & ~(1 << b)];
            return w;
        };
        auto def = factor(face);
        result = result + def / 6;
        for (int i = 0; i < 4; ++i) {
            int a = i ^ sh;
            int k = (i + 1) % 4;
            while (k == a)
                k = (k + 1) % 4;
            int l = 6 - i - a - k;
            auto t = sq(c[i][a] + (bit(face, i) ? offset[i] : -offset[i])) /
                     ((c[i][a] + c[i][k]) * (c[i][a] + c[i][l]));
            result = result + (factor(face ^ (1 << i)) - def) * t / 12;
        }
    }
    // Join the corrections on the two faces sharing each cube edge before
    // differentiating. The absolute cubic has a continuous second derivative.
    auto edge = [&](int face, int i, int j) {
        int sh = face_shift(face), a = i ^ sh, b = j ^ sh;
        Jet di = inverse(c[i][a] + c[i][b]), dj = inverse(c[j][a] + c[j][b]);
        return EdgeTerm{(c[i][a] + offset[i]) * di + (c[j][b] - offset[j]) * dj - I,
                        (c[i][a] + c[i][j]) * di + (c[j][b] + c[j][i]) * dj};
    };
    for (int f : faces)
        for (int i = 0; i < 4; ++i)
            if (bit(f, i))
                for (int j = 0; j < 4; ++j)
                    if (!bit(f, j)) {
                        int k = f ^ (1 << i) ^ (1 << j);
                        if (k < f)
                            continue;
                        auto t = edge(f, i, j), v = edge(k, j, i);
                        auto di = inverse(t.slope), dj = inverse(v.slope);
                        auto square = sq(t.delta), cubic = abs3(t.delta) * di * dj;
                        for (int b = 0; b < 4; ++b)
                            if (b != i && b != j) {
                                auto w = [&](int mask) { return parent[b][mask & ~(1 << b)]; };
                                auto coef = w(f) + w(k) - w(f ^ (1 << i)) - w(f ^ (1 << j));
                                result = result +
                                         coef *
                                             ((heights[f][b] * di + heights[k][b] * dj) * square / 2 -
                                              cubic) /
                                             12;
                            }
                    }
    return result;
}
} // namespace
int main(int argc, char **argv) {
    Int den = argc > 1 ? std::atoll(argv[1]) : 100;
    if (argc > 2 || den < 25)
        return 2;
    array<Range, D> box{};
    auto center = merit(box);
    assert(center.v.lo <= 0 && center.v.hi >= 0);
    for (auto g : center.d)
        assert(g.lo <= 0 && g.hi >= 0);
    // Integrate the Taylor remainder on sixteen equal radial intervals.
    constexpr Int segments = 16;
    std::array<Wide, T> l{}, u{};
    for (Int step = 1; step <= segments; ++step) {
        Int r = ceil_ratio(Wide(scale) * step, den * segments);
        for (auto &x : box)
            x = {-r, r};
        auto value = merit(box);
        Int a = segments - step + 1, b = a - 1, weight = a * a * a - b * b * b;
        for (int i = 0; i < T; ++i) {
            l[i] += weight * value.t[i].lo;
            u[i] += weight * value.t[i].hi;
        }
    }
    for (int i = 0; i < D; ++i)
        for (int j = 0; j <= i; ++j) {
            auto c = center.h[hh(i, j)];
            std::printf("H %d %d %lld %lld\n", i, j, (long long)c.lo, (long long)c.hi);
        }
    int t = 0;
    for (int i = 0; i < D; ++i)
        for (int j = 0; j <= i; ++j)
            for (int k = 0; k <= j; ++k, ++t)
                std::printf("T %d %d %d %lld %lld\n", i, j, k,
                            (long long)floor_ratio(l[t], segments * segments * segments),
                            (long long)ceil_ratio(u[t], segments * segments * segments));
}
