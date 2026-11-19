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
constexpr int D = 8;
#include "affine_angular_merit_weights.hpp"
#include "affine_merit_jet.hpp"
constexpr array<int, 12> masks{3, 5, 9, 12, 10, 6, 7, 11, 13, 8, 4, 2};
// The local chart orders the normals as 0,3,2,1. Geometry uses 0,1,2,3.
int reorder(int m) { return (m & 5) | ((m & 2) << 2) | ((m & 8) >> 2); }
bool bit(int m, int k) { return (m >> k) & 1; }
int sign(int row, int coord) { return (row == 0 || row == 1 + coord) ? 1 : -1; }

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
    array<array<Jet, 4>, 4> c{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (i != j) {
                for (int k = 0; k < 3; ++k)
                    c[i][j] = c[i][j] - (sign(j, k) > 0 ? n[i][k] : -n[i][k]);
                assert(c[i][j].v.lo > 0);
            }
    // Pair-symmetric merit uses seven centered representatives and no triple-point heights.
    std::array<Jet, 7> w{};
    Jet result;
    for (int z = 0; z < 12; ++z) {
        int mask = reorder(masks[z]), opp = 15 ^ mask, index = -1;
        constexpr int reps[]{1, 2, 4, 8, 3, 5, 9};
        for (int r = 0; r < 7; ++r)
            if (mask == reps[r] || opp == reps[r])
                index = r;
        assert(index >= 0);
        Jet value(z < 6 ? wA : wB);
        for (int j = 0; j < 8; ++j)
            value = value + p[j] * ratio(weight_slopes[z][j], 8192);
        assert(value.v.lo > 0);
        w[index] = w[index] + value;
        result = result - value * m;
    }
    for (int r = 4; r < 7; ++r)
        result = result + w[r] / 6;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (i != j) {
                int k = 0;
                while (k == i || k == j)
                    ++k;
                int l = 6 - i - j - k;
                auto g = c[i][j];
                auto triangle = sq(g) / ((g + c[i][k]) * (g + c[i][l]));
                result = result + (w[j] - w[3 + (i ^ j)]) * triangle / 12;
            }
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) {
            int k = 0;
            while (k == i || k == j)
                ++k;
            int l = 6 - i - j - k;
            auto invk = inverse(c[k][i] + c[k][j]), invl = inverse(c[l][i] + c[l][j]);
            auto delta = c[k][i] * invk - c[l][i] * invl;
            auto di = inverse((c[k][i] + c[k][l]) * invk + (c[l][j] + c[l][k]) * invl),
                 dj = inverse((c[l][i] + c[l][k]) * invl + (c[k][j] + c[k][l]) * invk);
            auto weight = w[4] + w[5] + w[6] - w[3 + (i ^ j)] - w[i] - w[j];
            // d_i - d_j = 2*delta. Hence the subtracted absolute-cubic term.
            auto overlap = (di + dj) * sq(delta) / 2 - di * dj * abs3(delta);
            result = result + weight * overlap / 12;
        }
    return result;
}
} // namespace
int main(int argc, char **argv) {
    Int den = argc > 1 ? std::atoll(argv[1]) : 100;
    if (argc > 3 || den < 15)
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
        int s = argc > 2 ? std::atoi(argv[2]) : -1;
        for (int k = 0; k < D; ++k)
            box[k] = s < 0 ? Range{-r, r} : (s >> k & 1) ? Range{-r, 0} : Range{0, r};
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
