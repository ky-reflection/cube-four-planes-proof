// Direction and cap cover for mixed whole-face, three-vertex and star types.
#include "affine_box_linear.hpp"
#include "affine_face_integral.hpp"
#include "affine_section_columns.hpp"
#include <cassert>
#include <iostream>
#include <numeric>
#include <unordered_map>
using namespace cube_affine_box;
#include "affine_face_pair_moments.hpp"
using namespace cube_face_pairs;
constexpr int rep[]{1, 2, 4, 8, 3, 5, 9};
constexpr int pairid[][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
struct Pair {
    R value;
    std::array<R, 8> grad{};
    R ac[3]{};
};
constexpr Int physical = (116 * U + 999) / 1000;
bool use_interpolation = false, use_transport = false, use_simplex = false, use_curvature = false;
struct Type {
    int axis[4]{}, orient[4]{}, bits[4]{}, count[4]{};
    bool face[4]{};
};
Type type;
struct BoxInfo {
    Box b;
    Caps e, H, density;
};
// Symmetric uniform tails. Inputs are fixed-point full projection widths, z is
// twice the threshold.
R half_tail(Int a, Int b, Int z) {
    if (a < b)
        std::swap(a, b);
    if (z >= a + b)
        return zero;
    if (z <= a - b)
        return pt(a - z) / times(pt(a), 2);
    R den = times(pt(a) * pt(b), 8);
    if (den.lo <= 0)
        return {0, U / 2};
    return square(pt(a + b - z)) / den;
}
R band_density(R A, R B, R height) {
    R x = (A - B) / 2, y = (A + B) / 2;
    Int a = std::min(U, maxabs(x)), b = std::min(U, maxabs(y));
    Int h = maxabs(height);
    Int a0 = std::min(U, minabs(x)), b0 = std::min(U, minabs(y));
    Int lower = 0;
    constexpr int N = 4;
    for (int j = 1; j <= N; ++j) {
        // s increasing moves the section outwards, for fixed transverse widths.
        Int shiftlo = 2 * h * j / N, shifthi = ceil_ratio(Wide(2) * h * j, N);
        lower += U - half_tail(a, b, U - shifthi).hi - half_tail(a, b, U + shiftlo).hi;
    }
    return {std::max(745 * U / 1000, lower / N), U - 2 * half_tail(a0, b0, U).lo};
}
bool trim(BoxInfo &o) {
    for (int i = 0; i < 4; ++i) {
        int k = type.axis[i], x = (k + 1) % 3, y = (k + 2) % 3;
        R h{-physical, physical};
        auto &a = o.b[2 * i], &b = o.b[2 * i + 1];
        int s = type.orient[i], bits = type.bits[i];
        if (type.face[i]) {
            a = unit(a);
            b = unit(b);
        }
        for (int sweep = 0; sweep < 3; ++sweep) {
            a = intersect(a, {-2 * U + minabs(b), 2 * U - minabs(b)});
            b = intersect(b, {-2 * U + minabs(a), 2 * U - minabs(a)});
            if (h.lo > h.hi || a.lo > a.hi || b.lo > b.hi)
                return false;
            for (int v = 0; v < 8; ++v) {
                int X = 2 * ((v >> x) & 1) - 1, Y = 2 * ((v >> y) & 1) - 1, Z = 2 * ((v >> k) & 1) - 1;
                bool plus = bits >> v & 1;
                R &p = X == Y ? a : b;
                R n = times(p, Y) + pt(Z * s * U);
                auto val = n / 2;
                if (plus)
                    h.hi = std::min(h.hi, val.hi);
                else
                    h.lo = std::max(h.lo, val.lo);
                R vmin = times(h, 2) - pt(Z * s * U);
                if (plus == (Y > 0))
                    p.lo = std::max(p.lo, Y > 0 ? vmin.lo : -vmin.hi);
                else
                    p.hi = std::min(p.hi, Y > 0 ? vmin.hi : -vmin.lo);
            }
        }
        if (h.lo > h.hi || a.lo > a.hi || b.lo > b.hi)
            return false;
        R density = type.face[i] ? oneU : use_transport ? band_density(a, b, h) : R{745 * U / 1000, U};
        o.density[i] = density;
        o.H[i] = h;
        o.e[i] = intersect({-G, G + T}, density * h);
        if (o.e[i].lo > o.e[i].hi)
            return false;
    }
    return true;
}
Int overflow(R u, R v) {
    Int a = maxabs(u), b = maxabs(v);
    if (a < b)
        std::swap(a, b);
    if (a + b <= U)
        return 0;
    if (!b || a - b >= U)
        return U - (oneU / pt(a)).lo;
    return std::min(U, (square(pt(a + b - U)) / (times(pt(a) * pt(b), 4))).hi);
}
struct PairKey {
    std::array<Int, 8> vals;
    bool operator==(const PairKey &) const = default;
};
struct KeyHash {
    size_t operator()(const PairKey &k) const {
        size_t h = 0;
        for (auto x : k.vals)
            h = h * 131 + size_t(x + 3 * U) / 65536;
        return h;
    }
};
std::unordered_map<PairKey, Pair, KeyHash> pair_tables[6];
Pair mixed_pair(int i, int j, V const *n, V const *center) {
    Pair p;
    std::array<seven_integer::Vector, 2> raw{};
    for (int w = 0; w < 2; ++w) {
        int r = w ? j : i;
        for (int k = 0; k < 3; ++k) {
            auto x = center[r][k];
            assert(x.lo == x.hi && x.lo % (U / 131072) == 0);
            raw[w][k] = x.lo / (U / 131072) * (w ? -1 : 1);
        }
    }
    // For zero heights, six cube-face pyramids give twice the volume of (+,-).
    Wide lo = 0, hi = 0;
    for (int k = 0; k < 3; ++k)
        for (int s : {-1, 1}) {
            cube_face_integral::L line[2];
            for (int w = 0; w < 2; ++w)
                line[w] = {2 * raw[w][(k + 1) % 3], 2 * raw[w][(k + 2) % 3], -s * raw[w][k]};
            auto area = cube_face_integral::projected_area(line, 2);
            lo += area.first;
            hi += area.second;
        }
    constexpr Int denom = 6 * cube_face_integral::coord_grid * cube_face_integral::coord_grid;
    p.value = {floor_ratio(lo * U, denom), ceil_ratio(hi * U, denom)};
    for (int w = 0; w < 2; ++w) {
        int r = w ? j : i, v = w ? i : j, k = type.axis[r], x = (k + 1) % 3, y = (k + 2) % 3;
        auto u = n[r][x], b = n[r][y], t = times(n[v][k], type.orient[r]);
        auto A = n[v][x] - t * u, B = n[v][y] - t * b;
        auto g1 = -grad(A, B), g2 = -grad(B, A);
        Int miss = type.face[r] ? 0 : overflow(u, b);
        // |(X +- Y)/2| <= 1/2 on the square of side one.
        R error{-(miss + 1) / 2, (miss + 1) / 2};
        p.grad[2 * r] = (g1 + g2) / 2 + error;
        p.grad[2 * r + 1] = (g2 - g1) / 2 + error;
    }
    return p;
}

struct TripleKey {
    std::array<Int, 18> values;
    bool operator==(const TripleKey &) const = default;
};
struct TripleHash {
    size_t operator()(const TripleKey &k) const {
        size_t h = 0;
        for (auto x : k.values)
            h = h * 131 + size_t(x + 3 * U) / 65536;
        return h;
    }
};
using OddRows = std::array<std::array<R, 4>, 4>;
std::unordered_map<TripleKey, std::array<R, 4>, TripleHash> slice_cache[4];
struct Strip {
    R a, b, c;
};
R keep_fraction(R p) { return {std::max(Int(0), p.lo), std::min(U, p.hi)}; }
std::array<R, 4> slice_rows(int omit, const BoxInfo &info, const std::array<V, 4> &normals) {
    TripleKey key{};
    int next = 0;
    for (int i = 0; i < 4; ++i)
        if (i != omit) {
            for (int k = 0; k < 2; ++k) {
                auto x = info.b[2 * i + k];
                key.values[next++] = x.lo;
                key.values[next++] = x.hi;
            }
            key.values[next++] = info.H[i].lo;
            key.values[next++] = info.H[i].hi;
        }
    auto &cache = slice_cache[omit];
    auto it = cache.find(key);
    if (it != cache.end())
        return it->second;
    std::array<R, 4> output{};
    constexpr Int grid2 = 1 << 17;
    for (int j = 0; j < 4; ++j)
        if (j != omit) {
            int k = type.axis[j], x = (k + 1) % 3, y = (k + 2) % 3, s = type.orient[j];
            int a = 0;
            while (a == j || a == omit)
                ++a;
            int b = 6 - j - omit - a;
            auto &n = normals[j];
            Strip line[4]{};
            line[0] = {n[x], n[y], {}};
            line[1] = {-n[x], -n[y], {}};
            R h[4];
            h[0] = info.H[j];
            h[1] = -h[0];
            for (int v = 0; v < 2; ++v) {
                int r = v ? b : a;
                auto z = times(normals[r][k], s);
                line[2 + v] = {normals[r][x] - z * n[x], normals[r][y] - z * n[y], {}};
                h[2 + v] = info.H[r] - z * info.H[j];
            }
            Wide mass[2]{};
            constexpr int parts = 2;
            for (int part = 0; part < parts; ++part) {
                R radius{U * part / parts, U * (part + 1) / parts};
                for (int f = 0; f < 4; ++f)
                    line[f].c = radius * h[f] + pt(f < 2 ? -U / 2 : 0);
                for (int m = 0; m < 4; ++m)
                    for (int X : {-1, 1})
                        for (int Y : {-1, 1}) {
                            cube_face_integral::L lines[4];
                            int count = 0;
                            for (int f = type.face[j] ? 2 : 0; f < 4; ++f) {
                                int sg = f < 2 || ((m >> (f - 2)) & 1) ? 1 : -1;
                                Int A = floor_ratio(Wide(times(line[f].a, sg * X).lo) * grid2, U);
                                Int B = floor_ratio(Wide(times(line[f].b, sg * Y).lo) * grid2, U);
                                Int C = ceil_ratio(Wide(times(line[f].c, sg).hi) * grid2, U);
                                lines[count++] = {2 * A, 2 * B, 4 * C - A - B};
                            }
                            mass[m == 0 || m == 3] += cube_face_integral::projected_area(lines, count).first;
                        }
            }
            constexpr Int denom =
                2 * cube_face_integral::coord_grid * cube_face_integral::coord_grid * 4 * parts;
            auto f = [&](int s) { return std::max(Int(0), floor_ratio(mass[s] * U, denom)); };
            output[j] = {f(1), U - f(0)};
            assert(output[j].lo <= output[j].hi);
        }
    if (cache.size() > 60000)
        cache.clear();
    cache.emplace(key, output);
    return output;
}

#include "affine_mixed_interpolation.hpp"

bool use_pair_caps = true, use_sections = true;
constexpr Int mass = U;
template <int K> std::array<int, K> proposal(const std::vector<std::array<Int, K>> &constraints) {
    constexpr int lim = 38, cols = K + 1 + lim + 1;
    double a[lim + 1][cols]{};
    int basis[lim]{};
    int nr = 1 + (int)constraints.size();
    assert(nr <= lim);
    int rhs = K + 1 + nr;
    for (int r = 0; r < nr; ++r) {
        basis[r] = K + 1 + r;
        a[r][basis[r]] = 1.;
    }
    a[0][rhs] = 1.;
    for (int k = 0; k < K; ++k)
        a[0][k] = 1.;
    for (int j = 1; j < nr; ++j) {
        a[j][K] = 1.;
        for (int k = 0; k < K; ++k)
            a[j][k] = -double(constraints[j - 1][k]) / mass;
    }
    a[nr][K] = -1.;
    for (int step = 0; step < 180; ++step) {
        int enter = 0;
        while (enter < rhs && a[nr][enter] >= -1e-10)
            ++enter;
        if (enter == rhs) {
            std::array<int, K> w{};
            for (int r = 0; r < nr; ++r)
                if (basis[r] < K) {
                    double t = a[r][rhs];
                    if (!std::isfinite(t) || t > 2.)
                        return {};
                    w[basis[r]] = std::lround(32768 * std::max(0., t));
                }
            if (std::accumulate(w.begin(), w.end(), 0) > 65536)
                return {};
            return w;
        }
        int leave = -1;
        double best = 1e99;
        for (int r = 0; r < nr; ++r)
            if (a[r][enter] > 1e-10) {
                auto v = std::max(0., a[r][rhs]) / a[r][enter];
                if (v < best - 1e-10 ||
                    (std::abs(v - best) <= 1e-10 && leave >= 0 && basis[r] < basis[leave])) {
                    leave = r;
                    best = v;
                }
            }
        if (leave < 0)
            return {};
        double inv = 1. / a[leave][enter];
        for (int k = 0; k <= rhs; ++k)
            a[leave][k] *= inv;
        basis[leave] = enter;
        for (int r = 0; r <= nr; ++r)
            if (r != leave) {
                double u = a[r][enter];
                if (!u)
                    continue;
                for (int k = 0; k <= rhs; ++k)
                    a[r][k] -= u * a[leave][k];
                a[r][enter] = 0.;
            }
    }
    return {};
}

#include "affine_pair_transfer.hpp"
struct Model {
    PairTransfer transfer[6];
    std::array<std::array<Int, 7>, 256> latticeVertex{};
    R latticeErr[6]{};
    bool haveLattice = false;
    OddRows prob{};
    bool haveOdd = false;
    std::array<V, 4> n{};
    Caps cap, lambda, hlimit, density;
    R den{};
    Pair pairs[6];
    R D0[6];
    bool use = true;
    std::array<Int, 128> lower{};
    R offset[8];
    int direction = -1;
    std::array<std::array<Int, 8>, 128> slopes{};
    static inline int co[128][7]{};
    static void init() {
        for (int m = 1; m < 128; ++m) {
            int r = std::countr_zero((unsigned)m);
            std::copy(co[m & (m - 1)], co[m & (m - 1)] + 7, co[m]);
            co[m][6] += r < 4 ? 1 : -1;
            for (int h = 0; h < 6; ++h) {
                auto [i, j] = pairid[h];
                co[m][h] += r < 4 ? -(r != i && r != j) : int(((rep[r] >> i) ^ (rep[r] >> j)) & 1);
            }
        }
    }
    Model(const BoxInfo &info, int depth) : hlimit(info.H), density(info.density) {
        const auto &box = info.b;
        V centers[4];
        for (int i = 0; i < 4; ++i) {
            R a = box[2 * i], b = box[2 * i + 1];
            int k = type.axis[i];
            n[i][k] = pt(type.orient[i] * U);
            n[i][(k + 1) % 3] = (a - b) / 2;
            n[i][(k + 2) % 3] = (a + b) / 2;
            cap[i] = info.e[i];
            centers[i][k] = n[i][k];
            centers[i][(k + 1) % 3] = (pt(mid(a)) - pt(mid(b))) / 2;
            centers[i][(k + 2) % 3] = (pt(mid(a)) + pt(mid(b))) / 2;
        }
        for (int h = 0; h < 6; ++h) {
            auto [i, j] = pairid[h];
            auto &p = pairs[h];
            V gr[2]{};
            int k = type.axis[i], l = type.axis[j];
            PairKey key;
            int nkey = 0;
            for (int r : {i, j})
                for (int a = 0; a < 2; ++a) {
                    auto v = box[2 * r + a];
                    key.vals[nkey++] = v.lo;
                    key.vals[nkey++] = v.hi;
                }
            auto &table = pair_tables[h];
            auto it = table.find(key);
            if (it != table.end()) {
                p = it->second;
                continue;
            }
            if (!type.face[i] || !type.face[j]) {
                p = mixed_pair(i, j, n.data(), centers);
                if (table.size() > 50000)
                    table.clear();
                table.emplace(key, p);
                continue;
            }
            if (k == l) {
                int x = (k + 1) % 3, y = (k + 2) % 3;
                R a = n[i][x] + n[j][x], b = n[i][y] + n[j][y];
                // All coordinates have dyadic endpoints and exact midpoints here.
                auto A = centers[i][x] + centers[j][x], B = centers[i][y] + centers[j][y];
                p.value = oneU - jvalue(A, B);
                gr[0][x] = gr[1][x] = -grad(a, b);
                gr[0][y] = gr[1][y] = -grad(b, a);
                p.ac[0] = {std::max(minabs(a), minabs(b)), std::max(maxabs(a), maxabs(b))};
            } else {
                int t = 3 - k - l;
                auto load = [&](auto &n) {
                    return std::array<R, 4>{times(n[i][l], -type.orient[j]), -n[i][t],
                                            times(n[j][k], -type.orient[i]), -n[j][t]};
                };
                auto [A, B, C, D] = load(centers);
                R div = oneU - A * C;
                if (div.lo <= 0)
                    p.value = {0, U};
                else
                    p.value = pt(U / 2) + (A + C) / 4 -
                              (C * square(B) + times(B * D, 2) + A * square(D)) / (times(div, 12));
                auto [a, b, c, d] = load(n);
                R z = oneU - a * c;
                auto g1 = ratio1(c * b + d, z), g2 = ratio1(a * d + b, z);
                gr[0][l] = times(pt(U / 4) - square(g1) / 12, -type.orient[j]);
                gr[0][t] = g1 / 6;
                gr[1][k] = times(pt(U / 4) - square(g2) / 12, -type.orient[i]);
                gr[1][t] = g2 / 6;
                p.ac[0] = a;
                p.ac[1] = c;
                p.ac[2] = z;
            }
            for (int w = 0; w < 2; ++w) {
                int v = w ? j : i, ax = type.axis[v];
                R x = gr[w][(ax + 1) % 3], y = gr[w][(ax + 2) % 3];
                p.grad[2 * v] = (x + y) / 2;
                p.grad[2 * v + 1] = (y - x) / 2;
            }
            if (table.size() > 50000)
                table.clear();
            table.emplace(key, p);
        }
        for (int i = 0; i < 4; ++i) {
            auto a = det(n[(i + 1) % 4], n[(i + 2) % 4], n[(i + 3) % 4]);
            lambda[i] = i % 2 ? -a : a;
        }
        for (int j = 0; j < 8; ++j)
            offset[j] = box[j] - pt(mid(box[j]));
        if (use_pair_caps)
            for (int h = 0; h < 6; ++h) {
                auto [i, j] = pairid[h];
                R dot{}, ni{}, nj{};
                for (int k = 0; k < 3; ++k) {
                    dot = dot + n[i][k] * n[j][k];
                    ni = ni + square(n[i][k]);
                    nj = nj + square(n[j][k]);
                }
                auto sqdot = square(dot), norm = ni * nj;
                if (dot.lo > 0 && Wide(sqdot.lo) * 1000000000000 > Wide(norm.hi) * 143814 * 143814)
                    use = false;
                if (dot.hi < 0 && Wide(sqdot.lo) * 1000000000000 > Wide(norm.hi) * 573370 * 573370)
                    use = false;
                auto &p = pairs[h];
                D0[h] = p.value;
                for (int k = 0; k < 8; ++k)
                    D0[h] = D0[h] + p.grad[k] * offset[k];
                if (D0[h].hi * 1000000 < 499000 * U || D0[h].lo * 1000000 > 654538 * U)
                    use = false;
            }
        if (use_sections && use && depth >= 9) {
            haveOdd = true;
            for (int i = 0; i < 4; ++i)
                prob[i] = slice_rows(i, info, n);
        }
        if (!use)
            return;
        if (use_curvature && depth >= 8)
            for (int h = 0; h < 6; ++h) {
                auto [i, j] = pairid[h];
                transfer[h] = transfer_coefficients(i, j, n, hlimit, density);
            }

        for (int m = 1; m < 128; ++m) {
            R v = pt(co[m][6] * U);
            std::array<R, 8> g{};
            for (int h = 0; h < 6; ++h)
                if (int c = co[m][h]) {
                    v = v + times(pairs[h].value, c) / 2;
                    for (int j = 0; j < 8; ++j)
                        g[j] = g[j] + times(pairs[h].grad[j], c) / 2;
                }
            for (int j = 0; j < 8; ++j) {
                v = v + g[j] * offset[j];
                slopes[m][j] = mid(g[j]);
            }
            lower[m] = v.lo;
        }
        if (use_interpolation && depth >= 4) {
            std::array<PairLattice, 6> lattice;
            R err[6]{};
            for (int h = 0; h < 6; ++h) {
                auto [i, j] = pairid[h];
                auto &p = pairs[h];
                auto &f = lattice[h];
                f = face_lattice(h, box, n);
                for (int w = 0; w < 2; ++w)
                    for (int z = 0; z < 2; ++z) {
                        int k = 2 * (w ? j : i) + z;
                        Int r = maxabs(offset[k]);
                        R g = p.grad[k];
                        Int lip = ceil_ratio(Wide(g.hi - g.lo) * r, 2 * U);
                        R bounds{-lip, lip};
                        if (f.smooth) {
                            R bend = square(pt(r)) * f.curvature[w][z] / 2;
                            bounds =
                                intersect(bounds, {-std::max(Int(0), bend.hi), -std::min(Int(0), bend.lo)});
                        }
                        err[h] = err[h] + bounds;
                    }
            }
            haveLattice = true;
            std::copy(err, err + 6, latticeErr);
            std::array<Int, 128> minima;
            minima.fill(std::numeric_limits<Int>::max());
            minima[0] = 0;
            for (int m = 0; m < 256; ++m) {
                R D[6];
                for (int h = 0; h < 6; ++h) {
                    auto [i, j] = pairid[h];
                    D[h] = lattice[h].corner[((m >> (2 * i)) & 3) + 4 * ((m >> (2 * j)) & 3)];
                }
                Int atom[7];
                for (int r = 0; r < 7; ++r) {
                    R f = pt((r < 4 ? 1 : -1) * U);
                    for (int h = 0; h < 6; ++h)
                        if (int c = co[1 << r][h])
                            f = f + times(D[h], c) / 2;
                    atom[r] = f.lo;
                }
                for (int r = 0; r < 7; ++r)
                    latticeVertex[m][r] = 6 * atom[r];
                Int z[128]{};
                for (int mask = 1; mask < 128; ++mask) {
                    z[mask] = z[mask & (mask - 1)] + atom[std::countr_zero(unsigned(mask))];
                    minima[mask] = std::min(minima[mask], z[mask]);
                }
            }
            for (int m = 1; m < 128; ++m) {
                Int v = minima[m];
                for (int h = 0; h < 6; ++h)
                    v += (times(err[h], co[m][h]) / 2).lo;
                lower[m] = std::max(lower[m], v);
            }
        }
        if (!use_interpolation)
            return;
        Int best = std::numeric_limits<Int>::min(), maxwidth = 0;
        for (auto &x : box)
            maxwidth = std::max(maxwidth, x.hi - x.lo);
        for (int k = 0; k < 8; ++k) {
            Int width = box[k].hi - box[k].lo;
            if (4 * width < maxwidth || width <= U / 2048)
                continue;
            Int hi = -U, lo = -U;
            for (int m = 1; m < 127; ++m) {
                int count = std::popcount(unsigned(m));
                Int base = lower[m] - 2 * count * Q, shift = (pt(slopes[m][k]) * pt(width / 2)).lo;
                lo = std::max(lo, (base - std::min(Int(0), shift)) / count);
                hi = std::max(hi, (base + std::max(Int(0), shift)) / count);
            }
            Int score = 3 * std::min(hi, lo) + 2 * std::max(hi, lo);
            if (score > best) {
                best = score;
                direction = k;
            }
        }
    }
    bool init_feasible() {
        if (!use)
            return false;
        bool pos = true, neg = true;
        for (auto v : lambda) {
            pos &= v.hi > 0;
            neg &= v.lo < 0;
        }
        if (!pos && !neg)
            return false;
        den = pt(6 * U);
        for (auto &v : lambda) {
            v = {minabs(v), maxabs(v)};
            den = den * v;
        }
        return true;
    }
    // Center/deficit comparison, with the cap and section constraints as
    // additional rows.
#include "affine_mixed_transport.hpp"
    bool possible(Caps &e) const {
        for (int sweep = 0; sweep < 2; ++sweep) {
            R massT{0, T}, B{};
            for (int j = 0; j < 4; ++j)
                B = B + lambda[j] * intersect(hlimit[j], e[j] / density[j]);
            if (B.hi < 0)
                return false;
            B = nonnegative(B);
            if (den.lo > U / 10000) {
                massT = intersect(massT, (square(B) * B) / den);
                if (massT.lo > massT.hi)
                    return false;
            }
            R corr[6];
            for (int h = 0; h < 6; ++h) {
                auto [i, j] = pairid[h];
                const auto &p = pairs[h];
                R bm{-ceil_ratio(Wide(square(e[i] - e[j]).hi) * 1000000, 3 * 379464),
                     ceil_ratio(Wide(square(e[i] + e[j]).hi) * 1000000, 3 * 345462)};
                if (!type.face[i] || !type.face[j])
                    corr[h] = bm;
                else if (type.axis[i] == type.axis[j]) {
                    R z = {0, bm.hi};
                    if (p.ac[0].lo > U / 1000)
                        z.hi = std::min(z.hi, (square(e[i] + e[j]) / p.ac[0]).hi);
                    corr[h] = z;
                } else if (p.ac[2].lo > U / 1000)
                    corr[h] = intersect(
                        bm,
                        (p.ac[1] * square(e[i]) + times(e[i] * e[j], 2) + p.ac[0] * square(e[j])) / p.ac[2]);
                else
                    corr[h] = bm;
                if (transfer[h].valid)
                    corr[h] = intersect(corr[h], transfer[h].at(e[i], e[j]));
                if (corr[h].lo > corr[h].hi)
                    return false;
                if (use_pair_caps) {
                    auto D = D0[h] - corr[h];
                    auto sums = e[i] + e[j], diff = e[i] - e[j];
                    if (D.lo + minabs(diff) > 8 * Q ||
                        D.hi + std::min(sums.hi, 2 * massT.hi - sums.lo) < U - 6 * Q)
                        return false;
                    R sum = intersect(sums, {U - 6 * Q - D.hi, D.hi + 6 * Q + 2 * massT.hi - U});
                    R dif = intersect(diff, {D.lo - 8 * Q, 8 * Q - D.lo});
                    if (sum.lo > sum.hi || dif.lo > dif.hi)
                        return false;
                    e[i] = intersect(e[i], intersect(sum - e[j], dif + e[j]));
                    e[j] = intersect(e[j], intersect(sum - e[i], e[i] - dif));
                    if (e[i].lo > e[i].hi || e[j].lo > e[j].hi)
                        return false;
                }
            }
            Int slack[128]{}, lo[128]{};
            for (int m = 1; m < 128; ++m) {
                R v = times(massT, co[m][6]);
                for (int h = 0; h < 6; ++h)
                    v = v + times(corr[h], co[m][h]) / 2;
                lo[m] = lower[m] - v.hi;
                slack[m] = 2 * std::popcount(unsigned(m)) * Q - lo[m];
                if (slack[m] < 0)
                    return false;
            }
            for (int m = 1; m < 128; ++m)
                slack[m] = std::min(slack[m], U - massT.lo - lo[127 ^ m]);
            if (haveOdd) {
                std::array<std::array<R, 4>, 7> odd{};
                for (int i = 0; i < 4; ++i)
                    odd[i] = prob[i];
                for (int r = 4; r < 7; ++r)
                    for (int j = 0; j < 4; ++j) {
                        R u{};
                        for (int i = 0; i < 4; ++i) {
                            int sg = rep[r] >> i & 1 ? -1 : 1;
                            u = u + times(odd[i][j] + (i == j ? oneU : zero), sg);
                        }
                        odd[r][j] = u / 2;
                    }
                Int abslo[7]{};
                for (int r = 0; r < 7; ++r) {
                    R a = r < 4 ? -massT : zero;
                    for (int j = 0; j < 4; ++j)
                        a = a + odd[r][j] * e[j];
                    abslo[r] = minabs(a);
                    Int room = slack[1 << r];
                    if (abslo[r] > room)
                        return false;
                    for (int j = 0; j < 4; ++j) {
                        auto c = odd[r][j];
                        if (minabs(c) < U / 4096)
                            continue;
                        auto term = c * e[j];
                        R other{a.lo - term.lo, a.hi - term.hi};
                        auto rhs = R{-room, room} - other;
                        auto v = intersect(e[j], c.lo > 0 ? rhs / c : (-rhs) / (-c));
                        if (v.lo > v.hi)
                            return false;
                        e[j] = v;
                        a = other + c * v;
                    }
                }
                Int sum[128]{};
                for (int m = 1; m < 128; ++m) {
                    sum[m] = sum[m & (m - 1)] + abslo[std::countr_zero((unsigned)m)];
                    if (sum[m] > slack[m])
                        return false;
                }
            }
            for (int w = 0; w < 81; ++w) {
                int z = w, a[4], sum = 0;
                R v{};
                for (int j = 0; j < 4; ++j) {
                    a[j] = z % 3 - 1;
                    z /= 3;
                    sum += a[j];
                    v = v + times(e[j], -2 * a[j]);
                }
                v = v + times(massT, sum);
                int masks[4]{};
                for (int r = 0; r < 7; ++r) {
                    int g = 0;
                    for (int j = 0; j < 4; ++j)
                        g += a[j] * ((rep[r] >> j & 1) ? 1 : -1);
                    for (int k = 0; k < std::abs(g); ++k)
                        masks[k] |= 1 << r;
                }
                Int room = 0;
                for (int m : masks)
                    room += slack[m];
                if (v.lo > room || v.hi < -room)
                    return false;
                for (int j = 0; j < 4; ++j)
                    if (a[j]) {
                        R x = times(e[j], -2 * a[j]), other{v.lo - x.lo, v.hi - x.hi};
                        auto b = times(R{-room - other.hi, room - other.lo}, -a[j]) / 2;
                        auto n = intersect(e[j], b);
                        if (n.lo > n.hi)
                            return false;
                        v = other + times(n, -2 * a[j]);
                        e[j] = n;
                    }
            }
        }
        return true;
    }
};
long long visits = 0, capnodes = 0, contractions = 0;
int maxdepth = 0;
bool caps(const Model &m, Caps e, int depth) {
    ++capnodes;
    if (!m.possible(e))
        return true;
    if (!depth)
        return false;
    int ax = 0;
    for (int i = 1; i < 4; ++i)
        if (e[i].hi - e[i].lo > e[ax].hi - e[ax].lo)
            ax = i;
    auto r = e;
    Int t = mid(r[ax]);
    e[ax].hi = r[ax].lo = t;
    return caps(m, e, depth - 1) && caps(m, r, depth - 1);
}
int limit = 0;
bool visit(Box b, int depth) {
    ++visits;
    if (limit && visits > limit)
        return false;
    BoxInfo info{b};
    if (!trim(info)) {
        return true;
    }
    for (auto &v : info.b) {
        constexpr Int step = U / 4096;
        v.lo = floor_ratio(v.lo, step) * step;
        v.hi = ceil_ratio(v.hi, step) * step;
    }
    b = info.b;
    maxdepth = std::max(maxdepth, depth);

    Model m(info, depth);
    Caps e = m.cap;
    if (!m.init_feasible())
        return true;
    if (use_transport) {
        if (!m.possible(e))
            return true;
        if (use_simplex && depth >= 9 && m.simplex_cover(e))
            return true;
        if (depth >= 9) {
            bool shorter = false;
            if (m.vertex_projection(e, b, shorter))
                return true;
            if (shorter) {
                ++contractions;
                return visit(b, depth);
            }
        }
    }
    if (caps(m, e, depth < 15 ? 1 : depth < 24 ? 4 : depth < 36 ? 6 : 10))
        return true;
    if (depth >= 70) {
        for (auto v : b)
            std::cout << v.lo << ":" << v.hi << " ";
        std::cout << "open" << std::endl;
        return false;
    }
    int ax = 0;
    for (int j = 1; j < 8; ++j)
        if (b[j].hi - b[j].lo > b[ax].hi - b[ax].lo)
            ax = j;
    if (m.direction >= 0)
        ax = m.direction;
    if (b[ax].hi - b[ax].lo <= U / 2048) {
        std::cout << "precision OPEN " << depth << std::endl;
        return false;
    }
    auto r = b;
    Int t = mid(r[ax]);
    b[ax].hi = r[ax].lo = t;
    return visit(b, depth + 1) && visit(r, depth + 1);
}
int main(int argc, char **argv) {
    if (argc < 6 || argc > 8)
        return 2;
    if (argc == 8) {
        auto mode = std::string_view(argv[7]);
        if (mode == "basic")
            use_pair_caps = use_sections = false;
        else if (mode == "pairs")
            use_sections = false;
        else if (mode == "interpolate")
            use_interpolation = true;
        else if (mode == "transport" || mode == "simplex" || mode == "curvature") {
            use_transport = use_interpolation = true;
            use_simplex = (mode == "simplex");
            use_curvature = (mode == "curvature");
        } else if (mode != "sections")
            return 2;
    }
    int input[4];
    for (int i = 0; i < 4; ++i)
        input[i] = std::atoi(argv[i + 1]);
    Model::init();
    int charts = 0;
    long long total = 0;
    int selected = std::atoi(argv[5]);
    bool ran = false;
    limit = argc > 6 ? std::atoi(argv[6]) : 0;
    for (int axes = 0; axes < 81; ++axes) {
        bool good = true;
        int a = axes;
        for (int i = 0; i < 4; ++i) {
            int k = a % 3;
            a /= 3;
            type.axis[i] = k;
            int bits = input[i];
            if (std::popcount((unsigned)bits) > 4)
                bits ^= 255;
            type.count[i] = std::popcount((unsigned)bits);
            int same = 7, v = std::countr_zero((unsigned)bits);
            for (int j = 0; j < 8; ++j)
                if (bits >> j & 1)
                    same &= ~(j ^ v);
            if (same && !(same >> k & 1))
                good = false;
            type.face[i] = same && type.count[i] == 4;
        }
        if (!good)
            continue;
        for (int g = 0; g < 16; ++g) {
            bool valid = true;
            for (int i = 0; i < 4; ++i) {
                int k = type.axis[i], bits = input[i];
                if (std::popcount((unsigned)bits) > 4)
                    bits ^= 255;
                if (!(g >> i & 1))
                    bits ^= 255;
                type.bits[i] = bits;
                int sign = 0;
                for (int v = 0; v < 8; ++v)
                    if (!(v >> k & 1) && (((bits >> v) ^ (bits >> (v ^ (1 << k)))) & 1))
                        sign = (bits >> v & 1) ? -1 : 1;
                assert(sign);
                type.orient[i] = sign;
                for (int j = 0; j < i; ++j)
                    if ((type.face[i] || type.count[i] == 3) && (type.face[j] || type.count[j] == 3) &&
                        type.axis[j] == k && type.orient[j] == sign)
                        valid = false;
            }
            for (int v = 0; v < 8; ++v) {
                int word = 0;
                for (int i = 0; i < 4; ++i)
                    word |= ((type.bits[i] >> v) & 1) << i;
                valid &= word != 0 && word != 15;
            }
            if (!valid)
                continue;
            if (selected >= 0 && charts++ != selected)
                continue;
            if (selected < 0)
                ++charts;
            ran = true;
            for (auto &t : lattice_cache)
                t.clear();
            for (auto &t : polygon_cache)
                t.clear();
            for (auto &t : pair_tables)
                t.clear();
            for (auto &t : slice_cache)
                t.clear();
            visits = capnodes = contractions = maxdepth = 0;
            Box b;
            b.fill({-2 * U, 2 * U});
            bool ok = visit(b, 0);
            total += visits;
            std::cout << "chart " << charts - 1 << " axes " << axes << " g " << g << " "
                      << (ok ? "proved" : "OPEN") << " nodes " << visits << " depth " << maxdepth;
            if (ok && use_transport) {
                assert((visits - contractions) % 2 == 1);
                std::cout << " leaves " << (visits - contractions + 1) / 2 << " contractions "
                          << contractions;
            }
            std::cout << std::endl;
            if (!ok)
                return 1;
        }
    }
    if (!ran)
        return 2;
    if (selected < 0)
        std::cout << "proved " << charts << " charts " << total << " nodes" << std::endl;
    return 0;
}
