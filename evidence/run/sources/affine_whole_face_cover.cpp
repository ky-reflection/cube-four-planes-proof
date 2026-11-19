// Four whole-face cuts with primary axes x, y, z, z. Outward sides on z are opposite.
#include "affine_box_linear.hpp"
#include <cassert>
#include <iostream>

using namespace cube_affine_box;
#include "affine_face_pair_moments.hpp"
using namespace cube_face_pairs;
constexpr int rep[]{1, 2, 4, 8, 3, 5, 9}, axis[]{0, 1, 2, 2}, orient[]{1, 1, 1, -1};
constexpr int pairid[][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
struct Pair {
    R value;
    std::array<R, 8> grad{};
    R ac[3]{};
};
struct Model {
    std::array<V, 4> n{};
    Caps cap, lambda;
    R den{};
    Pair pairs[6];
    std::array<Int, 128> lower{};
    R offset[8];
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
    Model(const Box &box) {
        V centers[4];
        for (int i = 0; i < 4; ++i) {
            R a = box[2 * i], b = box[2 * i + 1];
            int k = axis[i];
            n[i][k] = pt(orient[i] * U);
            n[i][(k + 1) % 3] = (a - b) / 2;
            n[i][(k + 2) % 3] = (a + b) / 2;
            cap[i] = {-std::min(G, (U - std::max(minabs(a), minabs(b)) + 1) / 2),
                      std::min(G + T, (U - std::max(minabs(a), minabs(b)) + 1) / 2)};
            centers[i][k] = n[i][k];
            centers[i][(k + 1) % 3] = (pt(mid(a)) - pt(mid(b))) / 2;
            centers[i][(k + 2) % 3] = (pt(mid(a)) + pt(mid(b))) / 2;
        }
        for (int h = 0; h < 6; ++h) {
            auto [i, j] = pairid[h];
            auto &p = pairs[h];
            V gr[2]{};
            int k = axis[i], l = axis[j];
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
                    return std::array<R, 4>{times(n[i][l], -orient[j]), -n[i][t],
                                            times(n[j][k], -orient[i]), -n[j][t]};
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
                gr[0][l] = times(pt(U / 4) - square(g1) / 12, -orient[j]);
                gr[0][t] = g1 / 6;
                gr[1][k] = times(pt(U / 4) - square(g2) / 12, -orient[i]);
                gr[1][t] = g2 / 6;
                p.ac[0] = a;
                p.ac[1] = c;
                p.ac[2] = z;
            }
            for (int w = 0; w < 2; ++w) {
                int v = w ? j : i, ax = axis[v];
                R x = gr[w][(ax + 1) % 3], y = gr[w][(ax + 2) % 3];
                p.grad[2 * v] = (x + y) / 2;
                p.grad[2 * v + 1] = (y - x) / 2;
            }
        }
        for (int i = 0; i < 4; ++i) {
            auto a = det(n[(i + 1) % 4], n[(i + 2) % 4], n[(i + 3) % 4]);
            lambda[i] = i % 2 ? -a : a;
        }
        for (int j = 0; j < 8; ++j)
            offset[j] = box[j] - pt(mid(box[j]));
        for (int m = 1; m < 128; ++m) {
            R v = pt(co[m][6] * U);
            std::array<R, 8> g{};
            for (int h = 0; h < 6; ++h)
                if (int c = co[m][h]) {
                    v = v + times(pairs[h].value, c) / 2;
                    for (int j = 0; j < 8; ++j)
                        g[j] = g[j] + times(pairs[h].grad[j], c) / 2;
                }
            for (int j = 0; j < 8; ++j)
                v = v + g[j] * offset[j];
            lower[m] = v.lo;
        }
    }
    bool init_feasible() {
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
    bool possible(Caps &e) const {
        for (int sweep = 0; sweep < 2; ++sweep) {
            R massT{0, T}, B{};
            for (int j = 0; j < 4; ++j)
                B = B + lambda[j] * e[j];
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
                if (axis[i] == axis[j]) {
                    R z = {0, bm.hi};
                    if (p.ac[0].lo > U / 1000)
                        z.hi = std::min(z.hi, (square(e[i] + e[j]) / p.ac[0]).hi);
                    corr[h] = z;
                } else if (p.ac[2].lo > U / 1000)
                    corr[h] = intersect(
                        bm, (p.ac[1] * square(e[i]) + times(e[i] * e[j], 2) + p.ac[0] * square(e[j])) /
                                p.ac[2]);
                else
                    corr[h] = bm;
                if (corr[h].lo > corr[h].hi)
                    return false;
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
long long visits = 0, capnodes = 0;
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
bool visit(Box b, int depth) {
    ++visits;
    maxdepth = std::max(maxdepth, depth);
    Model m(b);
    if (!m.init_feasible() || caps(m, m.cap, depth < 15 ? 1 : depth < 24 ? 4 : depth < 36 ? 6 : 10))
        return true;
    if (depth >= 70) {
        for (auto v : b)
            std::cout << v.lo << ":" << v.hi << " ";
        std::cout << "open" << std::endl;
        return false;
    }
    int ax = depth % 8;
    auto r = b;
    Int t = mid(r[ax]);
    b[ax].hi = r[ax].lo = t;
    return visit(b, depth + 1) && visit(r, depth + 1);
}
int main() {
    Model::init();
    Box b;
    b.fill({-U, U});
    bool ok = visit(b, 0);
    std::cout << (ok ? "proved" : "OPEN") << " nodes " << visits << " caps " << capnodes << " depth "
              << maxdepth << std::endl;
    return !ok;
}
