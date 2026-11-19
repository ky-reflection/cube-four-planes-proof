// All independent shifts in the angular neighborhood certified by the centered merit.
#include "affine_alternating_moments.hpp"
#include <iostream>
using namespace cube_moments;
using CapBox = std::array<I, 4>;
long long ntest = 0;
int radial_depth = 0;
#include "affine_angular_merit_weights.hpp"
std::array<I, 7> weights_at(const std::array<I, 8> &p) {
    I inv = divi({mass, mass}, {1554944255987, 1554944255988});
    I A = divide(inv, 6), B = divide(sub({mass, mass}, inv), 6);
    std::array<I, 7> out{};
    constexpr int masks[]{3, 5, 9, 12, 10, 6, 7, 11, 13, 8, 4, 2};
    for (int z = 0; z < 12; ++z) {
        int m = masks[z];
        m = (m & 5) | ((m & 2) << 2) | ((m & 8) >> 2);
        I w = z < 6 ? A : B;
        for (int j = 0; j < 8; ++j)
            w = add(w, divide(scalar(p[j], weight_slopes[z][j]), 8192));
        assert(w.l > 0);
        for (int r = 0; r < 7; ++r)
            if (m == reps[r] || m == (15 ^ reps[r]))
                out[r] = add(out[r], divide(w, 2));
    }
    return out;
}

struct OffsetTransfer {
    Geometry dom;
    std::array<I, 7> weight{};
    OffsetTransfer(const Matrix &c, std::array<I, 7> w) : dom(c), weight(w) {}
    bool ok(CapBox u, I S) const {
        CapBox hh;
        I B{};
        for (int i = 0; i < 4; ++i) {
            hh[i] = divi(u[i], dom.ratios[i]);
            B = add(B, mul(dom.tree[i], hh[i]));
        }
        if (B.h < 0)
            return true;
        I ub = positive(B), t;
        if (dom.den.l <= mass / 100000)
            return false;
        t = divi(mul(mul(sq(ub), ub), sq(S)), dom.den);
        std::array<I, 7> O{}, pure{};
        for (int i = 0; i < 4; ++i) {
            auto &f = dom.odd->triples[i];
            if (!f.valid)
                return false;
            pure[i] = linear_value(f.coeff, hh);
            I cubic{};
            for (auto &term : f.cubic) {
                auto A = linear_value(term.lin, hh), B = linear_value(term.argument, hh);
                cubic = add(cubic, mul(term.coefficient, mul(A, sq(B))));
            }
            pure[i] = divide(add(pure[i], mul(sq(S), cubic)), 12);
            O[i] = sub(pure[i], t);
        }
        for (int r = 4; r < 7; ++r) {
            I v{};
            for (int j = 0; j < 4; ++j)
                v = add(v, scalar(add(u[j], pure[j]), (reps[r] >> j & 1) ? -1 : 1));
            O[r] = divide(v, 2);
        }
        I score{};
        I ct{};
        for (int r = 0; r < 7; ++r) {
            score = add(score, mul(weight[r], {absmin(O[r]), absmax(O[r])}));
            ct = add(ct, scalar(weight[r], r < 4 ? 1 : -1));
        }
        Int best = score.l;
        for (auto &test : supports) {
            int free = 0;
            I lhs = scalar(t, test.weight);
            for (int j = 0; j < 4; ++j) {
                lhs = sub(lhs, scalar(u[j], 2 * test.w[j]));
                free += test.w[j] * (j == 0 ? 1 : -1);
            }
            lhs = sub(lhs, scalar(O[0], free));
            Int loc = absmin(lhs);
            int coef[7]{};
            for (int r = 1; r < 7; ++r) {
                for (int j = 0; j < 4; ++j)
                    coef[r] += test.w[j] * (reps[r] >> j & 1 ? 1 : -1);
                coef[r] = std::abs(coef[r]);
            }
            for (int p = 1; p < 7; ++p)
                if (coef[p]) {
                    Int multiplier = weight[p].l / coef[p];
                    Int l = mul({multiplier, multiplier}, {loc, loc}).l;
                    for (int r = 1; r < 7; ++r) {
                        Int a = weight[r].l - coef[r] * multiplier;
                        auto z = a < 0 ? absmax(O[r]) : absmin(O[r]);
                        l += mul({a, a}, {z, z}).l;
                    }
                    best = std::max(best, l);
                }
        }
        score.l = best;
        score = sub(score, mul(ct, t));
        for (int h = 0; h < 6; ++h) {
            auto [i, j] = edges[h];
            I w{};
            for (int r = 0; r < 7; ++r) {
                int c = r < 4 ? -(r != i && r != j) : int(((reps[r] >> i) ^ (reps[r] >> j)) & 1);
                w = add(w, scalar(weight[r], c));
            }
            I upper = mul(sq(add(u[i], u[j])),
                          {0, ceil_ratio(Wide(ceil_ratio(Wide(micro) * micro, 3 * 345462)) * mass, micro)});
            I lower =
                neg(mul(sq(sub(u[i], u[j])),
                        {0, ceil_ratio(Wide(ceil_ratio(Wide(micro) * micro, 3 * 379464)) * mass, micro)}));
            I d{lower.l, upper.h};
            auto &p = dom.pair[h];
            if (p.usable) {
                I g = sub(mul(hh[i], p.share[0]), mul(hh[j], p.share[1]));
                I f{2 * mass / 3, mass};
                if (p.minor.l > 0)
                    f.l = std::max(f.l, mass - ceil_ratio(Wide(mul({absmax(g), absmax(g)}, S).h) * mass,
                                                          3 * p.minor.l));
                if (2 * mul({absmax(g), absmax(g)}, S).h + p.minor.h <= p.major.l)
                    f.l = mass;
                auto v = sub(add(mul(p.loss[0], sq(hh[i])), mul(p.loss[1], sq(hh[j]))),
                             mul(f, divi(sq(g), p.major)));
                d = intersection(d, v);
            }
            score = sub(score, divide(mul(w, mul(S, d)), 2));
        }
        return score.l > 0;
    }
    bool walk(CapBox u, I s, int depth, int limit) const {
        ++ntest;
        radial_depth = std::max(radial_depth, depth);
        if (ok(u, s))
            return true;
        if (depth == limit)
            return false;
        int k = 0;
        for (int j = 1; j < 4; ++j)
            if (u[j].h - u[j].l > u[k].h - u[k].l)
                k = j;
        if (16 * (s.h - s.l) > u[k].h - u[k].l) {
            auto r = s;
            Int m = (s.l + s.h) / 2;
            s.h = r.l = m;
            return walk(u, s, depth + 1, limit) && walk(u, r, depth + 1, limit);
        }
        auto r = u;
        Int m = (u[k].l + u[k].h) / 2;
        u[k].h = r[k].l = m;
        return walk(u, s, depth + 1, limit) && walk(r, s, depth + 1, limit);
    }
    bool all(int limit) const {
        for (int i = 0; i < 4; ++i)
            for (int sign : {-1, 1}) {
                CapBox u;
                u.fill({-mass, mass});
                u[i] = {sign * mass, sign * mass};
                if (!walk(u, {0, sign < 0 ? 7 * Q - mass / 2 : 7 * Q - mass / 2 + TMAX}, 0, limit))
                    return false;
            }
        return true;
    }
};
using Angular = std::array<I, 8>;
bool localtest(const Angular &angles, int limit) {
    auto w = weights_at(angles);
    I one{mass, mass}, root{1554944255987, 1554944255988};
    I n[4][3];
    n[0][0] = add(one, angles[0]);
    n[0][1] = add(one, angles[1]);
    n[0][2] = sub(sub(one, angles[0]), angles[1]);
    n[1][0] = root;
    n[1][1] = neg(add(one, angles[6]));
    n[1][2] = neg(add(one, angles[7]));
    n[2][0] = neg(add(one, angles[4]));
    n[2][1] = root;
    n[2][2] = neg(add(one, angles[5]));
    n[3][0] = neg(add(one, angles[2]));
    n[3][1] = neg(add(one, angles[3]));
    n[3][2] = root;
    I scores[4][4]{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (i != j) {
                for (int k = 0; k < 3; ++k)
                    scores[i][j] = add(scores[i][j], scalar(n[i][k], j == 0 || j == 1 + k ? -1 : 1));
            }
    Matrix z{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (i != j) {
                if (!i) {
                    scores[i][j] = sub({3 * mass, 3 * mass}, scalar(n[i][j - 1], 2));
                }
                auto x = divide(scores[i][j], 2);
                z[i][j] = {int(floor_ratio(Wide(x.l) * grid, mass)),
                           int(ceil_ratio(Wide(x.h) * grid, mass))};
            }
    OffsetTransfer m(z, w);
    for (int i = 0; i < 4; ++i) {
        I c[3];
        int v = 0;
        for (int j = 0; j < 4; ++j)
            if (i != j)
                c[v++] = m.dom.c[i][j];
        auto slope = [&](bool lo) {
            I a[3];
            for (int j = 0; j < 3; ++j) {
                Int x = lo ? c[j].h : c[j].l;
                a[j] = {x, x};
            }
            auto aa = a[0], bb = a[1], cc = a[2];
            I num = add(add(mul(aa, bb), mul(aa, cc)), mul(bb, cc));
            auto den = mul(add(aa, bb), mul(add(aa, cc), add(bb, cc)));
            if (lo) {
                Int H = 3 * (7 * Q - mass / 2 + TMAX);
                num = sub(num, divide(sq({H, H}), 3));
            }
            return divi(num, den);
        };
        m.dom.ratios[i] = {slope(true).l, slope(false).h};
    }
    return m.all(limit);
}
int boxes = 0, leaves = 0, maxangular = 0, denominator;
bool certify(Angular b, int maximum, int depth = 0) {
    ++boxes;
    maxangular = std::max(maxangular, depth);
    if (localtest(b, 13)) {
        ++leaves;
        return true;
    }
    if (depth >= maximum) {
        if (localtest(b, 32)) {
            ++leaves;
            return true;
        }
        return false;
    }
    int k = depth % 8;
    auto r = b;
    Int mid = (b[k].l + b[k].h) / 2;
    r[k].l = b[k].h = mid;

    return certify(b, maximum, depth + 1) && certify(r, maximum, depth + 1);
}
int main(int argc, char **argv) {
    int den = argc > 1 ? std::atoi(argv[1]) : 35, limit = argc > 2 ? std::atoi(argv[2]) : 20;
    if (argc > 3 || den < 35 || limit < 0 || limit > 30)
        return 2;
    Angular p;
    p.fill(den ? I{-mass / den - 1, mass / den + 1} : I{});
    bool yes = certify(p, limit);
    std::cout << den << " " << (yes ? "proved" : "OPEN") << " angular " << boxes << " leaves " << leaves
              << " angular_depth " << maxangular << " radial " << ntest << " radial_depth " << radial_depth
              << std::endl;
    return !yes;
}
