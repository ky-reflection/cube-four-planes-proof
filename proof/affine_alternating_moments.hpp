#ifndef CUBE_AFFINE_ALTERNATING_MOMENTS_HPP
#define CUBE_AFFINE_ALTERNATING_MOMENTS_HPP
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>
namespace cube_moments {
using Int = std::int64_t;
using Wide = __int128_t;
constexpr int grid = 1 << 15;
constexpr Int mass = Int{1} << 40;
constexpr Int infinity = std::numeric_limits<Int>::max();

struct Coordinate {
    int lo, hi;
};
struct Bound {
    Int lo, hi;
};

Int floor_ratio(Wide numerator, Int denominator) {
    if (denominator <= 0)
        throw std::logic_error("nonpositive denominator");
    return numerator < 0 ? -((-numerator + denominator - 1) / denominator) : numerator / denominator;
}

Int ceil_ratio(Wide numerator, Int denominator) { return -floor_ratio(-numerator, denominator); }

Bound share(Coordinate a, Coordinate b) {
    const int lo_den = a.lo + b.hi;
    const int hi_den = a.hi + b.lo;
    return {lo_den ? floor_ratio(Wide{mass} * a.lo, lo_den) : 0,
            hi_den ? ceil_ratio(Wide{mass} * a.hi, hi_den) : mass};
}

Bound linked_ratio(Coordinate a, Coordinate b, Coordinate e) {
    // Enclose (a+e)/(a+b), taking the same a in numerator and denominator.
    Bound result{infinity, 0};
    for (const int x : {a.lo, a.hi}) {
        result.lo = std::min(result.lo, x + b.hi ? floor_ratio(Wide{mass} * (x + e.lo), x + b.hi) : 0);
        result.hi = std::max(result.hi, x + b.lo ? ceil_ratio(Wide{mass} * (x + e.hi), x + b.lo) : infinity);
    }
    return result;
}

using Matrix = std::array<std::array<Coordinate, 4>, 4>;
constexpr std::array<std::array<int, 2>, 6> edges{{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
struct I {
    Int l = 0, h = 0;
};
I neg(I v) { return {-v.h, -v.l}; }
I add(I a, I b) { return {a.l + b.l, a.h + b.h}; }
I sub(I a, I b) { return add(a, neg(b)); }
I scalar(I a, Int b) { return b >= 0 ? I{a.l * b, a.h * b} : I{a.h * b, a.l * b}; }
I mul(I a, I b) {
    auto p = std::minmax({Wide(a.l) * b.l, Wide(a.h) * b.h, Wide(a.l) * b.h, Wide(a.h) * b.l});
    return {floor_ratio(p.first, mass), ceil_ratio(p.second, mass)};
}
I divi(I a, I b) {
    assert(b.l > 0);
    return {floor_ratio(Wide(a.l) * mass, a.l < 0 ? b.l : b.h),
            ceil_ratio(Wide(a.h) * mass, a.h < 0 ? b.h : b.l)};
}
I divide(I a, Int b) { return {floor_ratio(a.l, b), ceil_ratio(a.h, b)}; }
I sq(I a) {
    if (a.l >= 0)
        return mul(a, a);
    if (a.h <= 0)
        return mul(neg(a), neg(a));
    return {0, ceil_ratio(std::max(Wide(a.l) * a.l, Wide(a.h) * a.h), mass)};
}
Int absmax(I a) { return std::max(a.h, -a.l); }
Int absmin(I a) { return std::max({Int(0), a.l, -a.h}); }
I intersection(I a, I b) { return {std::max(a.l, b.l), std::min(a.h, b.h)}; }
I positive(I a) { return {std::max(Int(0), a.l), std::max(Int(0), a.h)}; }
constexpr Int Q = (79759 * mass + 999999) / 1000000, TMAX = (8452 * mass + 999999) / 1000000;
constexpr int reps[] = {1, 2, 4, 8, 3, 5, 9};
struct Support {
    std::array<int, 4> w{};
    int weight = 0;
    std::array<int, 7> masks{}, multipliers{};
    int n = 0;
};
static const auto supports = [] {
    std::vector<Support> s;
    for (int a = -1; a <= 1; ++a)
        for (int b = -1; b <= 1; ++b)
            for (int c = -1; c <= 1; ++c)
                for (int d = -1; d <= 1; ++d) {
                    Support r;
                    r.w = {a, b, c, d};
                    r.weight = a + b + c + d;
                    if (a * 27 + b * 9 + c * 3 + d <= 0)
                        continue;
                    std::array<int, 7> dots{};
                    int maxval = 0;
                    for (int z = 0; z < 7; ++z) {
                        for (int i = 0; i < 4; ++i)
                            dots[z] += r.w[i] * ((reps[z] >> i & 1) ? 1 : -1);
                        dots[z] = std::abs(dots[z]);
                        maxval = std::max(maxval, dots[z]);
                    }
                    for (int v = 1; v <= maxval; ++v) {
                        int m = 0;
                        for (int z = 0; z < 7; ++z)
                            if (dots[z] >= v)
                                m |= 1 << z;
                        if (r.n && r.masks[r.n - 1] == m)
                            ++r.multipliers[r.n - 1];
                        else {
                            r.masks[r.n] = m;
                            r.multipliers[r.n++] = 1;
                        }
                    }
                    s.push_back(r);
                }
    return s;
}();
struct PairRegion {
    bool usable = false;
    I slope[2][2], share[2], loss[2];
    I major, minor;
};
using Linear4 = std::array<I, 4>;
I linear_value(const Linear4 &a, const std::array<I, 4> &H) {
    I v{};
    for (int j = 0; j < 4; ++j)
        v = add(v, mul(a[j], H[j]));
    return v;
}
struct CubicBound {
    Linear4 lin, argument;
    I coefficient;
};
struct TripleOdd {
    Linear4 coeff{};
    std::vector<CubicBound> cubic;
    bool valid = false;
};
struct OddModel {
    std::array<TripleOdd, 4> triples;
    OddModel(const std::array<std::array<I, 4>, 4> &c, const I (&cof)[4]) {
        const I one{mass, mass};
        I n[4][3];
        auto sign = [](int i, int k) { return i == 0 || i == k + 1 ? 1 : -1; };
        for (int i = 0; i < 4; ++i)
            for (int k = 0; k < 3; ++k) {
                I component{};
                for (int j = 0; j < 4; ++j)
                    if (j != i && j != (i ^ (k + 1)))
                        component = add(component, c[i][j]);
                n[i][k] = scalar(divide(component, 2), sign(i, k));
            }
        for (int o = 0; o < 4; ++o) {
            auto &tri = triples[o];
            I det = divide(cof[o], 4);
            if (det.l < mass / 80)
                continue;
            bool ok = true;
            Linear4 point[3]{};
            for (int i = 0; i < 4; ++i)
                if (i != o) {
                    int j = 0;
                    while (j == o || j == i)
                        ++j;
                    int k = 6 - o - i - j;
                    int dsign = 0;
                    I cross[3];
                    for (int x = 0; x < 3; ++x) {
                        int y = (x + 1) % 3, z = (x + 2) % 3;
                        cross[x] = sub(mul(n[j][y], n[k][z]), mul(n[j][z], n[k][y]));
                        dsign += sign(i, x) * (sign(j, y) * sign(k, z) - sign(j, z) * sign(k, y));
                    }
                    for (int x = 0; x < 3; ++x)
                        point[x][i] = scalar(divi(cross[x], det), dsign > 0 ? 1 : -1);
                }
            auto face = [&](int i, int j) {
                Linear4 p = point[(i ^ j) - 1];
                for (auto &t : p)
                    t = scalar(t, (i == 0 || j == 0) ? 1 : -1);
                return p;
            };
            for (int k = 0; k < 4; ++k)
                if (k != o) {
                    int i = 0;
                    while (i == o || i == k)
                        ++i;
                    int j = 6 - o - i - k;
                    I g = c[k][o];
                    auto d = mul(add(g, c[k][i]), add(g, c[k][j]));
                    if (d.l < mass / 512) {
                        ok = false;
                        break;
                    }
                    auto inv = divi(one, d);
                    auto w = face(i, j);
                    tri.coeff[k] = add(tri.coeff[k], mul(scalar(g, 4), inv));
                    for (int t = 0; t < 4; ++t)
                        tri.coeff[t] = add(tri.coeff[t], mul(scalar(mul(sq(g), inv), 2), w[t]));
                    Linear4 arg{};
                    arg[k] = one;
                    tri.cubic.push_back({w, arg, scalar(inv, 2)});
                }
            for (int i = 0; i < 4 && ok; ++i)
                for (int j = i + 1; j < 4 && ok; ++j)
                    if (i != o && j != o) {
                        int k = 6 - o - i - j;
                        I l = add(c[i][k], c[i][o]), r = add(c[j][k], c[j][o]);
                        if (std::min(l.l, r.l) < mass / 32) {
                            ok = false;
                            break;
                        }
                        l = divi(one, l);
                        r = divi(one, r);
                        I x = sub(add(mul(c[i][k], l), mul(c[j][o], r)), one);
                        I df = add(mul(add(c[i][k], c[i][j]), l), mul(add(c[j][o], c[j][i]), r));
                        I dk = add(mul(add(c[j][k], c[j][i]), r), mul(add(c[i][o], c[i][j]), l));
                        if (std::min(df.l, dk.l) < mass / 32) {
                            ok = false;
                            break;
                        }
                        df = divi(one, df);
                        dk = divi(one, dk);
                        I product = mul(df, dk), a = divide(add(df, dk), 2);
                        I xabs{absmin(x), absmax(x)};
                        I first = sub(scalar(mul(a, x), 4), scalar(mul(product, mul(x, xabs)), 6));
                        Linear4 arg{};
                        arg[i] = l;
                        arg[j] = neg(r);
                        auto fi = face(i, k), fj = face(j, k);
                        Linear4 twice_b;
                        for (int t = 0; t < 4; ++t) {
                            twice_b[t] = add(mul(fi[t], df), mul(fj[t], dk));
                            tri.coeff[t] = add(tri.coeff[t], sub(mul(first, arg[t]), mul(twice_b[t], sq(x))));
                        }
                        tri.cubic.push_back({twice_b, arg, {-mass, -mass}});
                        // The odd part of |x+y|^3 beyond its linear term has absolute value <= 2 |y|^3.
                        tri.cubic.push_back({arg, arg, {-2 * product.h, 2 * product.h}});
                    }
            if (ok) {
                for (auto &v : tri.coeff)
                    v = {std::max(Int(0), v.l), std::max(Int(0), v.h)};
                tri.valid = true;
            }
        }
    }
    using Caps = std::array<I, 4>;
    bool contract(Caps &e, const I (&H)[4], const I (&ratio)[4], I t, const std::array<Int, 128> &slack,
                  std::array<Linear4, 7> &L, std::array<I, 7> &offset) const {
        std::array<I, 4> hh;
        for (int i = 0; i < 4; ++i)
            hh[i] = H[i];

        for (int i = 0; i < 4; ++i) {
            auto &f = triples[i];
            if (!f.valid) {
                for (int j = 0; j < 4; ++j)
                    if (j != i)
                        L[i][j] = {0, 2 * mass};
                offset[i] = neg(t);
                continue;
            }
            for (int j = 0; j < 4; ++j)
                L[i][j] = divi(divide(f.coeff[j], 12), ratio[j]);
            I remainder{};
            for (auto &z : f.cubic) {
                I a = linear_value(z.lin, hh), b = linear_value(z.argument, hh);
                remainder = add(remainder, mul(z.coefficient, mul(a, sq(b))));
            }
            offset[i] = sub(divide(remainder, 12), t);
        }
        for (int k = 4; k < 7; ++k) {
            I r{};
            for (int i = 0; i < 4; ++i) {
                int s = (reps[k] >> i & 1) ? -1 : 1;
                // sum of signs vanishes, so the common -t cancels exactly.
                auto rem = add(offset[i], t);
                if (triples[i].valid) {
                } else
                    rem = {};
                r = add(r, scalar(rem, s));
                L[k][i] = add(L[k][i], I{s * mass, s * mass});
                for (int j = 0; j < 4; ++j)
                    L[k][j] = add(L[k][j], scalar(L[i][j], s));
            }
            offset[k] = divide(r, 2);
            for (auto &v : L[k])
                v = divide(v, 2);
        }
        Int minimum[7]{};
        for (int k = 0; k < 7; ++k) {
            I val = add(linear_value(L[k], e), offset[k]);
            minimum[k] = absmin(val);
            Int room = slack[1 << k];
            if (minimum[k] > room)
                return false;
            for (int j = 0; j < 4; ++j) {
                auto coeff = L[k][j];
                if (absmin(coeff) < mass / 16384)
                    continue;
                I term = mul(coeff, e[j]), other{val.l - term.l, val.h - term.h},
                  rhs = {-room - other.h, room - other.l};
                I bound = coeff.l > 0 ? divi(rhs, coeff) : divi(neg(rhs), neg(coeff));
                auto v = intersection(e[j], bound);
                if (v.l > v.h)
                    return false;
                // Keep the bounds used to form val consistent while contracting each coordinate.
                val = add(other, mul(coeff, v));
                e[j] = v;
            }
        }
        Int minsum[128] = {};
        for (int m = 1; m < 128; ++m) {
            minsum[m] = minsum[m & (m - 1)] + minimum[std::countr_zero(unsigned(m))];
            if (minsum[m] > slack[m])
                return false;
        }
        return true;
    }
};

constexpr Int micro = 1000000;
struct AngularErrors {
    std::array<std::array<Int, 2>, 6> rates{};
    AngularErrors(const Matrix &raw, Int volume_target = Q) {
        const I one{mass, mass};
        for (int h = 0; h < 6; ++h) {
            auto [i, j] = edges[h];
            Int bmpos = (micro * micro + 3 * 345462 - 1) / (3 * 345462),
                bmneg = (micro * micro + 3 * 379464 - 1) / (3 * 379464);
            rates[h] = {bmneg, bmpos};
            I bias[2], s[2], a[2]{};
            bool ok = true;
            int shared = (i ^ j) - 1;
            for (int l = 0; l < 2; ++l) {
                int r = l ? j : i, v = l ? i : j, u = r ^ (((shared + 1) % 3) + 1),
                    w = r ^ (((shared + 2) % 3) + 1);
                auto c = [&](int k) { return I{raw[r][k].lo * (mass / grid), raw[r][k].hi * (mass / grid)}; };
                I g = c(v), b = c(u), d = c(w), sum = add(b, d), den = mul(add(g, b), add(g, d));
                if (sum.l < mass / 64 || den.l < mass / 1024) {
                    ok = false;
                    break;
                }
                auto share1 = share(raw[r][v], raw[r][u]), share2 = share(raw[r][v], raw[r][w]);
                I p{share1.lo, share1.hi}, q{share2.lo, share2.hi};
                I tau = mul(p, q);
                Int H = std::min({g.h, b.h, d.h, 3 * (7 * volume_target - mass / 2 + TMAX)});
                s[l] = {mass - tau.h - divi(divide(sq({H, H}), 3), den).h, mass - tau.l};
                if (s[l].l < mass / 64) {
                    ok = false;
                    break;
                }
                bias[l] = {mass, 0};
                for (Int x : {p.l, p.h})
                    for (Int y : {q.l, q.h}) {
                        auto xy = mul({x, x}, {y, y});
                        bias[l].l = std::min(bias[l].l, x + y - 2 * xy.h);
                        bias[l].h = std::max(bias[l].h, x + y - 2 * xy.l);
                    }
                for (int z = 0; z < 2; ++z) {
                    // |n_{axis}|/|n_shared|. The corresponding numerator omits the score indexed by that
                    // axis.
                    int kept = z ? u : w, dropped = z ? w : u;
                    auto q = linked_ratio(raw[r][kept], raw[r][dropped], raw[r][v]);
                    a[z] = add(a[z], {q.lo, q.hi});
                }
            }
            if (!ok)
                continue;
            I A{std::max(a[0].l, a[1].l), std::max(a[0].h, a[1].h)};
            Int Bmin = std::min(a[0].l, a[1].l), Bmax = std::min(a[0].h, a[1].h);
            if (A.l <= 0 || A.h > 64 * mass)
                continue;
            Int hcap = 7 * volume_target - mass / 2 + TMAX;
            Int maxdiff = divi({hcap, hcap}, s[0]).h + divi({hcap, hcap}, s[1]).h;
            Int factor = 2 * mass / 3;
            if (Bmin > 0)
                factor = std::max(factor, mass - ceil_ratio(Wide(maxdiff) * mass, 3 * Bmin));
            if (2 * maxdiff + Bmax <= A.l)
                factor = mass;
            for (int up = 0; up <= 1; ++up) {
                Int k = up ? floor_ratio(Wide(factor) * mass, A.h) : ceil_ratio(Wide(mass) * mass, A.l);
                Int bi = up ? bias[0].h : bias[0].l, bj = up ? bias[1].h : bias[1].l;
                Int ai = up ? (bi - k) : (k - bi), aj = up ? (bj - k) : (k - bj),
                    den = std::numeric_limits<Int>::max();
                // det(C*11' - Q) after division by the two squared scale factors.
                for (Int x : {s[0].l, s[0].h})
                    for (Int y : {s[1].l, s[1].h}) {
                        I xx = sq({x, x}), yy = sq({y, y}), xy = mul({x, x}, {y, y});
                        I v = add(add(mul({-aj, -aj}, xx), mul({-ai, -ai}, yy)), mul({2 * k, 2 * k}, xy));
                        den = std::min(den, v.l);
                    }
                if (den <= mass / 100000)
                    continue;
                Int n = mul({k, k}, {bi + bj, bi + bj}).h - mul({bi, bi}, {bj, bj}).l;
                Int C = std::max({Int(0), divi({n, n}, {den, den}).h, divi({ai, ai}, sq(s[0])).h,
                                  divi({aj, aj}, sq(s[1])).h});
                rates[h][up] = std::min(rates[h][up], ceil_ratio(Wide(C) * micro, mass));
            }
        }
    }
};
struct Geometry {
    std::optional<OddModel> odd;
    I ratios[4], tree[4], max_h[4], den{};
    PairRegion pair[6];
    std::array<std::array<I, 4>, 4> c{};
    Geometry(const Matrix &raw) {
        const I one{mass, mass};
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                c[i][j] = {raw[i][j].lo * (mass / grid), raw[i][j].hi * (mass / grid)};
        for (int i = 0; i < 4; ++i) {
            // Tree cofactor: one of the 16 labeled directed trees directed to the root i.
            for (int a = 0; a < 4; ++a)
                if (a != ((i + 1) % 4))
                    for (int b = 0; b < 4; ++b)
                        if (b != ((i + 2) % 4))
                            for (int d = 0; d < 4; ++d)
                                if (d != ((i + 3) % 4)) {
                                    int f[4];
                                    f[i] = i;
                                    f[(i + 1) % 4] = a;
                                    f[(i + 2) % 4] = b;
                                    f[(i + 3) % 4] = d;
                                    bool ok = true;
                                    for (int j = 0; j < 4; ++j)
                                        ok &= f[f[f[j]]] == i;
                                    if (ok) {
                                        I v = one;
                                        for (int j = 0; j < 4; ++j)
                                            if (j != i)
                                                v = mul(v, c[j][f[j]]);
                                        tree[i] = add(tree[i], v);
                                    }
                                }
        }
        odd.emplace(c, tree);
        // Ratios relative to a common tree scale keep all intermediate fixed-point quantities small.
        Int scale = mass;
        for (auto p : tree)
            scale = std::max(scale, p.h);
        I ref{scale, scale};
        for (auto &u : tree)
            u = divi(u, ref);
        den = scalar(ref, 12);
        for (auto v : tree)
            den = mul(den, v);
        for (unsigned e = 0; e < 6; ++e) {
            auto [i, j] = edges[e];
            int x = (i ^ j) - 1, y = (x + 1) % 3, z = (x + 2) % 3;
            auto &p = pair[e];
            p.usable = true;
            I axes[2]{};
            for (int t = 0; t < 2; ++t) {
                int r = t ? j : i, s = t ? i : j;
                auto component = [&](int axis) {
                    I v{};
                    for (int h = 0; h < 4; ++h)
                        if (h != r && h != (r ^ (axis + 1)))
                            v = add(v, c[r][h]);
                    return v;
                };
                I common = component(x), nx = component(y), ny = component(z);
                if (common.l < mass / 32 || nx.l < mass / 32 || ny.l < mass / 32) {
                    p.usable = false;
                    break;
                }
                p.share[t] = divi(one, common);
                axes[0] = add(axes[0], divi(nx, common));
                axes[1] = add(axes[1], divi(ny, common));
                p.loss[t] = divi(c[r][s], mul(common, mul(nx, ny)));
            }
            if (p.usable) {
                p.major = {std::max(axes[0].l, axes[1].l), std::max(axes[0].h, axes[1].h)};
                p.minor = {std::min(axes[0].l, axes[1].l), std::min(axes[0].h, axes[1].h)};
            }
        }
    }
};
} // namespace cube_moments
#endif
