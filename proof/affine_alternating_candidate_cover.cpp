// Affine cover for alternating star directions meeting the candidate neighborhood.
#include "affine_alternating_moments.hpp"
#include "affine_box_linear.hpp"
#include "affine_face_integral.hpp"
#include "affine_section_columns.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
using namespace cube_moments;
enum class ProjectionMode { plain, project, columns, flow, inverse };
ProjectionMode mode = ProjectionMode::plain;
bool use_columns() {
    return mode == ProjectionMode::columns || mode == ProjectionMode::flow || mode == ProjectionMode::inverse;
}
struct Box {
    std::array<Coordinate, 8> variables;
    int depth;
};
using Pivots = std::array<int, 4>;
using Pattern = std::array<int, 3>;
using RowBounds = std::array<Int, 27>;

constexpr std::array<std::array<int, 4>, 3> splits{{{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}}};

#include "centered_chart_terms.hpp"
std::array<Int, 128> thresholds;
const Int volume_target = cube_affine_box::candidate_upper();
#include "affine_alternating_supports.hpp"
#include "affine_angular_interpolation.hpp"
using QForm = std::array<std::array<I, 4>, 4>;
QForm outer_form(const Linear4 &a, const Linear4 &b, I scale) {
    QForm f{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j <= i; ++j)
            f[i][j] = mul(scale, divide(add(mul(a[i], b[j]), mul(a[j], b[i])), 2));
    return f;
}
void include_form(QForm &dest, const QForm &a, int n) {
    if (!n)
        return;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j <= i; ++j)
            dest[i][j] = add(dest[i][j], scalar(a[i][j], n));
}
struct QuadraticOffsets {
    QForm absoluteOdd[7]{}, transfer[2][6][2]{}, tet{};
    bool haveT = false;
    I t{};
    static Int cmax(I a) { return std::max(-a.l, a.h); }
    QuadraticOffsets(const std::array<I, 4> &e, const I (&H)[4], const I (&ratio)[4], const I (&tree)[4],
                     I den, I tetra, const OddModel &odd, const AngularErrors &errors,
                     const PairRegion (&pair)[6])
        : t(tetra) {
        I one{mass, mass}, B{};
        Linear4 g{};
        std::array<I, 4> hh;
        for (int j = 0; j < 4; ++j) {
            g[j] = divi(tree[j], ratio[j]);
            hh[j] = H[j];
            B = add(B, mul(tree[j], H[j]));
        }
        B = positive(B);
        I tOdd{};
        if (den.l > mass / 512) {
            auto w = divi(B, den);
            haveT = true;
            tet = outer_form(g, g, w);
            tOdd = divi(sq(B), den);
        }
        std::array<Linear4, 7> L{};
        bool valid[7]{};
        for (int i = 0; i < 4; ++i) {
            auto &f = odd.triples[i];
            valid[i] = f.valid;
            if (!valid[i]) {
                for (int j = 0; j < 4; ++j)
                    if (j != i)
                        L[i][j] = {0, 2 * mass};
                continue;
            }
            for (int j = 0; j < 4; ++j)
                L[i][j] = divi(divide(f.coeff[j], 12), ratio[j]);
            for (auto &z : f.cubic) {
                I ua = linear_value(z.argument, hh);
                I v = divide(mul(z.coefficient, sq(ua)), 12);
                for (int j = 0; j < 4; ++j)
                    L[i][j] = add(L[i][j], divi(mul(v, z.lin[j]), ratio[j]));
            }
        }
        for (int r = 4; r < 7; ++r) {
            valid[r] = true;
            for (int i = 0; i < 4; ++i) {
                int s = (reps[r] >> i & 1) ? -1 : 1;
                L[r][i] = add(L[r][i], {s * mass, s * mass});
                for (int j = 0; j < 4; ++j)
                    L[r][j] = add(L[r][j], scalar(L[i][j], s));
            }
            for (auto &x : L[r])
                x = divide(x, 2);
        }
        for (int r = 0; r < 7; ++r) {
            if (r < 4) {
                if (!haveT)
                    continue;
                for (int j = 0; j < 4; ++j)
                    L[r][j] = sub(L[r][j], mul(tOdd, g[j]));
            }
            Linear4 proposed{};
            for (int i = 0; i < 4; ++i) {
                Int b = (L[r][i].l + L[r][i].h) / 2;
                proposed[i] = {b, b};
            }
            Int maxval = std::max(mass / 64, absmax(linear_value(proposed, e)));
            I scale = divi(one, {maxval, maxval});
            // |beta.e|<=1. Whichever sign O takes, |O|>=(beta.e)*O.
            absoluteOdd[r] = outer_form(proposed, L[r], scale);
        }
        for (int h = 0; h < 6; ++h) {
            auto [i, j] = edges[h];
            auto &p = pair[h];
            for (int pos = 0; pos < 2; ++pos) {
                int sg = pos ? 1 : -1;
                Linear4 d{};
                d[i] = one;
                d[j] = {sg * mass, sg * mass};
                Int k = ceil_ratio(Wide(errors.rates[h][pos]) * mass, micro);
                transfer[0][h][pos] = outer_form(d, d, {k, k});
                transfer[1][h][pos] = transfer[0][h][pos];
            }
            if (!p.usable || p.major.l < mass / 24)
                continue;
            Linear4 v{};
            v[i] = divi(p.share[0], ratio[i]);
            v[j] = neg(divi(p.share[1], ratio[j]));
            I delta = linear_value(v, e);
            I factor{2 * mass / 3, mass};
            if (p.minor.l > 0)
                factor.l = std::max(factor.l, mass - ceil_ratio(Wide(absmax(delta)) * mass, 3 * p.minor.l));
            if (2 * absmax(delta) + p.minor.h <= p.major.l)
                factor.l = mass;
            QForm actual{};
            actual[i][i] = divi(p.loss[0], sq(ratio[i]));
            actual[j][j] = divi(p.loss[1], sq(ratio[j]));
            auto convex = outer_form(v, v, divi(factor, p.major));
            include_form(actual, convex, -1);
            transfer[1][h][1] = actual;
            for (auto &row : actual)
                for (auto &x : row)
                    x = neg(x);
            transfer[1][h][0] = actual;
        }
    }
    // Compare the minimum quadratic on the four-dimensional cap box with -gap.
    static bool enough(const QForm &f, const std::array<I, 4> &e, Int gap) {
        if (gap <= 0)
            return false;
        Int R2 = 0;
        for (auto v : e)
            R2 += sq(v).h;
        if (!R2)
            return true;
        Int shift = std::max(Int(0), floor_ratio(Wide(gap - 1) * mass, R2));
        for (int s = 0; s < 8; ++s) {
            I a[4][4];
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j <= i; ++j) {
                    Int x = ((s >> i) ^ (s >> j)) & 1 ? -f[i][j].h : f[i][j].l;
                    if (i == j)
                        x += shift;
                    a[j][i] = a[i][j] = {x, x};
                }
            // Exact Schur complements. Positive off-diagonal entries may be removed on this orthant.
            for (int i = 0; i < 4; ++i) {
                bool off = false;
                for (int j = i + 1; j < 4; ++j)
                    off |= a[j][i].l < 0;
                if (a[i][i].l < 0)
                    return false;
                if (!off)
                    continue;
                if (a[i][i].l < mass / 10000)
                    return false;
                for (int j = i + 1; j < 4; ++j)
                    for (int k = i + 1; k <= j; ++k)
                        a[j][k] = sub(a[j][k], divi(mul(a[j][i], a[k][i]), a[i][i]));
            }
        }
        return true;
    }
    bool excludes(const std::array<I, 4> &e, const std::array<Int, 128> &L) const {
        std::array<std::pair<Int, int>, 126> order;
        for (int m = 1; m < 127; ++m)
            order[m - 1] = {floor_ratio(L[m], 6) - 2 * volume_target * std::popcount(unsigned(m)), m};
        constexpr int maxw = 8;
        std::partial_sort(order.begin(), order.begin() + maxw, order.end(), std::greater<>());
        for (int n = 0; n < maxw; ++n) {
            auto [gap, m] = order[n];
            if (gap <= 0)
                break;
            QForm base{};
            int cl = 0, coeff[6]{};
            for (int r = 0; r < 7; ++r)
                if (m >> r & 1) {
                    include_form(base, absoluteOdd[r], 2);
                    cl += r < 4 ? 1 : -1;
                    for (int h = 0; h < 6; ++h) {
                        auto [i, j] = edges[h];
                        coeff[h] += r < 4 ? -(r != i && r != j) : int(((reps[r] >> i) ^ (reps[r] >> j)) & 1);
                    }
                }
            gap *= 2;
            if (cl > 0) {
                if (haveT)
                    include_form(base, tet, -2 * cl);
                else
                    gap -= 2 * cl * t.h;
            }
            for (int strategy = 0; strategy < 2; ++strategy) {
                auto f = base;
                for (int h = 0; h < 6; ++h)
                    if (coeff[h])
                        include_form(f, transfer[strategy][h][coeff[h] > 0], -std::abs(coeff[h]));
                if (enough(f, e, gap))
                    return true;
            }
        }
        return false;
    }
};
struct MomentDomain : Geometry {
    std::optional<AngularAtoms> dual;
    std::optional<AngularErrors> errors;
    const DirectionPairs::Interpolation *lattice = nullptr;
    std::array<Int, 128> lower{};
    MomentDomain(const Matrix &raw) : Geometry(raw) {
        const I one{mass, mass};
        for (int i = 0; i < 4; ++i) {
            int p = -1;
            for (int j = 0; j < 4; ++j)
                if (j != i && c[i][j].l == mass)
                    p = j;
            int a = -1, b = -1;
            for (int j = 0; j < 4; ++j)
                if (j != p && j != i) {
                    if (a < 0)
                        a = j;
                    else
                        b = j;
                }
            I x = c[i][a], y = c[i][b];
            Int maxH = std::min({x.h, y.h, 3 * (TMAX + 7 * volume_target - mass / 2)});
            max_h[i] = {-maxH, maxH};
            auto slope = [&](Int x, Int y, Int H) {
                if (x + y == 0)
                    return one;
                I X{x, x}, Y{y, y}, h{H, H};
                return divi(add(one, divi(sub(mul(X, Y), divide(sq(h), 3)), add(X, Y))),
                            mul(add(one, X), add(one, Y)));
            };
            ratios[i] = {std::max(mass / 3, slope(x.h, y.h, maxH).l), std::min(mass, slope(x.l, y.l, 0).h)};
        }
    }
    using Caps = std::array<I, 4>;
    bool inverse_caps(Caps &, const I (&)[4], const std::array<Int, 128> &) const;
    bool possible(Caps &caps, int sweeps = 2) const {
        for (int pass = 0; pass < sweeps; ++pass) {
            I H[4], B{};
            for (int i = 0; i < 4; ++i) {
                caps[i] = intersection(caps[i], mul(ratios[i], max_h[i]));
                if (caps[i].l > caps[i].h)
                    return false;
                H[i] = intersection(divi(caps[i], ratios[i]), max_h[i]);
                B = add(B, mul(tree[i], H[i]));
            }
            if (B.h < 0)
                return false;
            B = positive(B);
            I t{0, TMAX};
            if (den.l > mass / 1000000) {
                auto v = divi(mul(sq(B), B), den);
                t = intersection(t, v);
                if (t.l > t.h)
                    return false;
            }
            for (auto &c : caps) {
                c.h = std::min(c.h, 7 * volume_target - mass / 2 + t.h);
                if (c.l > c.h)
                    return false;
            }
            I corrections[6]{};
            for (int e = 0; e < 6; ++e) {
                auto [i, j] = edges[e];
                I d{-ceil_ratio(Wide(sq(sub(caps[i], caps[j])).h) * 1000000, 3 * 379464),
                    ceil_ratio(Wide(sq(add(caps[i], caps[j])).h) * 1000000, 3 * 345462)};
                auto &p = pair[e];
                if (p.usable) {
                    I v = sub(mul(p.share[0], H[i]), mul(p.share[1], H[j]));
                    // For |v| <= (A+B)/2, 2v^2/(3A) <= J(v)-J(0) <= v^2/A.
                    auto quad = divi(sq(v), p.major);
                    I factor{2 * mass / 3, mass};
                    if (p.minor.l > 0) {
                        factor.l =
                            std::max(factor.l, mass - ceil_ratio(Wide(absmax(v)) * mass, 3 * p.minor.l));
                    }
                    if (2 * absmax(v) + p.minor.h <= p.major.l)
                        factor.l = mass;
                    I exact = sub(add(mul(p.loss[0], sq(H[i])), mul(p.loss[1], sq(H[j]))), mul(quad, factor));
                    d = intersection(d, exact);
                }
                int m = 0;
                for (int r = 0; r < 7; ++r)
                    if (((reps[r] >> i) ^ (reps[r] >> j)) & 1)
                        m |= 1 << r;
                I d0{floor_ratio(lower[m], 6), mass - floor_ratio(lower[127 ^ m], 6)};
                I D = sub(d0, d), sums = add(caps[i], caps[j]), diff = sub(caps[i], caps[j]);
                if (D.l + absmin(diff) > 8 * volume_target ||
                    D.h + std::min(sums.h, 2 * t.h - sums.l) < mass - 6 * volume_target)
                    return false;
                I su = intersection(
                      sums, {mass - 6 * volume_target - D.h, D.h + 6 * volume_target + 2 * t.h - mass}),
                  di = intersection(diff, {D.l - 8 * volume_target, 8 * volume_target - D.l});
                if (su.l > su.h || di.l > di.h)
                    return false;
                caps[i] = intersection(caps[i], intersection(sub(su, caps[j]), add(di, caps[j])));
                caps[j] = intersection(caps[j], intersection(sub(su, caps[i]), sub(caps[i], di)));
                if (caps[i].l > caps[i].h || caps[j].l > caps[j].h)
                    return false;
                if (d.l > d.h)
                    return false;
                corrections[e] = d;
            }
            static const auto coeff = [] {
                std::array<std::array<int, 7>, 128> out{};
                for (int m = 1; m < 128; ++m) {
                    int k = std::countr_zero(unsigned(m));
                    out[m] = out[m ^ (1 << k)];
                    out[m][6] += k < 4 ? 1 : -1;
                    for (int h = 0; h < 6; ++h) {
                        auto [i, j] = edges[h];
                        if (k < 4)
                            out[m][h] -= i != k && j != k;
                        else
                            out[m][h] += ((reps[k] >> i) ^ (reps[k] >> j)) & 1;
                    }
                }
                return out;
            }();
            std::array<Int, 128> slack{}, Ylo{};
            for (int m = 1; m < 128; ++m) {
                I offset = scalar(t, coeff[m][6]);
                I d{};
                for (int h = 0; h < 6; ++h)
                    d = add(d, scalar(corrections[h], coeff[m][h]));
                offset = add(offset, divide(d, 2));
                Ylo[m] = floor_ratio(lower[m], 6) - offset.h;
                Int capacity = std::popcount(unsigned(m)) * 2 * volume_target;
                slack[m] = capacity - Ylo[m];
                if (slack[m] < 0)
                    return false;
            }
            for (int m = 1; m < 128; ++m) {
                slack[m] = std::min(slack[m], Wupper(m, t, Ylo));
            }
            if (mode == ProjectionMode::inverse) {
                for (int k = 0; k < 7; ++k)
                    for (int z = 0; z < 128; z += 2 << k)
                        for (int l = z; l < z + (1 << k); ++l)
                            slack[l] = std::min(slack[l], slack[l + (1 << k)]);
                if (pass == 0 && sweeps == 3 && lattice) {
                    for (int z = 1; z < 128; ++z)
                        for (int sub = (z - 1) & z; sub > z / 2; sub = (sub - 1) & z)
                            slack[z] = std::min(slack[z], slack[sub] + slack[z ^ sub]);
                    if (!inverse_caps(caps, H, slack))
                        return false;
                }
            }
            std::array<Linear4, 7> oddform{};
            std::array<I, 7> oddoffset{};
            if (pass == 0 && sweeps == 3 && errors) {
                QuadraticOffsets model(caps, H, ratios, tree, den, t, *odd, *errors, pair);
                if (model.excludes(caps, lower))
                    return false;
            }
            if (!odd->contract(caps, H, ratios, t, slack, oddform, oddoffset))
                return false;
            for (auto &r : supports) {
                Int room = 0;
                for (int k = 0; k < r.n; ++k)
                    room += r.multipliers[k] * slack[r.masks[k]];
                I value = scalar(t, r.weight);
                for (int i = 0; i < 4; ++i)
                    value = sub(value, scalar(caps[i], 2 * r.w[i]));
                if (value.l > room || value.h < -room)
                    return false;
                for (int i = 0; i < 4; ++i)
                    if (r.w[i]) {
                        I term = scalar(caps[i], -2 * r.w[i]);
                        I other{value.l - term.l, value.h - term.h};
                        I needed = divide(scalar({-room - other.h, room - other.l}, -r.w[i]), 2);
                        // value used the bounds at the start of this support test, so update simultaneously
                        // below.
                        pending[i] = intersection(caps[i], needed);
                    } else
                        pending[i] = caps[i];
                caps = pending;
                for (auto e : caps)
                    if (e.l > e.h)
                        return false;
            }
            if (pass + 1 == sweeps && dual && dual->exclude(caps, t, corrections, oddform, oddoffset))
                return false;
        }
        return true;
    }
    static Int Wupper(int m, I t, const std::array<Int, 128> &l) { return mass - t.l - l[127 ^ m]; }
    mutable Caps pending{};
};

using OddCaps = std::array<I, 4>;
using PureRows = std::array<Linear4, 4>;
// Integrate each partial derivative of the odd triple mass along the physical offset ray.
PureRows section_coefficients(const OddModel &odd, const I (&rho)[4], const I (&hh)[4]) {
    PureRows out{};
    OddCaps H{};
    std::copy(hh, hh + 4, H.begin());
    for (int i = 0; i < 4; ++i) {
        const auto &f = odd.triples[i];
        auto &row = out[i];
        for (int j = 0; j < 4; ++j)
            if (j != i)
                row[j] = {0, mass};
        if (!f.valid)
            continue;
        auto coeff = f.coeff;
        for (auto &z : f.cubic) {
            auto a = linear_value(z.lin, H), b = linear_value(z.argument, H);
            I ab = mul(a, b), bb = sq(b);
            for (int j = 0; j < 4; ++j)
                if (j != i) {
                    I val;
                    val = divide(add(mul(z.lin[j], bb), scalar(mul(z.argument[j], ab), 2)), 3);
                    coeff[j] = add(coeff[j], mul(z.coefficient, val));
                }
        }
        for (int j = 0; j < 4; ++j)
            if (j != i)
                row[j] = intersection(row[j], divi(divide(coeff[j], 12), rho[j]));
    }
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i)
            if (i != j) {
                I others{};
                for (int k = 0; k < 4; ++k)
                    if (k != i)
                        others = add(others, out[k][j]);
                out[i][j] = intersection(out[i][j], {mass - others.h, 3 * mass - others.l});
                if (out[i][j].l > out[i][j].h)
                    throw std::logic_error("inconsistent section coefficient");
            }
    return out;
}
struct SectionColumn {
    cube_section::Column poly;
    int j = 0;
    void init(int col, const PureRows &P, I r) {
        j = col;
        cube_section::Triple lo{}, hi{};
        int n = 0;
        for (int i = 0; i < 4; ++i)
            if (i != j) {
                lo[n] = P[i][j].l;
                hi[n++] = P[i][j].h;
            }
        poly.init(mass, lo, hi, r.l, r.h);
        if (poly.vertices.empty())
            poly.init(mass, {0, 0, 0}, {mass, mass, mass}, 0, mass);
    }
    std::array<Int, 4> lowest(const std::array<int, 4> &cost) const {
        std::array<int, 3> c;
        int n = 0;
        for (int i = 0; i < 4; ++i)
            if (i != j)
                c[n++] = cost[i];
        auto v = poly.lowest(c);
        std::array<Int, 4> z{};
        n = 0;
        for (int i = 0; i < 4; ++i)
            if (i != j)
                z[i] = v[n++];
        return z;
    }
};
struct ProbabilityColumns {
    std::array<SectionColumn, 4> cols;
    ProbabilityColumns(const MomentDomain &m, const OddCaps &cap, const PureRows &P) {
        I B{};
        I b[4];
        for (int j = 0; j < 4; ++j) {
            b[j] = divi(m.tree[j], m.ratios[j]);
            B = add(B, mul(b[j], cap[j]));
        }
        B = positive(B);
        I common{0, mass};
        if (m.den.l > mass / 32768)
            common = divi(sq(B), m.den);
        for (int j = 0; j < 4; ++j) {
            I r{0, ceil_ratio(Wide(224) * mass, 1000)};
            if (m.den.l > mass / 32768)
                r = intersection(r, mul(common, b[j]));
            cols[j].init(j, P, r);
        }
    }
    using BySign = std::array<std::array<std::array<Int, 4>, 2>, 4>;
    BySign supports(const std::array<int, 4> &objective) const {
        BySign out;
        for (int j = 0; j < 4; ++j) {
            out[j][1] = cols[j].lowest(objective);
            auto cost = objective;
            for (auto &v : cost)
                v = -v;
            out[j][0] = cols[j].lowest(cost);
        }
        return out;
    }
};

#include "affine_odd_marginal_bounds.hpp"
struct JointCenter {
    using Vec = std::array<Int, 7>;
    using Weights = std::array<int, 14>;
    using Pay = std::array<Int, 14>;
    Vec base{};
    std::array<std::array<Vec, 16>, 4> hull{};
    std::array<std::array<Vec, 2>, 6> overlap{};
    static Int checked(Wide x) {
        if (x > infinity || x < -infinity)
            throw std::overflow_error("oracle input");
        return (Int)x;
    }
    JointCenter(const AngularAtoms &d) {
        for (int r = 0; r < 7; ++r)
            base[r] = d.constants[r] - 12 * volume_target;
        for (int r = 0; r < 7; ++r)
            for (int j = 0; j < 6; ++j) {
                int a = 12 + j, c = d.coefficients[r][a];
                base[r] += c * d.atoms[a].center;
                overlap[j][0][r] = c * d.atoms[a].error.lo;
                overlap[j][1][r] = c * d.atoms[a].error.hi;
            }
        for (int i = 0; i < 4; ++i) {
            auto &row = d.rows[i];
            Int slope[7][2]{};
            for (int r = 0; r < 7; ++r)
                for (int j = 0; j < 6; ++j)
                    for (int u = 0; u < 2; ++u)
                        slope[r][u] += d.coefficients[r][12 + j] * d.atoms[12 + j].slope[2 * i + u];
            for (int x = 0; x <= 3; ++x)
                for (int y = 0; y <= 3; ++y) {
                    Wide numerator[7]{};
                    Int den = 0;
                    for (int a = 0; a < 2; ++a)
                        for (int b = 0; b < 2; ++b) {
                            int w = (a ? x : 3 - x) * (b ? y : 3 - y);
                            if (!w)
                                continue;
                            int k = 3 * (x - a) + y - b;
                            Int D = row.denominator[k];
                            den += w * D;
                            for (int r = 0; r < 7; ++r) {
                                Int n = 0;
                                for (int u = 0; u < 3; ++u)
                                    n += d.coefficients[r][3 * i + u] * row.numerator[u][k];
                                numerator[r] += w * (Wide(n) * mass + Wide(D) * ((2 * a - 1) * slope[r][0] +
                                                                                 (2 * b - 1) * slope[r][1]));
                            }
                        }
                    for (int r = 0; r < 7; ++r)
                        hull[i][4 * x + y][r] = den ? floor_ratio(numerator[r], den) : 0;
                    if (!den)
                        hull[i][4 * x + y][0] = infinity;
                }
        }
    }
    static Wide dot(const Vec &v, const Weights &w) {
        Wide z = 0;
        for (int r = 0; r < 7; ++r)
            z += Wide(v[r]) * (w[2 * r] + w[2 * r + 1]);
        return z;
    }
    Vec center(const Weights &w) const {
        Vec z = base;
        auto extremum = [&](auto &points) {
            int bid = -1;
            Wide best = 0;
            for (int j = 0; j < (int)points.size(); ++j) {
                if (points[j][0] == infinity)
                    continue;
                auto a = dot(points[j], w);
                if (bid < 0 || a < best) {
                    best = a;
                    bid = j;
                }
            }
            assert(bid >= 0);
            for (int r = 0; r < 7; ++r)
                z[r] += points[bid][r];
        };
        for (auto &r : hull)
            extremum(r);
        for (auto &e : overlap)
            extremum(e);
        return z;
    }
};
// Proposal only. The integer oracle separately rechecks the candidate separating vector.
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
#include "affine_deficit_projection.hpp"
#include "affine_offset_simplex.hpp"
class Checker {
  public:
    explicit Checker(Pivots pivots) {
        active = pivots;
        int coordinate = 0;
        for (int i = 0; i < 4; ++i) {
            columns_[i][0] = pivots[i];
            int pos = 1;
            for (int j = 0; j < 4; ++j)
                if (j != i && j != pivots[i])
                    columns_[i][pos++] = j;
            for (int k = 0; k < 3; ++k) {
                coordinate_[i][columns_[i][k]] = k ? coordinate++ : -1;
            }
        }
        for (int i = 0; i < 7; ++i)
            cells_[i] = combination(1U << i);
        combinations_.push_back(combination(112)); // All three pair volumes.
        for (int count = 1; count <= 6; ++count) {
            for (unsigned mask = 1; mask < 127; ++mask) {
                if (mask != 112 && std::popcount(mask) == count)
                    combinations_.push_back(combination(mask));
            }
        }
    }

    DirectionPairs angular_pairs;
    Pivots active{};
    std::array<Int, 128> physical{};
    std::optional<DirectionPairs::Interpolation> lattice;
    int exclude(const Box &box) {
        Matrix c{};
        std::array<const RowBounds *, 4> rows{};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j)
                if (i != j) {
                    const int axis = coordinate_[i][j];
                    c[i][j] = axis < 0 ? Coordinate{grid, grid} : box.variables[axis];
                }
            rows[i] = &row_bound(box.variables[2 * i], box.variables[2 * i + 1]);
        }
        std::array<Bound, 6> overlaps{};
        for (unsigned e = 0; e < edges.size(); ++e) {
            const auto [i, j] = edges[e];
            overlaps[e] = edge_overlap(c, i, j);
        }
        I normal[4][3], norm[4]{};
        for (int i = 0; i < 4; ++i)
            for (int k = 0; k < 3; ++k) {
                I t{};
                for (int j = 0; j < 4; ++j)
                    if (j != i && j != (i ^ (k + 1)))
                        t = add(t, {c[i][j].lo * (mass / grid), c[i][j].hi * (mass / grid)});
                t = divide(t, 2);
                normal[i][k] = scalar(t, i == 0 || i == k + 1 ? 1 : -1);
                norm[i] = add(norm[i], sq(t));
            }
        for (auto [i, j] : edges) {
            I d{};
            for (int k = 0; k < 3; ++k)
                d = add(d, mul(normal[i][k], normal[j][k]));
            auto d2 = sq(d), den = mul(norm[i], norm[j]);
            if (d.l > 0 && Wide(d2.l) * 1000000000000 > Wide(den.h) * 143814 * 143814)
                return 255;
            if (d.h < 0 && Wide(d2.l) * 1000000000000 > Wide(den.h) * 573370 * 573370)
                return 255;
        }
        lattice.reset();
        if (box.depth >= 6) {
            auto state = angular_pairs.prepare(c, active);
            physical = DirectionPairs::bounds(state);
            for (int m = 1; m < 127; ++m)
                if (physical[m] > thresholds[m])
                    return m;
            lattice.emplace(state);
            if (mode == ProjectionMode::inverse) {
                std::array<Int, 128> minima;
                minima.fill(infinity);
                minima[0] = 0;
                for (auto &row : lattice->corners) {
                    Int z[128]{};
                    for (int m = 1; m < 127; ++m) {
                        z[m] = z[m & (m - 1)] + row[std::countr_zero(unsigned(m))];
                        minima[m] = std::min(minima[m], z[m]);
                    }
                }
                static const auto cs = [] {
                    std::array<std::array<int, 6>, 128> a{};
                    for (int m = 1; m < 128; ++m) {
                        a[m] = a[m & (m - 1)];
                        int r = std::countr_zero(unsigned(m));
                        for (int h = 0; h < 6; ++h)
                            a[m][h] += 3 * DirectionPairs::coeff(r, h);
                    }
                    return a;
                }();
                I E[6]{};
                for (int h = 0; h < 6; ++h)
                    for (auto e : lattice->err[h])
                        E[h] = add(E[h], e);
                for (int m = 1; m < 127; ++m) {
                    Int v = minima[m];
                    for (int h = 0; h < 6; ++h) {
                        int c = cs[m][h];
                        v += c * (c > 0 ? E[h].l : E[h].h);
                    }
                    physical[m] = std::max(physical[m], v);
                    if (physical[m] > thresholds[m])
                        return m;
                }
            } else {
                constexpr int early[]{112, 1, 2, 4, 8, 16, 32, 64, 48, 80, 96, 7, 11, 13, 14};
                for (int m : early) {
                    std::array<int, 14> w{};
                    for (int r = 0; r < 7; ++r)
                        w[2 * r] = (m >> r) & 1;
                    auto z = lattice->worst(w);
                    Int v = 0;
                    for (int r = 0; r < 7; ++r)
                        if (w[2 * r])
                            v += z[r];
                    physical[m] = std::max(physical[m], v);
                    if (physical[m] > thresholds[m])
                        return m;
                }
            }
        } else
            physical.fill(-mass);
        std::optional<AngularErrors> rates;
        for (const auto &choice : combinations_) {
            Int diff = margin(choice, rows, overlaps);
            if (diff > 0)
                return choice.mask;
            Int lower = std::max(physical[choice.mask], diff + thresholds[choice.mask]);
            if (box.depth >= 6 && Wide(lower) * micro > Wide(convex_support[choice.mask].base) * mass) {
                if (!rates)
                    rates.emplace(c, volume_target);
                if (excludes(*rates, choice.mask, lower))
                    return choice.mask;
            }
        }
        return 0;
    }

    MomentDomain model(const Box &box) {
        Matrix c{};
        std::array<const RowBounds *, 4> rows{};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j)
                if (i != j) {
                    int a = coordinate_[i][j];
                    c[i][j] = a < 0 ? Coordinate{grid, grid} : box.variables[a];
                }
            rows[i] = &row_bound(box.variables[2 * i], box.variables[2 * i + 1]);
        }
        MomentDomain result(c);
        if (lattice)
            result.lattice = &*lattice;
        std::array<Bound, 6> overlaps;
        for (int h = 0; h < 6; ++h) {
            auto [i, j] = edges[h];
            overlaps[h] = edge_overlap(c, i, j);
        }
        for (auto &choice : combinations_)
            result.lower[choice.mask] =
                std::max(physical[choice.mask], margin(choice, rows, overlaps) + thresholds[choice.mask]);
        if (box.depth >= 12) {
            auto &d = result.dual.emplace();
            result.errors.emplace(c, volume_target);
            for (int i = 0; i < 4; ++i) {
                auto x = box.variables[2 * i], y = box.variables[2 * i + 1];
                d.rows[i] = row_controls(x, y);
                const auto [it, inserted] = affine_cache_.try_emplace(row_key(x, y));
                if (inserted)
                    it->second = affine_row(x, y);
                for (int k = 0; k < 3; ++k) {
                    auto t = it->second[k];
                    auto &a = d.atoms[3 * i + k];
                    a.center = t.center;
                    a.error = t.error;
                    a.slope.fill(0);
                    a.slope[2 * i] = t.slope[0];
                    a.slope[2 * i + 1] = t.slope[1];
                }
            }
            for (int h = 0; h < 6; ++h) {
                auto q = overlaps[h];
                auto &a = d.atoms[12 + h];
                a.center = (q.lo + q.hi) / 2;
                a.error = {q.lo - a.center, q.hi - a.center};
                a.slope.fill(0);
            }
            if (box.depth >= 18 && std::all_of(box.variables.begin(), box.variables.end(),
                                               [](Coordinate t) { return t.lo >= grid / 8; })) {
                auto refined = affine_overlaps(box);
                for (int h = 0; h < 6; ++h) {
                    auto &a = refined[h];
                    Int radius = 0;
                    for (auto t : a.slope)
                        radius += std::abs(t);
                    if (a.error.hi - a.error.lo + 2 * radius < 3 * (overlaps[h].hi - overlaps[h].lo))
                        d.atoms[12 + h] = a;
                }
            }
            for (int r = 0; r < 7; ++r) {
                auto &choice = cells_[r];
                d.constants[r] = choice.constant;
                for (int i = 0; i < 4; ++i) {
                    int digits = choice.row_patterns[i];
                    for (int j = 2; j >= 0; --j) {
                        d.coefficients[r][3 * i + j] = digits % 3 - 1;
                        digits /= 3;
                    }
                }
                for (int h = 0; h < 6; ++h)
                    d.coefficients[r][12 + h] = choice.overlap_coefficients[h];
            }
        }
        result.lower[127] = 6 * mass;
        return result;
    }
    // Trim all eight parameter intervals at once using a lower affine function
    // of any one parameter from the two corresponding faces of the corner grid.
    bool shave(Box &box, const std::array<Int, 128> &L) {
        if (!lattice || box.depth < 5)
            return false;
        std::array<std::pair<Int, int>, 126> priority;
        for (int m = 1; m < 127; ++m)
            priority[m - 1] = {(L[m] - thresholds[m]) / std::popcount(unsigned(m)), m};
        constexpr int K = 14;
        std::partial_sort(priority.begin(), priority.begin() + K, priority.end(), std::greater<>());
        auto base = box.variables;
        std::array<Int, 256> part{};
        bool shrunk = false;
        for (int a = 0; a < K; ++a) {
            int id = priority[a].second;
            Int err = 0;
            for (int h = 0; h < 6; ++h) {
                int k = 0;
                for (int r = 0; r < 7; ++r)
                    if (id >> r & 1)
                        k += 3 * DirectionPairs::coeff(r, h);
                if (k)
                    for (int v = 0; v < 4; ++v)
                        err += k * (k > 0 ? lattice->err[h][v].l : lattice->err[h][v].h);
            }
            for (int mask = 0; mask < 256; ++mask) {
                Int z = err;
                for (int r = 0; r < 7; ++r)
                    if (id >> r & 1)
                        z += lattice->corners[mask][r];
                part[mask] = z;
            }
            Int cap = thresholds[id];
            for (int j = 0; j < 8; ++j) {
                Int val[2]{infinity, infinity};
                int count = 128 >> j;
                for (int w = 0; w < count; ++w) {
                    val[0] = std::min(val[0], part[2 * w]);
                    val[1] = std::min(val[1], part[2 * w + 1]);
                    part[w] = std::min(part[2 * w], part[2 * w + 1]);
                }
                int width = base[j].hi - base[j].lo;
                if (width < 4)
                    continue;
                if (std::min(val[0], val[1]) > cap) {
                    box.variables[j] = {1, 0};
                    return true;
                }
                int step = std::max(1, width / 8);
                if (val[0] > cap) {
                    Int cut = floor_ratio(Wide(val[0] - cap) * width, val[0] - val[1]);
                    int trim = (int(cut) / step) * step;
                    if (trim > 0) {
                        box.variables[j].lo = std::max(box.variables[j].lo, base[j].lo + trim);
                        shrunk = true;
                    }
                }
                if (val[1] > cap) {
                    Int cut = floor_ratio(Wide(val[1] - cap) * width, val[1] - val[0]);
                    int trim = (int(cut) / step) * step;
                    if (trim > 0) {
                        box.variables[j].hi = std::min(box.variables[j].hi, base[j].hi - trim);
                        shrunk = true;
                    }
                }
            }
        }
        return shrunk;
    }

    int lattice_axis = -1;
    int probe_axis(const Box &box, const std::array<Int, 128> &L) const {
        std::array<std::pair<Int, int>, 126> rank;
        for (int m = 1; m < 127; ++m)
            rank[m - 1] = {(L[m] - thresholds[m]) / std::popcount(unsigned(m)), m};
        constexpr int K = 16;
        std::partial_sort(rank.begin(), rank.begin() + K, rank.end(), std::greater<>());
        Int widest = 0;
        for (auto t : box.variables)
            widest = std::max(widest, Int(t.hi - t.lo));
        Int value[8][2];
        for (auto &r : value)
            r[0] = r[1] = -8 * mass;
        for (int q = 0; q < K; ++q) {
            int mask = rank[q].second;
            Int errs[8]{}, cost = -thresholds[mask];
            for (int h = 0; h < 6; ++h) {
                int c = 0;
                for (int r = 0; r < 7; ++r)
                    if (mask >> r & 1)
                        c += 3 * DirectionPairs::coeff(r, h);
                for (int j = 0; j < 4; ++j) {
                    auto x = lattice->err[h][j];
                    int k = 2 * edges[h][j / 2] + j % 2;
                    errs[k] += c * (c > 0 ? x.l : x.h);
                }
            }
            for (Int e : errs)
                cost += e;
            Int f[256];
            for (int v = 0; v < 256; ++v) {
                Int z = cost;
                for (int r = 0; r < 7; ++r)
                    if (mask >> r & 1)
                        z += lattice->corners[v][r];
                f[v] = z;
            }
            for (int k = 0; k < 8; ++k) {
                Int l = infinity, r = infinity, m = infinity;
                for (int s = 0; s < 256; ++s)
                    if (!(s >> k & 1)) {
                        auto a = f[s], b = f[s ^ (1 << k)];
                        l = std::min(l, a);
                        r = std::min(r, b);
                        m = std::min(m, (a + b) / 2);
                    }
                int count = std::popcount(unsigned(mask));
                value[k][0] = std::max(value[k][0], (std::min(l, m) - errs[k] / 2) / count);
                value[k][1] = std::max(value[k][1], (std::min(r, m) - errs[k] / 2) / count);
            }
        }
        int best = -1;
        Int score = -1000 * mass;
        for (int k = 0; k < 8; ++k) {
            int width = box.variables[k].hi - box.variables[k].lo;
            if (width <= 1 || width * 6 < widest)
                continue;
            auto a = std::min(value[k][0], value[k][1]), b = std::max(value[k][0], value[k][1]);
            Int v = (a > 0 ? mass : 0) + (b > 0 ? mass / 8 : 0) + 3 * std::min<Int>(0, a) +
                    2 * std::min<Int>(0, b);
            if (v > score) {
                score = v;
                best = k;
            }
        }
        return best;
    }

    int pick_axis(const Box &box) {
        if (lattice_axis >= 0)
            return lattice_axis;
        Matrix c{};
        std::array<const RowBounds *, 4> row{};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j)
                if (j != i) {
                    int a = coordinate_[i][j];
                    c[i][j] = a < 0 ? Coordinate{grid, grid} : box.variables[a];
                }
            row[i] = &row_bound(box.variables[2 * i], box.variables[2 * i + 1]);
        }
        std::array<Bound, 6> ov;
        for (int h = 0; h < 6; ++h) {
            auto [i, j] = edges[h];
            ov[h] = edge_overlap(c, i, j);
        }
        struct Test {
            Int val;
            unsigned index;
        };
        std::array<Test, 126> tests;
        int num = 0;
        for (unsigned i = 0; i < combinations_.size(); ++i) {
            auto &ch = combinations_[i];
            tests[num++] = {(margin(ch, row, ov) + thresholds[ch.mask]) / ch.cells - 12 * volume_target, i};
        }
        constexpr int K = 6;
        std::partial_sort(tests.begin(), tests.begin() + K, tests.begin() + num,
                          [](auto a, auto b) { return a.val > b.val; });
        int axis = -1, maxwidth = 0;
        for (auto t : box.variables)
            maxwidth = std::max(maxwidth, t.hi - t.lo);
        Int best = std::numeric_limits<Int>::min();
        for (int v = 0; v < 8; ++v) {
            auto base = box.variables[v];
            int width = base.hi - base.lo, mid = (base.hi + base.lo) / 2;
            if (width <= 1 || width * 4 < maxwidth)
                continue;
            int r = v / 2, j = columns_[r][1 + v % 2];
            auto original = row[r];
            auto originalov = ov;
            Int chscore[2]{};
            for (int side = 0; side < 2; ++side) {
                auto t = side ? Coordinate{mid, base.hi} : Coordinate{base.lo, mid};
                auto a = box.variables[2 * r], b = box.variables[2 * r + 1];
                if (v % 2)
                    b = t;
                else
                    a = t;
                row[r] = &row_bound(a, b);
                c[r][j] = t;
                for (int e = 0; e < 6; ++e)
                    if (edges[e][0] != r && edges[e][1] != r)
                        ov[e] = edge_overlap(c, edges[e][0], edges[e][1]);
                chscore[side] = std::numeric_limits<Int>::min();
                for (int z = 0; z < K; ++z) {
                    auto &choice = combinations_[tests[z].index];
                    chscore[side] = std::max(
                        chscore[side], (margin(choice, row, ov) + thresholds[choice.mask]) / choice.cells -
                                           12 * volume_target);
                }
            }
            row[r] = original;
            ov = originalov;
            c[r][j] = base;
            auto high = std::max(chscore[0], chscore[1]), low = std::min(chscore[0], chscore[1]);
            Int score = (low > 0 ? mass : 0) + (high > 0 ? mass / 4 : 0) + 3 * std::min<Int>(0, low) +
                        2 * std::min<Int>(0, high);
            if (score > best) {
                axis = v;
                best = score;
            }
        }
        return axis;
    }

  private:
    std::array<std::array<int, 3>, 4> columns_{};
    std::array<std::array<int, 4>, 4> coordinate_{};
    std::vector<Combination> combinations_;
    std::array<Combination, 7> cells_{};
    std::unordered_map<std::uint64_t, RowBounds> cache_;
    std::unordered_map<std::uint64_t, std::array<AffineTriangle, 3>> affine_cache_;
    Combination combination(unsigned mask) const {
        Combination result{};
        result.mask = mask;
        result.cells = std::popcount(mask);
        result.constant = 2 * mass * std::popcount(mask >> 4);
        for (int i = 0; i < 4; ++i) {
            Pattern pattern{};
            for (int k = 0; k < 3; ++k) {
                const int j = columns_[i][k];
                pattern[k] = (mask >> j) & 1;
                for (unsigned s = 0; s < splits.size(); ++s) {
                    const auto [a, b, c, d] = splits[s];
                    if ((mask & (1U << (4 + s))) && ((i == a || i == b) == (j == a || j == b)))
                        --pattern[k];
                }
            }
            result.row_patterns[i] = pattern_id(pattern);
        }
        for (unsigned e = 0; e < edges.size(); ++e) {
            const auto [i, j] = edges[e];
            int value = -int((mask >> i) & 1) - int((mask >> j) & 1);
            for (unsigned s = 0; s < splits.size(); ++s) {
                const auto [a, b, c, d] = splits[s];
                if ((mask & (1U << (4 + s))) && ((i == a || i == b) != (j == a || j == b)))
                    ++value;
            }
            result.overlap_coefficients[e] = value;
        }
        return result;
    }

    Int margin(const Combination &choice, const std::array<const RowBounds *, 4> &rows,
               const std::array<Bound, 6> &overlaps) const {
        Int lower = choice.constant - thresholds[choice.mask];
        for (int i = 0; i < 4; ++i)
            lower += (*rows[i])[choice.row_patterns[i]];
        for (unsigned e = 0; e < edges.size(); ++e) {
            const int coefficient = choice.overlap_coefficients[e];
            lower += coefficient * (coefficient < 0 ? overlaps[e].hi : overlaps[e].lo);
        }
        return lower;
    }

    static std::uint64_t row_key(Coordinate x, Coordinate y) {
        return (static_cast<std::uint64_t>(x.lo) << 48) | (static_cast<std::uint64_t>(x.hi) << 32) |
               (static_cast<std::uint64_t>(y.lo) << 16) | static_cast<std::uint64_t>(y.hi);
    }

    std::array<AffineOverlap, 6> affine_overlaps(const Box &box) const {
        std::array<std::array<FirstDerivative, 4>, 4> box_jets{}, center_jets{};
        constexpr Int step = mass / grid;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (i != j) {
                    const int axis = coordinate_[i][j];
                    auto &jet = box_jets[i][j];
                    auto &center = center_jets[i][j];
                    if (axis < 0) {
                        jet.value = center.value = {mass, mass};
                    } else {
                        const auto [lo, hi] = box.variables[axis];
                        const Int mid = (lo + hi) * step / 2;
                        const Int radius = (hi - lo) * step / 2;
                        jet.value = {lo * step, hi * step};
                        center.value = {mid, mid};
                        jet.slope[axis] = center.slope[axis] = {radius, radius};
                    }
                }
        std::array<AffineOverlap, 6> result{};
        for (int e = 0; e < 6; ++e) {
            const auto [i, j] = edges[e];
            const auto full = overlap_jet(box_jets, i, j);
            const auto center = overlap_jet(center_jets, i, j);
            auto &model = result[e];
            model.center = (center.value.lo + center.value.hi) / 2;
            Int radius = 0;
            for (int j = 0; j < 8; ++j) {
                const Int slope = (center.slope[j].lo + center.slope[j].hi) / 2;
                model.slope[j] = slope;
                radius += std::max(std::abs(full.slope[j].lo - slope), std::abs(full.slope[j].hi - slope));
            }
            model.error = {center.value.lo - model.center - radius, center.value.hi - model.center + radius};
        }
        return result;
    }

    const RowBounds &row_bound(Coordinate x, Coordinate y) {
        const auto [entry, inserted] = cache_.try_emplace(row_key(x, y));
        auto &values = entry->second;
        if (!inserted)
            return values;
        values.fill(infinity);
        // Degree-(2,2) controls of x, y and their squares, with respective
        // common denominators 2*grid and grid^2.
        const auto controls = row_controls(x, y);
        for (int k = 0; k < 9; ++k) {
            const Int denominator = controls.denominator[k];
            // A zero denominator control has all three numerators zero.
            if (!denominator)
                continue;
            int index = 0;
            for (int a : {-1, 0, 1})
                for (int b : {-1, 0, 1})
                    for (int c : {-1, 0, 1}) {
                        const Int lower = floor_ratio(Wide{mass} * (a * controls.numerator[0][k] +
                                                                    b * controls.numerator[1][k] +
                                                                    c * controls.numerator[2][k]),
                                                      denominator);
                        values[index] = std::min(values[index], lower);
                        ++index;
                    }
        }
        if (values[13] != 0)
            throw std::logic_error("empty row control net");
        return values;
    }
};

long long game_nodes = 0, cap_nodes = 0, angular = 0, bm = 0, exact = 0, pending = 0, maxnodes;
int maxdepth, deepest = 0;
std::array<int, 4> piv;
bool fit(const MomentDomain &m, MomentDomain::Caps c, int depth) {
    ++cap_nodes;
    if (!m.possible(c, 2))
        return true;
    if (!depth)
        return false;
    int k = 0;
    for (int j = 1; j < 4; ++j)
        if (c[j].h - c[j].l > c[k].h - c[k].l)
            k = j;
    Int mid = (c[k].h + c[k].l) / 2;
    if (mid == c[k].l)
        return false;
    auto r = c;
    c[k].h = r[k].l = mid;
    return fit(m, c, depth - 1) && fit(m, r, depth - 1);
}

bool in_local_chart(const Box &b) {
    Coordinate c[4][4]{};
    int next = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (i != j)
                c[i][j] = j == piv[i] ? Coordinate{grid, grid} : b.variables[next++];
    for (int o = 0; o < 4; ++o) {
        bool ok = true;
        for (int i = 0; i < 4 && ok; ++i) {
            int close = 7;
            int main = (i ^ o) - 1;
            for (int mask = 0; mask < 8 && ok; ++mask) {
                Int n[3]{}, total = 0;
                int index = 0;
                for (int j = 0; j < 4; ++j)
                    if (j != i) {
                        int x = (mask >> index++ & 1) ? c[i][j].hi : c[i][j].lo;
                        for (int k = 0; k < 3; ++k)
                            if (j != (i ^ (1 + k)))
                                n[k] += x;
                    }
                for (auto x : n)
                    total += x;
                if (i == o) {
                    for (int k = 0; k < 3; ++k)
                        if (105 * n[k] < 34 * total || 105 * n[k] > 36 * total)
                            close &= ~(1 << k);
                    ok = std::popcount(unsigned(close)) >= 2;
                } else {
                    for (int k = 0; k < 3; ++k)
                        if (k != main)
                            ok &= 2450 * n[k] * n[k] >= 1156 * n[main] * n[main] &&
                                  2450 * n[k] * n[k] <= 1296 * n[main] * n[main];
                }
            }
        }
        if (ok)
            return true;
    }
    return false;
}

long long shavings = 0, dcuts = 0;
bool visit(Checker &checker, Box b, MomentDomain::Caps caps) {
    if (piv == std::array<int, 4>{1, 0, 0, 0}) {
        b.variables[0].hi = std::min(b.variables[0].hi, b.variables[1].hi);
        b.variables[1].lo = std::max(b.variables[1].lo, b.variables[0].lo);
        if (b.variables[0].lo > b.variables[0].hi || b.variables[1].lo > b.variables[1].hi) {
            return true;
        }
    }
    if (piv == std::array<int, 4>{1, 2, 3, 0}) {
        // Opposite-cycle column score is smallest in row zero.
        constexpr int ax[]{0, 3, 4, 6}; // initialized explicitly in the generator
        for (int i = 1; i < 4; ++i) {
            b.variables[ax[0]].hi = std::min(b.variables[ax[0]].hi, b.variables[ax[i]].hi);
            b.variables[ax[i]].lo = std::max(b.variables[ax[i]].lo, b.variables[ax[0]].lo);
        }
        for (auto t : b.variables)
            if (t.lo > t.hi) {
                return true;
            }
    }

    if (++angular > maxnodes) {
        for (auto v : b.variables)
            std::cout << v.lo << ":" << v.hi << " ";
        std::cout << " depth " << b.depth << std::endl;
        return false;
    }
    deepest = std::max(deepest, b.depth);
    checker.lattice_axis = -1;
    if (in_local_chart(b) || checker.exclude(b) ||
        std::all_of(b.variables.begin(), b.variables.end(),
                    [](Coordinate x) { return x.lo * 10 >= 7 * grid; })) {
        ++bm;
        return true;
    }
    if (b.depth >= 8 && (b.depth % 2 == 0 || b.depth == maxdepth)) {
        auto m = checker.model(b);
        ++cap_nodes;
        if (mode != ProjectionMode::plain && checker.shave(b, m.lower)) {
            if (std::any_of(b.variables.begin(), b.variables.end(),
                            [](auto r) { return r.lo > r.hi || r.hi == 0; })) {
                ++exact;
                return true;
            }
            ++shavings;
            return visit(checker, b, caps);
        }
        if (!m.possible(caps, 3)) {
            ++exact;
            return true;
        }
        if (m.dual && m.errors && b.depth >= 14 && OffsetGame(m, caps).all()) {
            ++game_nodes;
            ++exact;
            return true;
        }
        if (mode != ProjectionMode::plain && m.dual && m.errors && m.lattice && b.depth >= 12) {
            DeficitCut oracle(m, caps);
            if (oracle.ok()) {
                ++dcuts;
                ++exact;
                return true;
            }
            auto original = b;
            bool trim = false, empty = false;
            auto adjust = [&]() {
                for (int j = 0; j < 8; ++j) {
                    int width = original.variables[j].hi - original.variables[j].lo;
                    int step = std::max(1, width / 8);
                    auto &v = b.variables[j];
                    auto o = original.variables[j];
                    Int dl = floor_ratio(Wide(width) * oracle.projection.lo[j], mass),
                        du = floor_ratio(Wide(width) * (mass - oracle.projection.hi[j]), mass);
                    if (oracle.projection.hi[j] < oracle.projection.lo[j])
                        empty = true;
                    v.lo = std::max(v.lo, o.lo + int(dl) / step * step);
                    v.hi = std::min(v.hi, o.hi - int(du) / step * step);
                    trim |= v.lo != o.lo || v.hi != o.hi;
                    empty |= v.lo > v.hi || v.hi == 0;
                }
            };
            adjust();
            if (!trim && b.depth % 4 == 0) {
                std::array<int, 16> order;
                for (int j = 0; j < 16; ++j)
                    order[j] = j;
                std::partial_sort(order.begin(), order.begin() + 3, order.end(), [&](int a, int b) {
                    return oracle.conditional[a] > oracle.conditional[b];
                });
                for (int j = 0; j < 3 && !trim; ++j) {
                    oracle.attempt(order[j], 7);
                    adjust();
                }
            }
            if (empty) {
                ++dcuts;
                ++exact;
                return true;
            }
            if (trim) {
                ++shavings;
                return visit(checker, b, caps);
            }
        }
        int capdepth = b.depth < 16 ? 2 : b.depth < 24 ? 4 : b.depth < 32 ? 6 : 9;
        if (b.depth == maxdepth)
            capdepth = 18;
        if (fit(m, caps, capdepth)) {
            ++exact;
            return true;
        }
        if (mode != ProjectionMode::plain && checker.lattice)
            checker.lattice_axis = checker.probe_axis(b, m.lower);
    }
    if (b.depth == maxdepth) {
        ++pending;
        for (auto x : b.variables)
            std::cout << x.lo << ":" << x.hi << " ";
        std::cout << std::endl;
        return false;
    }
    int axis = checker.pick_axis(b);
    if (axis < 0)
        return false;
    auto rhs = b;
    int mid = (b.variables[axis].lo + b.variables[axis].hi) / 2;
    b.variables[axis].hi = mid;
    rhs.variables[axis].lo = mid;
    ++b.depth;
    ++rhs.depth;
    return visit(checker, b, caps) && visit(checker, rhs, caps);
}
} // namespace
int main(int argc, char **argv) {
    if (argc < 2 || argc > 7)
        return 2;
    std::string_view word = argv[1];
    if (word != "1000" && word != "1001" && word != "1002" && word != "1200")
        return 2;
    for (int i = 0; i < 4; ++i)
        piv[i] = word[i] - '0';
    maxdepth = argc > 2 ? std::atoi(argv[2]) : 96;
    maxnodes = argc > 3 ? std::atoll(argv[3]) : 200000000;
    if (maxdepth < 0 || maxdepth > 112 || maxnodes < 1)
        return 2;
    if (argc >= 5) {
        std::string_view algorithm = argv[4];
        if (algorithm == "project")
            mode = ProjectionMode::project;
        else if (algorithm == "columns")
            mode = ProjectionMode::columns;
        else if (algorithm == "flow")
            mode = ProjectionMode::flow;
        else if (algorithm == "inverse")
            mode = ProjectionMode::inverse;
        else if (algorithm != "plain")
            return 2;
    }
    for (auto &x : thresholds)
        if (!(std::cin >> x))
            return 2;
    Checker checker(piv);
    Box b{};
    b.variables.fill({0, grid});
    if (argc >= 6) {
        // Sixty-four closed sectors cover the full angular domain (binary digits correspond to coordinates).
        int sector = -1;
        std::string_view arg = argv[5];
        auto result = std::from_chars(arg.begin(), arg.end(), sector);
        if (result.ec != std::errc{} || result.ptr != arg.end() || sector < 0 || sector >= 64)
            return 2;
        constexpr int axes[6]{0, 2, 4, 6, 3, 7};
        for (int j = 0; j < 6; ++j)
            b.variables[axes[j]] = (sector >> j & 1) ? Coordinate{grid / 2, grid} : Coordinate{0, grid / 2};
        b.depth = 6;
        std::cout << "sector " << sector << " ";
        if (argc == 7) {
            int patch = -1;
            std::string_view arg = argv[6];
            auto rr = std::from_chars(arg.begin(), arg.end(), patch);
            if (rr.ec != std::errc{} || rr.ptr != arg.end() || patch < 0 || patch >= 64)
                return 2;
            constexpr int extra[6]{1, 5, 0, 2, 4, 6};
            for (int j = 0; j < 6; ++j) {
                auto &v = b.variables[extra[j]];
                int mid = (v.lo + v.hi) / 2;
                if (patch >> j & 1)
                    v.lo = mid;
                else
                    v.hi = mid;
            }
            b.depth += 6;
            std::cout << "patch " << patch << " ";
        }
    }
    MomentDomain::Caps c;
    c.fill({-(7 * volume_target - mass / 2), 7 * volume_target - mass / 2 + TMAX});
    bool ok = visit(checker, b, c);
    std::cout << word << " " << (ok ? "complete" : "OPEN") << " angular " << angular << " cap " << cap_nodes
              << " bm " << bm << " moment " << exact << " deepest " << deepest << " game " << game_nodes
              << " unresolved " << pending << "\n";
    if (mode != ProjectionMode::plain)
        std::cout << "contractions " << shavings << " deficit " << dcuts << "\n";
    return !ok;
}
