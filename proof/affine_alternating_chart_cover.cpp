// Complete affine offset cover for the four-cycle of maximal columns.
#include "affine_alternating_moments.hpp"
#include "affine_box_linear.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
using namespace cube_moments;
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
constexpr Int volume_target = Q;
#include "affine_alternating_supports.hpp"
struct MomentDomain : Geometry {
    std::optional<AngularAtoms> dual;
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
            Int maxH = std::min({x.h, y.h, 3 * (TMAX + 7 * Q - mass / 2)});
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
                Int capacity = std::popcount(unsigned(m)) * 2 * Q;
                slack[m] = capacity - Ylo[m];
                if (slack[m] < 0)
                    return false;
            }
            for (int m = 1; m < 128; ++m) {
                slack[m] = std::min(slack[m], Wupper(m, t, Ylo));
            }
            std::array<Linear4, 7> oddform{};
            std::array<I, 7> oddoffset{};
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

class Checker {
  public:
    explicit Checker(Pivots pivots) {
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
        std::optional<AngularErrors> rates;
        for (const auto &choice : combinations_) {
            Int diff = margin(choice, rows, overlaps);
            if (diff > 0)
                return choice.mask;
            Int lower = diff + thresholds[choice.mask];
            if (box.depth >= 6 && Wide(lower) * micro > Wide(convex_support[choice.mask].base) * mass) {
                if (!rates)
                    rates.emplace(c);
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
        std::array<Bound, 6> overlaps;
        for (int h = 0; h < 6; ++h) {
            auto [i, j] = edges[h];
            overlaps[h] = edge_overlap(c, i, j);
        }
        for (auto &choice : combinations_)
            result.lower[choice.mask] = margin(choice, rows, overlaps) + thresholds[choice.mask];
        if (box.depth >= 12) {
            auto &d = result.dual.emplace();
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
            if (box.depth >= 28 && std::all_of(box.variables.begin(), box.variables.end(),
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
    int pick_axis(const Box &box) {
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
            tests[num++] = {margin(ch, row, ov) / ch.cells, i};
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
                    if (edges[e][0] == r || edges[e][1] == r)
                        ov[e] = edge_overlap(c, edges[e][0], edges[e][1]);
                chscore[side] = std::numeric_limits<Int>::min();
                for (int z = 0; z < K; ++z) {
                    auto &choice = combinations_[tests[z].index];
                    chscore[side] = std::max(chscore[side], margin(choice, row, ov) / choice.cells);
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

} // namespace
long long cap_nodes = 0, angular = 0, bm = 0, exact = 0, pending = 0, maxnodes;
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
bool visit(Checker &checker, Box b, MomentDomain::Caps caps) {
    if (piv == std::array<int, 4>{1, 2, 3, 0}) {
        // Opposite-cycle column score is smallest in row zero.
        constexpr int ax[]{0, 3, 4, 6};
        for (int i = 1; i < 4; ++i) {
            b.variables[ax[0]].hi = std::min(b.variables[ax[0]].hi, b.variables[ax[i]].hi);
            b.variables[ax[i]].lo = std::max(b.variables[ax[i]].lo, b.variables[ax[0]].lo);
        }
        for (auto t : b.variables)
            if (t.lo > t.hi) {
                return true;
            }
    }

    if (++angular > maxnodes)
        return false;
    deepest = std::max(deepest, b.depth);
    if (checker.exclude(b) || std::all_of(b.variables.begin(), b.variables.end(),
                                          [](Coordinate x) { return x.lo * 10 >= 7 * grid; })) {
        ++bm;
        return true;
    }
    if (b.depth >= 8 && (b.depth % 2 == 0 || b.depth == maxdepth)) {
        auto m = checker.model(b);
        ++cap_nodes;
        if (!m.possible(caps, 2)) {
            ++exact;
            return true;
        }
        int capdepth = b.depth < 16 ? 1 : b.depth < 24 ? 3 : b.depth < 32 ? 5 : 7;
        if (b.depth == maxdepth)
            capdepth = 18;
        if (fit(m, caps, capdepth)) {
            ++exact;
            return true;
        }
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

int main(int argc, char **argv) {
    if (argc > 3)
        return 2;
    maxdepth = argc > 1 ? std::atoi(argv[1]) : 80;
    maxnodes = argc > 2 ? std::atoll(argv[2]) : 200000000;
    if (maxdepth < 0 || maxdepth > 112 || maxnodes < 1)
        return 2;
    for (auto &x : thresholds)
        if (!(std::cin >> x))
            return 2;
    piv = {1, 2, 3, 0};
    Checker checker(piv);
    Box b{};
    b.variables.fill({0, grid});
    MomentDomain::Caps c;
    c.fill({-(7 * Q - mass / 2), 7 * Q - mass / 2 + TMAX});
    bool ok = visit(checker, b, c);
    std::cout << "1230 " << (ok ? "complete" : "OPEN") << " angular " << angular << " cap " << cap_nodes
              << " bm " << bm << " moment " << exact << " deepest " << deepest << " unresolved " << pending
              << "\n";
    return !ok;
}
