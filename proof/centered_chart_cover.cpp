// Exact direct comparison to (1+2 sqrt(2))/48 in a centered vertex chart.
// Build with c++ -O3 -std=c++20 centered_chart_cover.cpp -o centered_chart_cover.
// Usage: centered_chart_cover 1032 [depth] [node-budget] [triangle-start] [overlap-start]
// A nonzero exit code means the requested cover was not completed.

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
struct Box {
    std::array<Coordinate, 8> variables;
    int depth;
};

Int floor_ratio(Wide numerator, Int denominator) {
    if (denominator <= 0)
        throw std::logic_error("nonpositive denominator");
    return numerator < 0 ? -((-numerator + denominator - 1) / denominator) : numerator / denominator;
}

Int ceil_ratio(Wide numerator, Int denominator) { return -floor_ratio(-numerator, denominator); }

Int candidate_upper() {
    Int lower = mass / 2, upper = mass;
    while (upper - lower > 1) {
        const Int middle = (lower + upper) / 2;
        const Wide shifted = 4 * middle - mass;
        if (shifted * shifted > Wide{8} * mass * mass)
            upper = middle;
        else
            lower = middle;
    }
    return upper;
}

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
using Pivots = std::array<int, 4>;
using Pattern = std::array<int, 3>;
using RowBounds = std::array<Int, 27>;

constexpr std::array<std::array<int, 2>, 6> edges{{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
constexpr std::array<std::array<int, 4>, 3> splits{{{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}}};

#include "centered_chart_terms.hpp"

Weights propose_weights(const std::array<Int, 7> &intercepts, const std::array<LinearSlopes, 7> &slopes) {
    // Maximize w·intercepts minus the sum of |w·slope_j|. The auxiliary
    // variables enclose the absolute values. All right sides are nonnegative.
    constexpr int variables = 15, constraints = 17, columns = variables + constraints + 1;
    std::array<std::array<double, columns>, constraints + 1> tableau{};
    std::array<int, constraints> basic{};
    for (int row = 0; row < constraints; ++row) {
        tableau[row][variables + row] = 1.;
        basic[row] = variables + row;
    }
    tableau[0].back() = 1.;
    for (int i = 0; i < 7; ++i) {
        tableau[0][i] = 1.;
        tableau.back()[i] = -intercepts[i] / static_cast<double>(mass);
    }
    for (int j = 0; j < 8; ++j) {
        tableau.back()[7 + j] = 1.;
        for (int sign : {-1, 1}) {
            auto &row = tableau[1 + 2 * j + (sign + 1) / 2];
            row[7 + j] = -1.;
            for (int i = 0; i < 7; ++i)
                row[i] = sign * slopes[i][j] / static_cast<double>(mass);
        }
    }
    for (int step = 0; step < 100; ++step) {
        int enter = 0;
        while (enter < columns - 1 && tableau.back()[enter] >= -1e-11)
            ++enter;
        if (enter == columns - 1) {
            Weights result{};
            for (int row = 0; row < constraints; ++row)
                if (basic[row] < 7) {
                    const double value = tableau[row].back();
                    if (!std::isfinite(value) || value > 2.)
                        return {};
                    result[basic[row]] = std::llround(std::max(0., value) * 65536.);
                }
            return result;
        }
        int leave = -1;
        double best = std::numeric_limits<double>::infinity();
        for (int row = 0; row < constraints; ++row)
            if (tableau[row][enter] > 1e-11) {
                const double ratio = std::max(0., tableau[row].back()) / tableau[row][enter];
                if (ratio < best - 1e-11 ||
                    (leave >= 0 && std::abs(ratio - best) <= 1e-11 && basic[row] < basic[leave])) {
                    leave = row;
                    best = ratio;
                }
            }
        if (leave < 0)
            return {};
        auto &pivot = tableau[leave];
        const double reciprocal = 1. / pivot[enter];
        for (double &entry : pivot)
            entry *= reciprocal;
        basic[leave] = enter;
        for (int row = 0; row <= constraints; ++row)
            if (row != leave) {
                const double factor = tableau[row][enter];
                if (factor == 0.)
                    continue;
                for (int j = 0; j < columns; ++j)
                    tableau[row][j] -= factor * pivot[j];
                tableau[row][enter] = 0.;
            }
    }
    return {};
}

class Checker {
  public:
    explicit Checker(Pivots pivots, int affine_depth, int overlap_depth, bool moment_mode = false)
        : affine_depth_(affine_depth), overlap_depth_(overlap_depth), moment_mode_(moment_mode) {
        if (moment_mode_) {
            for (auto &v : moment_thresholds_)
                if (!(std::cin >> v))
                    throw std::invalid_argument("expected 128 certified moment thresholds");
            if (moment_thresholds_.front() != 0 || moment_thresholds_.back() != 6 * mass)
                throw std::invalid_argument("invalid moment thresholds");
        }
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
        for (const auto &choice : combinations_) {
            if (margin(choice, rows, overlaps) > 0)
                return choice.mask;
        }
        if (!moment_mode_ && box.depth >= affine_depth_ && exclude_affine(box, overlaps))
            return 128;
        return 0;
    }

  private:
    std::array<std::array<int, 3>, 4> columns_{};
    std::array<std::array<int, 4>, 4> coordinate_{};
    std::vector<Combination> combinations_;
    std::array<Combination, 7> cells_{};
    std::unordered_map<std::uint64_t, RowBounds> cache_;
    std::unordered_map<std::uint64_t, std::array<AffineTriangle, 3>> affine_cache_;
    const int affine_depth_;
    const int overlap_depth_;
    const bool moment_mode_;
    std::array<Int, 128> moment_thresholds_{};
    const Int target_ = candidate_upper();

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
        Int lower =
            choice.constant - (moment_mode_ ? moment_thresholds_[choice.mask] : target_ * choice.cells);
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

    bool exclude_affine(const Box &box, const std::array<Bound, 6> &overlaps) {
        std::array<const std::array<AffineTriangle, 3> *, 4> triangles{};
        for (int i = 0; i < 4; ++i) {
            const auto x = box.variables[2 * i], y = box.variables[2 * i + 1];
            const auto [entry, inserted] = affine_cache_.try_emplace(row_key(x, y));
            if (inserted)
                entry->second = affine_row(x, y);
            triangles[i] = &entry->second;
        }
        std::array<AffineOverlap, 6> models{};
        for (int e = 0; e < 6; ++e) {
            const auto [lo, hi] = overlaps[e];
            const Int center = (lo + hi) / 2;
            models[e] = {center, {}, {lo - center, hi - center}};
        }
        if (check_affine(triangles, models))
            return true;
        if (box.depth < overlap_depth_ || std::any_of(
                                              box.variables.begin(), box.variables.end(),
                                              [](Coordinate interval) { return interval.lo < grid / 8; }))
            return false;
        return check_affine(triangles, affine_overlaps(box));
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

    bool check_affine(const std::array<const std::array<AffineTriangle, 3> *, 4> &triangles,
                      const std::array<AffineOverlap, 6> &overlaps) const {
        std::array<Int, 7> intercepts{};
        std::array<LinearSlopes, 7> slopes{};
        for (int cell = 0; cell < 7; ++cell) {
            const auto &choice = cells_[cell];
            Int intercept = choice.constant - target_;
            for (int i = 0; i < 4; ++i) {
                int digits = choice.row_patterns[i];
                for (int j = 2; j >= 0; --j) {
                    const int coefficient = digits % 3 - 1;
                    digits /= 3;
                    const auto &term = (*triangles[i])[j];
                    intercept +=
                        coefficient * (term.center + (coefficient < 0 ? term.error.hi : term.error.lo));
                    for (int k = 0; k < 2; ++k)
                        slopes[cell][2 * i + k] += coefficient * term.slope[k];
                }
            }
            for (int e = 0; e < 6; ++e) {
                const int coefficient = choice.overlap_coefficients[e];
                const auto &atom = overlaps[e];
                intercept += coefficient * (atom.center + (coefficient < 0 ? atom.error.hi : atom.error.lo));
                for (int j = 0; j < 8; ++j)
                    slopes[cell][j] += coefficient * atom.slope[j];
            }
            intercepts[cell] = intercept;
            for (Int slope : slopes[cell])
                intercept -= std::abs(slope);
            if (intercept > 0)
                return true;
        }
        const Weights weights = propose_weights(intercepts, slopes);
        Wide lower = 0;
        std::array<Wide, 8> slopes_total{};
        for (int cell = 0; cell < 7; ++cell)
            lower += Wide{weights[cell]} * (cells_[cell].constant - target_);
        for (int i = 0; i < 4; ++i) {
            std::array<Wide, 2> slope{};
            int divisor = 9;
            for (int j = 0; j < 3; ++j) {
                Int coefficient = 0;
                for (int cell = 0; cell < 7; ++cell) {
                    coefficient += weights[cell] * (cells_[cell].row_patterns[i] / divisor % 3 - 1);
                }
                divisor /= 3;
                const auto &atom = (*triangles[i])[j];
                lower +=
                    Wide{coefficient} * (atom.center + (coefficient < 0 ? atom.error.hi : atom.error.lo));
                for (int k = 0; k < 2; ++k)
                    slope[k] += Wide{coefficient} * atom.slope[k];
            }
            for (int k = 0; k < 2; ++k)
                slopes_total[2 * i + k] += slope[k];
        }
        for (int e = 0; e < 6; ++e) {
            Int coefficient = 0;
            for (int cell = 0; cell < 7; ++cell)
                coefficient += weights[cell] * cells_[cell].overlap_coefficients[e];
            const auto &atom = overlaps[e];
            lower += Wide{coefficient} * (atom.center + (coefficient < 0 ? atom.error.hi : atom.error.lo));
            for (int j = 0; j < 8; ++j)
                slopes_total[j] += Wide{coefficient} * atom.slope[j];
        }
        for (const Wide term : slopes_total)
            lower -= term < 0 ? -term : term;
        return lower > 0;
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

template <class Number> Number parse(std::string_view word) {
    Number result{};
    const auto [end, error] = std::from_chars(word.data(), word.data() + word.size(), result);
    if (error != std::errc{} || end != word.data() + word.size())
        throw std::invalid_argument("invalid integer");
    return result;
}

int candidate_radius(const Box &box, const Pivots &pivots) {
    int root = 0, upper = grid;
    while (upper - root > 1) {
        const int mid = (root + upper) / 2;
        const Int sum = grid + mid;
        if (sum * sum < 2 * Int{grid} * grid)
            root = mid;
        else
            upper = mid;
    }
    int best = grid + 1;
    for (int small = 0; small < 4; ++small) {
        if (std::find(pivots.begin(), pivots.end(), small) != pivots.end())
            continue;
        int axis = 0, radius = 0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                if (i == j || j == pivots[i])
                    continue;
                const auto [lo, hi] = box.variables[axis++];
                const int a = j == small ? root : grid;
                const int b = j == small ? upper : grid;
                radius = std::max({radius, b - lo, hi - a});
            }
        best = std::min(best, radius);
    }
    return best;
}

} // namespace

int main(int argc, char **argv) {
    const bool moment_mode = argc > 1 && std::string_view(argv[1]) == "--affine-moments";
    if (moment_mode) {
        --argc;
        ++argv;
    }
    const std::string_view chart = argc > 1 ? argv[1] : "1032";
    const int max_depth = argc > 2 ? parse<int>(argv[2]) : 32;
    const auto limit = argc > 3 ? parse<std::uint64_t>(argv[3]) : 10000000;
    const int affine_depth = argc > 4 ? parse<int>(argv[4]) : max_depth + 1;
    const int overlap_depth = argc > 5 ? parse<int>(argv[5]) : max_depth + 1;
    if (chart.size() != 4 || max_depth < 0 || max_depth > 8 * 15)
        return 2;
    Pivots pivots{};
    for (int i = 0; i < 4; ++i) {
        pivots[i] = chart[i] - '0';
        if (pivots[i] < 0 || pivots[i] > 3 || pivots[i] == i)
            return 2;
    }
    Checker checker(pivots, affine_depth, overlap_depth, moment_mode);
    Box initial{};
    initial.variables.fill({0, grid});
    std::vector<Box> stack{initial};
    std::uint64_t visited = 0, individual = 0, combined = 0, affine = 0, unresolved = 0;
    int deepest = 0, root_mask = 0, open_radius = 0;
    Wide covered = 0;
    while (!stack.empty() && visited < limit) {
        Box box = stack.back();
        stack.pop_back();
        ++visited;
        deepest = std::max(deepest, box.depth);
        const bool regular = moment_mode && std::all_of(box.variables.begin(), box.variables.end(),
                                                        [](Coordinate v) { return 10 * v.lo >= 7 * grid; });
        const int verdict = regular ? 128 : checker.exclude(box);
        if (verdict) {
            if (verdict == 128)
                ++affine;
            else
                ++(std::popcount(static_cast<unsigned>(verdict)) == 1 ? individual : combined);
            if (box.depth == 0)
                root_mask = verdict;
            covered += Wide{1} << (max_depth - box.depth);
        } else if (box.depth == max_depth) {
            ++unresolved;
            open_radius = std::max(open_radius, candidate_radius(box, pivots));
        } else {
            const int axis = box.depth % 8;
            const auto [lo, hi] = box.variables[axis];
            const int midpoint = (lo + hi) / 2;
            Box right = box;
            ++box.depth;
            ++right.depth;
            box.variables[axis].hi = midpoint;
            right.variables[axis].lo = midpoint;
            stack.push_back(right);
            stack.push_back(box);
        }
    }
    const bool complete = unresolved == 0 && stack.empty() && covered == (Wide{1} << max_depth);
    std::cout << "chart=" << chart << " visited=" << visited << " individual=" << individual
              << " combined=" << combined << (moment_mode ? " regular=" : " affine=") << affine
              << " unresolved=" << unresolved << " pending=" << stack.size() << " depth=" << deepest
              << " root-mask=" << root_mask << " open-radius=" << open_radius << '/' << grid
              << (complete ? " complete\n" : " incomplete\n");
    return complete ? 0 : 1;
}
