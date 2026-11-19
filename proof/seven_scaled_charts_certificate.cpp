// Exact continuous-box verifier for the three-plane, seven-cell case.
// Build: c++ -std=c++20 -O3 seven_scaled_charts_certificate.cpp -o seven_scaled_charts_certificate
// Complete cover:
//   ./seven_scaled_charts_certificate --horizontal 30
//   ./seven_scaled_charts_certificate --vertical 36
// The two runs jointly prove M > 4/25 for every three-plane arrangement with
// seven positive-volume cells. Exit status 1 means the printed cover is incomplete.
// The default command checks a smaller horizontal subdomain quickly.
//
// Orient the three planes so the missing sign is +++. A minimizer of the Farkas
// witness function on the weight simplex can be chosen at two coordinate-zero
// break lines. Pairwise overlap excludes the boundary of that simplex. After
// absorbing the weights and applying cube symmetries, the rows sum to r*e_z and
// their largest absolute coefficient is 1. The two charts are the cases where
// this coefficient is horizontal and vertical, respectively.
//
// With seven cells and maximum <= 4/25, every cell has mass >= 1/25. If delta_i
// denotes 1/2 minus the positive cap mass, then .02 <= delta_i <= .14 and
// Delta=sum delta_i <= .3. The all-negative mass is at least .01+Delta/2.
// Barthe--Koldobsky's slab theorem implies h_i <= delta_i*||n_i||_2: first use
// monotonicity to put 2*h_i/||n_i||_2 below 3/4, then apply the slab bound.
// https://doi.org/10.1016/S0001-8708(02)00055-5
// Brunn--Minkowski applied along common pair-offset scaling implies the central
// +++ mass w >= .005+.099*Delta, hence w >= .01094. Each cell volume to the
// power 1/3 is also concave on the feasible offset polytope when the vertex
// cells have positive volume. Thus a common large cell at all vertices suffices.

#include "seven_integer_volume.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

using seven_integer::Bound;
using seven_integer::Bounds;
using seven_integer::Int;
using seven_integer::Vector;
using seven_integer::Wide;

using Interval = std::array<Int, 2>;
using Parameters = std::array<Interval, 6>;
using Normals = std::array<std::array<Interval, 3>, 3>;

constexpr Int S = seven_integer::parameter_grid;
constexpr Int G = seven_integer::space_grid;
constexpr Int V = seven_integer::volume_grid;
constexpr Int WG = Int{1} << 18;
constexpr Int R = Int{1} << 16;

constexpr Int w_min_numerator = 547;
constexpr Int w_min_denominator = 50000;

constexpr Bounds cube = {{{-G / 2, G / 2}, {-G / 2, G / 2}, {-G / 2, G / 2}}};

Int floor_div(Wide a, Wide b) { return Int(a / b - (a % b < 0)); }
Int ceil_div(Wide a, Wide b) { return -floor_div(-a, b); }

Int integer_sqrt(Int value) {
    if (value <= 1)
        return value;
    Int estimate = Int{1} << ((std::bit_width(static_cast<std::uint64_t>(value)) + 1) / 2);
    while (true) {
        const Int next = (estimate + value / estimate) / 2;
        if (next >= estimate)
            return estimate;
        estimate = next;
    }
}

Int sqrt_upper(Int value) {
    const Int root = integer_sqrt(value);
    return root + (root * root < value);
}

Interval add(Interval a, Interval b) { return {a[0] + b[0], a[1] + b[1]}; }
Interval negate(Interval a) { return {-a[1], -a[0]}; }
Interval subtract(Interval a, Interval b) { return add(a, negate(b)); }

Interval magnitude(Interval a) {
    return {std::max({Int{0}, a[0], -a[1]}), std::max(std::abs(a[0]), std::abs(a[1]))};
}

Interval multiply(Interval a, Interval b) {
    const std::array<Int, 4> terms = {a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]};
    return {*std::min_element(terms.begin(), terms.end()), *std::max_element(terms.begin(), terms.end())};
}

bool restrict_to(Interval &a, Interval b) {
    a = {std::max(a[0], b[0]), std::min(a[1], b[1])};
    return a[0] <= a[1];
}

bool sum_in(Interval &a, Interval &b, Interval target) {
    return restrict_to(a, subtract(target, b)) && restrict_to(b, subtract(target, a));
}

enum class Chart { Horizontal, Vertical };

bool unpack(Parameters &p, Chart chart, Normals &n, Interval &vertical) {
    // Σt_i <= (3/10)*sqrt(3) for any potential counterexample.
    const Int max_sum = ceil_div(3 * sqrt_upper(3 * S * S), 5);
    if (chart == Chart::Horizontal) {
        // p1=(1,a), p2=(b,c), p3=(-1-b,-a-c). Reflect z to make its sum nonnegative.
        if (!sum_in(p[0], p[2], {-S, S}))
            return false;
        const Interval target = {0, max_sum};
        if (!restrict_to(p[3], subtract(target, add(p[4], p[5]))) ||
            !restrict_to(p[4], subtract(target, add(p[3], p[5]))) ||
            !restrict_to(p[5], subtract(target, add(p[3], p[4]))))
            return false;
        n[0] = {Interval{S, S}, p[0], p[3]};
        n[1] = {p[1], p[2], p[4]};
        n[2] = {subtract({-S, -S}, p[1]), negate(add(p[0], p[2])), p[5]};
        if (!restrict_to(n[2][1], {-S, S}))
            return false;
        vertical = add(p[3], add(p[4], p[5]));
        return restrict_to(vertical, target);
    }
    // n1=(a,b,1), 0<=b<=a<=1; order z2>=z3 by swapping the remaining cuts.
    if (!restrict_to(p[1], {0, p[0][1]}) || !restrict_to(p[0], {p[1][0], S}) ||
        !sum_in(p[0], p[2], {-S, S}) || !sum_in(p[1], p[3], {-S, S}) || !restrict_to(p[4], {p[5][0], S}) ||
        !restrict_to(p[5], {-S, p[4][1]}))
        return false;
    const Interval target = {-max_sum, max_sum};
    if (!sum_in(p[4], p[5], subtract(target, {S, S})))
        return false;
    n[0] = {p[0], p[1], Interval{S, S}};
    n[1] = {p[2], p[3], p[4]};
    n[2] = {negate(add(p[0], p[2])), negate(add(p[1], p[3])), p[5]};
    if (!restrict_to(n[2][0], {-S, S}) || !restrict_to(n[2][1], {-S, S}))
        return false;
    vertical = add({S, S}, add(p[4], p[5]));
    return restrict_to(vertical, target);
}

struct Rational {
    Wide num;
    Wide den;
};

Rational cap(Vector coefficients, Int h) {
    Int total = 0;
    Wide denominator = 1;
    int dimension = 0;
    for (Int w : coefficients)
        if (w) {
            coefficients[dimension++] = w;
            total += w;
            denominator *= 2 * w * dimension;
        }
    if (!dimension)
        return {h < 0 ? 1 : 0, 1};
    const Int limit = total - 2 * h;
    if (limit <= 0)
        return {0, 1};
    if (limit >= 2 * total)
        return {1, 1};
    Wide numerator = 0;
    for (int mask = 0; mask < (1 << dimension); ++mask) {
        Int x = limit, sign = 1;
        for (int j = 0; j < dimension; ++j)
            if (mask & (1 << j)) {
                x -= 2 * coefficients[j];
                sign = -sign;
            }
        if (x > 0) {
            Wide contribution = sign;
            for (int j = 0; j < dimension; ++j)
                contribution *= x;
            numerator += contribution;
        }
    }
    return {numerator, denominator};
}

Interval quantile(Vector coefficients, Int q, Int denominator) {
    std::sort(coefficients.begin(), coefficients.end());
    std::uint64_t key = static_cast<std::uint64_t>(q * 32 + denominator);
    for (Int w : coefficients)
        key = (key << 17) | static_cast<std::uint64_t>(w);
    struct CacheEntry {
        std::uint64_t key;
        Interval value;
    };
    static std::array<CacheEntry, 1 << 17> cache{};
    const auto slot = static_cast<std::size_t>((key * 11400714819323198485ull) >> (64 - 17));
    if (cache[slot].key == key)
        return cache[slot].value;
    Int lo = 0, hi = (coefficients[0] + coefficients[1] + coefficients[2] + 1) / 2;
    while (hi - lo > 1) {
        const Int midpoint = (lo + hi) / 2;
        const Rational p = cap(coefficients, midpoint);
        if (p.num * denominator > p.den * q)
            lo = midpoint;
        else
            hi = midpoint;
    }
    cache[slot] = {key, {lo, hi}};
    return {lo, hi};
}

Int max_component(const Vector &a) { return *std::max_element(a.begin(), a.end()); }

Interval minor_at(const Normals &a, int i, int j, int axis) {
    const int x = (axis + 1) % 3, y = (axis + 2) % 3;
    return subtract(multiply(a[i][x], a[j][y]), multiply(a[i][y], a[j][x]));
}

Wide positive_cube(Int x) { return x <= 0 ? 0 : Wide(x) * x * x; }
Wide positive_square(Int x) { return x <= 0 ? 0 : Wide(x) * x; }

Rational positive_mean(Int a, Int b, Int threshold) {
    if (a < b)
        std::swap(a, b);
    const Int r = a + b - 2 * threshold;
    if (r <= 0)
        return {0, 1};
    if (b == 0)
        return {positive_square(r) - positive_square(r - 2 * a), Wide(8) * a};
    return {positive_cube(r) - positive_cube(r - 2 * a) - positive_cube(r - 2 * b) +
                positive_cube(r - 2 * a - 2 * b),
            Wide(48) * a * b};
}

Vector pair_strip_upper(const Normals &normals, const Vector &h_lo, Int planar_det_hi) {
    Vector bounds = {ceil_div(4 * V, 25), ceil_div(4 * V, 25), ceil_div(4 * V, 25)};
    for (int first = 0; first < 3; ++first)
        for (int second = first + 1; second < 3; ++second)
            for (int pivot = 0; pivot < 3; ++pivot) {
                int i = first, j = second;
                if (normals[i][pivot][1] < 0)
                    std::swap(i, j);
                if (normals[i][pivot][0] <= 0 || normals[j][pivot][1] >= 0)
                    continue;
                const Interval a = normals[i][pivot], b = negate(normals[j][pivot]);
                std::array<Int, 2> upper{};
                int count = 0;
                for (int k = 0; k < 3; ++k)
                    if (k != pivot) {
                        // For two horizontal coefficients the width uses the common
                        // projected determinant exactly, even near opposite normals.
                        if (pivot < 2 && k < 2)
                            upper[count++] = planar_det_hi;
                        else
                            upper[count++] =
                                magnitude(add(multiply(b, normals[i][k]), multiply(a, normals[j][k])))[1];
                    }
                const Rational r = positive_mean(upper[0], upper[1], b[0] * h_lo[i] + a[0] * h_lo[j]);
                const int missing = 3 - first - second;
                bounds[missing] = std::min(bounds[missing], ceil_div(r.num * V, r.den * a[0] * b[0]));
            }
    return bounds;
}

struct Inequality {
    Vector a;
    Int b;
};

struct Offset {
    Vector numerator;
    Int denominator;
    auto operator<=>(const Offset &) const = default;
};

Vector cross(const Vector &a, const Vector &b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

Int dot(const Vector &a, const Vector &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

std::vector<Offset> vertices(const std::vector<Inequality> &constraints) {
    std::vector<Offset> result;
    const int size = static_cast<int>(constraints.size());
    for (int i = 0; i < size - 2; ++i)
        for (int j = i + 1; j < size - 1; ++j) {
            const auto [a, l] = constraints[i];
            const auto [b, m] = constraints[j];
            const auto ab = cross(a, b);
            for (int k = j + 1; k < size; ++k) {
                const auto [c, n] = constraints[k];
                Int det = dot(ab, c);
                if (!det)
                    continue;
                const Vector bc = cross(b, c), ca = cross(c, a);
                Vector x{};
                for (int p = 0; p < 3; ++p)
                    x[p] = l * bc[p] + m * ca[p] + n * ab[p];
                if (det < 0) {
                    det = -det;
                    for (auto &p : x)
                        p = -p;
                }
                bool feasible = true;
                for (const auto &row : constraints)
                    if (Wide(row.a[0]) * x[0] + Wide(row.a[1]) * x[1] + Wide(row.a[2]) * x[2] >
                        Wide(row.b) * det) {
                        feasible = false;
                        break;
                    }
                if (feasible) {
                    Int divisor = det;
                    for (Int component : x)
                        divisor = std::gcd(divisor, component);
                    for (Int &component : x)
                        component /= divisor;
                    result.push_back({x, det / divisor});
                }
            }
        }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

template <Bound bound>
Int robust_volume(const Normals &normal, const Offset &offset, int mask, bool precise, int active = 7) {
    Int result = 0;
    for (int octant = 0; octant < (precise ? 8 : 1); ++octant) {
        Bounds bounds = cube;
        if (precise)
            for (int j = 0; j < 3; ++j)
                bounds[j] = octant & (1 << j) ? Interval{0, G / 2} : Interval{-G / 2, 0};
        std::array<Vector, 3> rows{};
        Vector h{};
        for (int i = 0; i < 3; ++i) {
            if (!(active & (1 << i))) {
                h[i] = -S;
                continue;
            }
            const Int sign = mask & (1 << i) ? 1 : -1;
            const Wide threshold = sign * offset.numerator[i];
            h[i] = bound == Bound::Lower ? ceil_div(threshold, offset.denominator)
                                         : floor_div(threshold, offset.denominator);
            Int radius = 0;
            for (int j = 0; j < 3; ++j) {
                const auto [lo, hi] = normal[i][j];
                if (precise) {
                    bool use_lo = (bool(octant & (1 << j)) == (sign > 0));
                    if constexpr (bound == Bound::Upper)
                        use_lo = !use_lo;
                    rows[i][j] = sign * (use_lo ? lo : hi);
                } else {
                    const Int mid = floor_div(lo + hi, 2);
                    rows[i][j] = sign * mid;
                    radius += std::max(mid - lo, hi - mid);
                }
            }
            if (!precise)
                h[i] += (bound == Bound::Lower ? 1 : -1) * ((radius + 1) / 2);
            if (h[i] < -S || h[i] > S)
                return bound == Bound::Lower ? 0 : V;
        }
        result += seven_integer::volume<bound>(rows, h, bounds);
    }
    return result;
}

constexpr std::array<Vector, 13> directions = {{{1, 0, 0},
                                                {0, 1, 0},
                                                {0, 0, 1},
                                                {1, 1, 0},
                                                {1, -1, 0},
                                                {1, 0, 1},
                                                {1, 0, -1},
                                                {0, 1, 1},
                                                {0, 1, -1},
                                                {1, 1, 1},
                                                {1, 1, -1},
                                                {1, -1, 1},
                                                {1, -1, -1}}};

Int central_cone_upper(const Normals &normals, const Interval &det, Int r_hi) {
    Int upper = WG / 2;
    const Wide full_det = Wide(det[1]) * r_hi;
    if (det[0])
        upper = std::min(upper, ceil_div(Wide(r_hi) * r_hi * WG, Wide(48) * det[0]));
    for (const Vector &e : directions) {
        std::array<Interval, 3> minors{};
        for (int i = 0; i < 3; ++i)
            for (int axis = 0; axis < 3; ++axis)
                if (e[axis]) {
                    auto m = minor_at(normals, (i + 1) % 3, (i + 2) % 3, axis);
                    minors[i] = add(minors[i], e[axis] > 0 ? m : negate(m));
                }
        Wide product = 1;
        if (minors[0][0] > 0 && minors[1][0] > 0 && minors[2][0] > 0)
            for (const auto m : minors)
                product *= m[0];
        else if (minors[0][1] < 0 && minors[1][1] < 0 && minors[2][1] < 0)
            for (const auto m : minors)
                product *= -m[1];
        else
            continue;
        const Int l1 = std::abs(e[0]) + std::abs(e[1]) + std::abs(e[2]);
        const Wide den = 48 * product;
        const Wide num = full_det * full_det * l1 * l1 * l1 * WG;
        if (num < Wide(upper) * den)
            upper = ceil_div(num, den);
    }
    return upper;
}

enum class Result {
    Open,
    Outside,
    Cone,
    Norm,
    Cap,
    Missing,
    Strip,
    Offset,
    CapSum,
    Pair,
    Central,
    Cell,
    Count,
};

Result check(Parameters &box, Chart chart, bool precise) {
    Normals normals;
    Interval vertical;
    if (!unpack(box, chart, normals, vertical))
        return Result::Outside;
    auto [r_lo, r_hi] = magnitude(vertical);
    const Interval det = magnitude(minor_at(normals, 0, 1, 2));
    const Int cone_threshold = floor_div(Wide(48) * det[0] * w_min_numerator, w_min_denominator);
    r_lo = std::max(r_lo, integer_sqrt(cone_threshold));
    if (r_hi < r_lo || Wide(r_hi) * r_hi * w_min_denominator < Wide(48) * det[0] * w_min_numerator)
        return Result::Cone;

    std::array<Vector, 3> min_weights{}, max_weights{};
    Vector min_offset{}, max_offset{}, min_maximum{}, norm_upper{}, support_twice{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const auto [lo, hi] = magnitude(normals[i][j]);
            min_weights[i][j] = lo;
            max_weights[i][j] = hi;
        }
        const Int max_norm = max_component(max_weights[i]);
        // A row this small forces the sum of the other two signed forms below
        // h_i+||n_i||_1/2 < (7/50*2+3/2)/50 < 1/25. One of those two normals
        // has a coefficient of magnitude 1, so its pair overlap is below 1/25.
        if (max_norm * 50 < S)
            return Result::Norm;
        min_maximum[i] = max_component(min_weights[i]);
        norm_upper[i] = sqrt_upper(dot(max_weights[i], max_weights[i]));
        support_twice[i] = max_weights[i][0] + max_weights[i][1] + max_weights[i][2];
        min_offset[i] = std::max(min_maximum[i] / 50, quantile(min_weights[i], 12, 25)[0]);
        max_offset[i] = std::min(ceil_div(7 * norm_upper[i], 50), quantile(max_weights[i], 9, 25)[1]);
        if (min_offset[i] > max_offset[i])
            return Result::Cap;
    }
    Vector sorted_norms = norm_upper;
    std::sort(sorted_norms.begin(), sorted_norms.end());
    r_hi = std::min(r_hi, ceil_div(7 * (sorted_norms[1] + sorted_norms[2]) + sorted_norms[0], 25));
    r_hi = std::min(r_hi, 2 * (max_offset[0] + max_offset[1] + max_offset[2]));
    if (r_hi < r_lo)
        return Result::Missing;
    if (det[0] > 0) {
        // The same offsets that exclude +++ obey r <= 2T <= 2N_max*Delta.
        // Combining this with the two cone bounds gives a quadratic in |r|.
        const Int longest = sorted_norms[2];
        const auto discriminant = [&](Int r) {
            return Wide(2000) * longest * r * r - Wide(48) * det[0] * (10 * longest + 99 * r);
        };
        if (discriminant(r_hi) < 0)
            return Result::Cone;
        if (discriminant(r_lo) < 0) {
            Int hi = r_hi;
            while (hi - r_lo > 1) {
                const Int midpoint = (r_lo + hi) / 2;
                if (discriminant(midpoint) < 0)
                    r_lo = midpoint;
                else
                    hi = midpoint;
            }
        }
    }
    // q_i = s_i + d_ij + d_ik, with s_i <= 4/25. Tighten h_i using
    // that cap upper bound, then reevaluate the two-halfspace strips.
    Int sum_upper_pairs = 3 * ceil_div(4 * V, 25);
    for (int iter = 0; iter < 3; ++iter) {
        const Vector upper = pair_strip_upper(normals, min_offset, det[1]);
        sum_upper_pairs = upper[0] + upper[1] + upper[2];
        if (*std::min_element(upper.begin(), upper.end()) * 25 < V)
            return Result::Strip;
        bool changed = false;
        for (int i = 0; i < 3; ++i) {
            const Int bound = ceil_div(4 * V, 25) + upper[(i + 1) % 3] + upper[(i + 2) % 3];
            if (bound * 25 < 9 * V)
                return Result::Strip;
            if (bound * 25 >= 12 * V)
                continue;
            Int lo = min_offset[i], hi = max_offset[i];
            const auto at_hi = cap(min_weights[i], hi);
            if (at_hi.num * V > at_hi.den * bound)
                return Result::Strip;
            while (hi - lo > 1) {
                const Int midpoint = (lo + hi) / 2;
                const auto q = cap(min_weights[i], midpoint);
                if (q.num * V > q.den * bound)
                    lo = midpoint;
                else
                    hi = midpoint;
            }
            changed |= lo > min_offset[i];
            min_offset[i] = lo;
        }
        if (!changed)
            break;
    }
    // The pair sum is at least 51/100-Delta/2. Smaller upper bounds on
    // the pairs force larger Delta, which forces a larger centered cone.
    const Int deficit_lower =
        std::max({3 * V, 51 * V - 100 * sum_upper_pairs, floor_div(Wide(25) * V * r_lo, sorted_norms[2])});
    const Int w_lower_num = 250 * V + 99 * deficit_lower;
    const Int w_lower_den = 50000 * V;
    r_lo = std::max(r_lo, integer_sqrt(floor_div(Wide(48) * det[0] * w_lower_num, w_lower_den)));
    if (r_hi < r_lo || Wide(r_hi) * r_hi * w_lower_den < Wide(48) * det[0] * w_lower_num)
        return Result::Cone;
    Int cone_upper = central_cone_upper(normals, det, r_hi);
    if (precise && Wide(cone_upper) * w_lower_den >= Wide(WG) * w_lower_num)
        cone_upper = std::min(
            cone_upper, ceil_div(Wide(robust_volume<Bound::Upper>(normals, {Vector{}, 1}, 7, true)) * WG, V));
    if (Wide(cone_upper) * w_lower_den < Wide(WG) * w_lower_num)
        return Result::Cone;

    std::vector<Inequality> constraints;
    for (int i = 0; i < 3; ++i) {
        Vector axis{};
        axis[i] = 1;
        constraints.push_back({axis, max_offset[i]});
        axis[i] = -1;
        constraints.push_back({axis, -min_offset[i]});
    }
    constraints.push_back({{-1, -1, -1}, -r_lo / 2});
    Vector inverse{};
    for (int i = 0; i < 3; ++i)
        inverse[i] = R * S / norm_upper[i];
    const Int deficit_rhs = std::min(ceil_div(Wide(3) * R * S, 10),
                                     ceil_div((Wide(200) * cone_upper - WG) * 5 * R * S, Wide(99) * WG));
    constraints.push_back({inverse, deficit_rhs});
    for (int k = 0; k < 3; ++k) {
        const int i = (k + 1) % 3, j = (k + 2) % 3;
        Vector row = {1, 1, 1};
        row[k] = 0;
        const Int bound =
            ceil_div(Wide(25) * (r_hi + support_twice[k]) - 2 * std::max(min_maximum[i], min_maximum[j]), 50);
        if (min_offset[i] + min_offset[j] > bound)
            return Result::Offset;
        constraints.push_back({row, bound});
    }

    Int pair_upper_sum = 0;
    const Offset lows = {min_offset, 1};
    const Offset highs = {max_offset, 1};
    for (int k = 0; k < 3; ++k) {
        const int mask = 7 ^ (1 << k);
        const Int upper = robust_volume<Bound::Upper>(normals, lows, mask, false, mask);
        if (25 * upper < V)
            return Result::Pair;
        const Int lower = robust_volume<Bound::Lower>(normals, highs, mask, false, mask);
        if (25 * lower > 4 * V)
            return Result::Pair;
        pair_upper_sum += upper;
    }
    if (25 * pair_upper_sum < 9 * V)
        return Result::Pair;

    const auto offsets = vertices(constraints);
    if (offsets.empty())
        return Result::Offset;
    bool cap_sum_possible = false;
    for (const auto &h : offsets) {
        Int upper = 0;
        for (int i = 0; i < 3; ++i) {
            const auto q = cap(max_weights[i], floor_div(h.numerator[i], h.denominator));
            upper += ceil_div(q.num * V, q.den);
        }
        // The cap sum is at least 6/5, and the bound on w also implies
        // sum(q_i) >= 307/198 - (1000/99)*w. Maximize the convex cap sum
        // over the relaxed offset polytope by evaluating its vertices.
        if (5 * upper >= 6 * V && 198 * WG * upper + 2000 * cone_upper * V >= 307 * WG * V) {
            cap_sum_possible = true;
            break;
        }
    }
    if (!cap_sum_possible)
        return Result::CapSum;
    if (det[0] > 0) {
        bool central_possible = false;
        for (const auto &h : offsets) {
            const Int T = h.numerator[0] + h.numerator[1] + h.numerator[2];
            const Wide d = h.denominator;
            // The all-negative slice is contained in its full horizontal triangle.
            const Wide upper_num = Wide(12) * T * T + Wide(r_hi) * r_hi * d * d;
            const Wide upper_den = Wide(24) * det[0] * d * d;
            // Its mass is at least 1/100 + (sum of cap deficits)/2.
            const Int upper = ceil_div(upper_num * V, upper_den);
            const Wide scale = Wide(R) * S * d;
            const Wide offset_sum = Wide(inverse[0]) * h.numerator[0] + Wide(inverse[1]) * h.numerator[1] +
                                    Wide(inverse[2]) * h.numerator[2];
            if (Wide(100) * upper * scale >= Wide(V) * (scale + 50 * offset_sum)) {
                central_possible = true;
                break;
            }
        }
        if (!central_possible)
            return Result::Central;
    }
    for (bool sharp : {false, true}) {
        if (sharp && !precise)
            break;
        for (int mask : {1, 2, 4, 0, 3, 5, 6}) {
            bool works = true;
            for (const auto &offset : offsets)
                if (25 * robust_volume<Bound::Lower>(normals, offset, mask, sharp) <= 4 * V) {
                    works = false;
                    break;
                }
            if (works)
                return Result::Cell;
        }
    }
    return Result::Open;
}

struct Node {
    Parameters box;
    int depth;
};

Parameters domain(Chart chart) {
    if (chart == Chart::Horizontal)
        return {{{0, S}, {-S / 2, 0}, {-S, S}, {-S, S}, {-S, S}, {-S, S}}};
    return {{{0, S}, {0, S}, {-S, S}, {-S, S}, {-S, S}, {-S, S}}};
}

Parameters near_antipodal_domain(bool wide) {
    if (wide)
        return {{{0, S}, {-S / 4, 0}, {-S / 2, S / 2}, {-S, -S / 2}, {-S, S}, {S / 2, S}}};
    return {{{0, S}, {-S / 8, 0}, {-S / 4, S / 4}, {-S, -3 * S / 4}, {-S, S}, {3 * S / 4, S}}};
}

int main(int argc, char **argv) {
    const std::string_view selection = argc > 1 ? argv[1] : "--wide-antipodal";
    if (selection != "--near-antipodal" && selection != "--wide-antipodal" && selection != "--horizontal" &&
        selection != "--vertical") {
        std::cerr << "usage: seven_scaled_charts_certificate "
                     "[--near-antipodal|--wide-antipodal|--horizontal|--vertical] [depth] [node-limit]\n";
        return 2;
    }
    const Chart chart = selection == "--vertical" ? Chart::Vertical : Chart::Horizontal;
    const int depth = argc > 2 ? std::atoi(argv[2]) : 30;
    const Int node_limit = argc > 3 ? std::atoll(argv[3]) : 1000000;
    Parameters initial = domain(chart);
    if (selection == "--near-antipodal" || selection == "--wide-antipodal")
        initial = near_antipodal_domain(selection == "--wide-antipodal");
    std::vector<Node> stack = {{initial, 0}};
    std::array<Int, static_cast<std::size_t>(Result::Count)> stats{};
    Int visited = 0, open = 0;
    while (!stack.empty() && visited < node_limit) {
        auto [box, level] = stack.back();
        stack.pop_back();
        ++visited;
        const auto result = check(box, chart, level >= 12);
        if (result != Result::Open) {
            ++stats[static_cast<std::size_t>(result)];
            continue;
        }
        int split = -1;
        Int score = 0;
        for (int i = 0; i < 6; ++i) {
            const auto width = box[i][1] - box[i][0];
            const Int priority = width * ((chart == Chart::Horizontal ? i < 3 : i < 4) ? 2 : 1);
            if (width > 1 && priority > score) {
                score = priority;
                split = i;
            }
        }
        if (level >= depth || split < 0) {
            if (++open == 1) {
                std::cout << "first open box (scale " << S << "):";
                for (auto pair : box)
                    std::cout << " [" << pair[0] << ',' << pair[1] << ']';
                std::cout << '\n';
            }
            continue;
        }
        auto next = box;
        const Int midpoint = floor_div(box[split][0] + box[split][1], 2);
        box[split][1] = midpoint;
        next[split][0] = midpoint;
        stack.push_back({next, level + 1});
        stack.push_back({box, level + 1});
    }
    const bool complete = open == 0 && stack.empty();
    std::cout << selection << ": " << (complete ? "certified" : "incomplete") << " visited=" << visited
              << " unresolved=" << open << " pending=" << stack.size() << '\n';
    constexpr std::array<std::string_view, static_cast<std::size_t>(Result::Count)> names = {
        "open",           "outside", "cone",    "small-normal", "cap",     "missing",
        "opposite-strip", "offset",  "cap-sum", "pair",         "central", "cell",
    };
    for (std::size_t i = 1; i < stats.size(); ++i)
        std::cout << names[i] << '=' << stats[i] << ' ';
    std::cout << '\n';
    return complete ? 0 : 1;
}
