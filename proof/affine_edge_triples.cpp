// Three affine cuts with two-vertex sides. Volumes are fractions of the unit cube.
#include "seven_integer_volume.hpp"
#include <bit>
#include <iostream>

using seven_integer::Bound;
using seven_integer::Bounds;
using seven_integer::Int;
using seven_integer::Vector;
using seven_integer::space_grid;
using seven_integer::volume_grid;
using std::array;

namespace {
constexpr Int grid = 4096;
constexpr Bounds cube{{{-space_grid / 2, space_grid / 2},
                       {-space_grid / 2, space_grid / 2},
                       {-space_grid / 2, space_grid / 2}}};

struct EdgeCut {
    int along, first, second, endpoint;
    // p0, p1 and w, in units of 1/grid.
    Vector lo{0, 0, grid}, hi{grid, grid, 2 * grid};
};

EdgeCut edge_chart(int split, bool swap) {
    const unsigned bits = std::popcount(unsigned(split)) == 2 ? split : 255 ^ split;
    assert(std::popcount(bits) == 2);
    const int a = std::countr_zero(bits), b = std::countr_zero(bits & (bits - 1));
    assert(std::has_single_bit(unsigned(a ^ b)));
    const int along = std::countr_zero(unsigned(a ^ b));
    const int first = (along + 1 + swap) % 3;
    return {along, first, 3 - along - first, a};
}

// The cap volume (p0^2 + p0*p1 + p1^2)/(6w) must be at least 9/25.
bool restrict_cap(EdgeCut &cut) {
    auto &lo = cut.lo;
    auto &hi = cut.hi;
    auto numerator = [](Int a, Int b) { return 25 * (a * a + a * b + b * b); };
    const Int needed = 54 * grid * lo[2];
    const Int largest = numerator(hi[0], hi[1]);
    if (largest < needed)
        return false;
    hi[2] = std::min(hi[2], (largest + 54 * grid - 1) / (54 * grid));
    for (int j = 0; j < 2; ++j) {
        if (numerator(lo[j], hi[1 - j]) >= needed)
            continue;
        Int a = lo[j], b = hi[j];
        while (b - a > 1) {
            const Int m = (a + b) / 2;
            if (numerator(m, hi[1 - j]) < needed)
                a = m;
            else
                b = m;
        }
        lo[j] = a; // The true lower endpoint can lie between a and b.
    }
    return true;
}

struct Plane {
    Vector normal;
    Int height;
};

Plane guaranteed_side(const EdgeCut &cut, bool inside) {
    const Int p0 = inside ? cut.lo[0] : cut.hi[0];
    const Int p1 = inside ? cut.lo[1] : cut.hi[1];
    const Int w = inside ? cut.hi[2] : cut.lo[2];
    // In [0,1]^3: distance_1 + w*distance_2 <= (1-x)*p0 + x*p1.
    // Translate to the centered cube and clear the factor 1/2.
    const int sign = inside ? 1 : -1;
    Vector n{};
    n[cut.along] = sign * 2 * (p1 - p0);
    n[cut.first] = sign * ((cut.endpoint >> cut.first) & 1 ? 2 : -2) * grid;
    n[cut.second] = sign * ((cut.endpoint >> cut.second) & 1 ? 2 : -2) * w;
    return {n, sign * (grid + w - p0 - p1)};
}

using Box = array<EdgeCut, 3>;
struct Stats {
    Int nodes = 0, leaves = 0;
    int depth = 0;
};

bool cover(Box box, Stats &stats, int depth = 0) {
    ++stats.nodes;
    stats.depth = std::max(stats.depth, depth);
    for (auto &cut : box) {
        if (!restrict_cap(cut)) {
            ++stats.leaves;
            return true;
        }
    }
    for (int mask = 0; mask < 8; ++mask) {
        array<Vector, 3> normals;
        Vector heights;
        for (int i = 0; i < 3; ++i) {
            const auto p = guaranteed_side(box[i], (mask >> i) & 1);
            normals[i] = p.normal;
            heights[i] = p.height;
        }
        const Int lower = seven_integer::volume<Bound::Lower>(normals, heights, cube);
        if (25 * lower > 4 * volume_grid) {
            ++stats.leaves;
            return true;
        }
    }
    Int width = 0;
    int row = 0, coordinate = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const Int length = box[i].hi[j] - box[i].lo[j];
            if (length > width) {
                width = length;
                row = i;
                coordinate = j;
            }
        }
    }
    if (width <= 1 || depth == 48)
        return false;
    const Int mid = (box[row].hi[coordinate] + box[row].lo[coordinate]) / 2;
    Box other = box;
    box[row].hi[coordinate] = other[row].lo[coordinate] = mid;
    return cover(box, stats, depth + 1) && cover(other, stats, depth + 1);
}
} // namespace

int main() {
    const std::array<std::array<int, 3>, 7> triples{
        {{3, 5, 10}, {3, 5, 17}, {3, 5, 34}, {3, 5, 48}, {3, 5, 119}, {3, 12, 80}, {3, 68, 95}}};
    for (const auto &triple : triples) {
        Stats stats;
        for (int chart = 0; chart < 8; ++chart) {
            Box box;
            for (int i = 0; i < 3; ++i)
                box[i] = edge_chart(triple[i], (chart >> i) & 1);
            if (!cover(box, stats)) {
                std::cerr << "Unresolved triple " << triple[0] << ' ' << triple[1] << ' ' << triple[2]
                          << ", chart " << chart << '\n';
                return 1;
            }
        }
        assert(2 * stats.leaves == stats.nodes + 8);
        std::cout << triple[0] << ' ' << triple[1] << ' ' << triple[2] << "   " << stats.nodes << ' '
                  << stats.leaves << ' ' << stats.depth << '\n';
    }
    return 0;
}
