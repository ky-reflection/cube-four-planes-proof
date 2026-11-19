// Two affine cuts inducing the same corner star on the cube vertices.
#include "seven_integer_volume.hpp"
#include <iostream>

using seven_integer::Bound;
using seven_integer::Int;
using seven_integer::Vector;
using seven_integer::space_grid;
using seven_integer::volume_grid;

namespace {
constexpr Int grid = 4096;
constexpr seven_integer::Bounds cube{{{-space_grid / 2, space_grid / 2},
                                      {-space_grid / 2, space_grid / 2},
                                      {-space_grid / 2, space_grid / 2}}};
struct Row {
    Vector lo{0, 0, 0}, hi{grid, grid, grid};
};
using Box = std::array<Row, 2>;
Int nodes = 0, leaves = 0;
int deepest = 0;

bool cover(Box box, int depth) {
    ++nodes;
    deepest = std::max(deepest, depth);
    // Simultaneous coordinate permutation puts the first row in ascending order.
    auto &first = box[0];
    for (int j = 1; j < 3; ++j)
        first.lo[j] = std::max(first.lo[j], first.lo[j - 1]);
    for (int j = 1; j >= 0; --j)
        first.hi[j] = std::min(first.hi[j], first.hi[j + 1]);
    for (auto &row : box) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (i != j)
                    row.lo[i] = std::max(row.lo[i], grid - row.hi[j]);
            }
            if (row.lo[i] > row.hi[i]) {
                ++leaves;
                return true;
            }
        }
    }

    for (int mask : {0, 3, 1, 2}) {
        std::array<Vector, 3> normals{};
        Vector heights{}; // Third inequality is 0 >= 0.
        for (int i = 0; i < 2; ++i) {
            const bool inside = (mask >> i) & 1;
            const auto &a = inside ? box[i].hi : box[i].lo;
            const int sign = inside ? 1 : -1;
            for (int j = 0; j < 3; ++j)
                normals[i][j] = -2 * sign * a[j];
            heights[i] = sign * (a[0] + a[1] + a[2] - 2 * grid);
        }
        const Int lower = seven_integer::volume<Bound::Lower>(normals, heights, cube);
        if (25 * lower > 8 * volume_grid) {
            ++leaves;
            return true;
        }
    }
    Int width = 0;
    int row = 0, coordinate = 0;
    for (int i = 0; i < 2; ++i) {
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
    const Int mid = (box[row].lo[coordinate] + box[row].hi[coordinate]) / 2;
    Box other = box;
    box[row].hi[coordinate] = other[row].lo[coordinate] = mid;
    return cover(box, depth + 1) && cover(other, depth + 1);
}
} // namespace

int main() {
    if (!cover(Box{}, 0)) {
        std::cerr << "Unresolved star pair\n";
        return 1;
    }
    assert(2 * leaves == nodes + 1);
    std::cout << "Affine star pair: " << nodes << " nodes, " << leaves << " leaves, depth " << deepest
              << ", unresolved 0\n";
}
