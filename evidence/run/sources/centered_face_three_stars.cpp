// Centered face and star direction covers. With no arguments check the original XBCD case.
// Build: c++ -std=c++20 -O3 centered_face_three_stars.cpp -o /tmp/cube-face-stars
// A nonzero exit code means the cover did not finish.
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {
using Int = std::int64_t;
using Wide = __int128_t;
using Vector = std::array<Int, 3>;
using Triangle = std::array<Vector, 3>;

constexpr Int area_grid = Int{1} << 25;
constexpr Int unit = 24 * area_grid; // grid denominator for a disagreement probability
constexpr Int cap_num = 2, cap_den = 25;
constexpr int max_depth = 32;
constexpr std::array<std::array<int, 2>, 6> edges{{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}};

Int floor_ratio(Wide a, Wide b) {
    assert(b > 0);
    return Int(a / b - (a % b < 0));
}

Vector cross(Vector a, Vector b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}

struct Vertex {
    Vector point; // (X,Y,W), W > 0, represents (X/W,Y/W)
    int incoming; // supporting line of the edge arriving here
};

// Clip one face of [-1,1]^3 and round twice its area down to the grid.
// Every vertex is the intersection of two original lines.
Int face_area_lower(const std::array<Vector, 6>& normals, int face) {
    std::array<Vector, 10> lines{{{0,1,1},{-1,0,1},{0,-1,1},{1,0,1}}};
    std::array<Vertex, 12> polygon{{{{-1,-1,1},3},{{1,-1,1},0},{{1,1,1},1},{{-1,1,1},2}}};
    int count = 4;
    const int axis = face % 3;
    for (int i = 0; i < 6; ++i) {
        const auto& n = normals[i];
        const int line_id = 4 + i;
        const Vector line{n[(axis+1)%3], n[(axis+2)%3], face < 3 ? n[axis] : -n[axis]};
        lines[line_id] = line;
        std::array<Int, 12> values;
        bool inside = false, outside = false;
        for (int j = 0; j < count; ++j) {
            auto [x, y, w] = polygon[j].point;
            values[j] = line[0]*x + line[1]*y + line[2]*w;
            inside |= values[j] > 0;
            outside |= values[j] < 0;
        }
        if (!outside) continue;
        if (!inside) return 0;
        std::array<Vertex, 12> clipped;
        int size = 0;
        for (int j = 0; j < count; ++j) {
            const Int before = values[j == 0 ? count-1 : j-1], after = values[j];
            const int old_line = polygon[j].incoming;
            if ((before < 0 && after > 0) || (before > 0 && after < 0)) {
                Vector p = cross(lines[old_line], line);
                if (p[2] < 0) for (Int& x : p) x = -x;
                assert(p[2] > 0);
                clipped[size++] = {p, before < 0 ? line_id : old_line};
            }
            if (after >= 0) {
                // If we enter at an existing vertex, the preceding edge is the new cut.
                const int incoming = before < 0 && after == 0 ? line_id : old_line;
                clipped[size++] = {polygon[j].point, incoming};
            }
        }
        assert(size <= line_id+1);
        count = size;
        std::copy_n(clipped.begin(), count, polygon.begin());
    }
    Int area = 0;
    for (int j = 0; j < count; ++j) {
        const Vector& p = polygon[j].point;
        const Vector& q = polygon[(j+1)%count].point;
        area += floor_ratio((Wide(p[0])*q[1] - Wide(p[1])*q[0])*area_grid,
                            Wide(p[2])*q[2]);
    }
    return std::max(Int{0}, area);
}

struct Region {
    Triangle vertices;
    std::uint32_t id;
};

Region child(const Region& parent, int which) {
    Region result{};
    result.id = 4*parent.id + which;
    // All new vertices are represented at double the parent's scale.
    for (int i = 0; i < 3; ++i) {
        const Vector& a = parent.vertices[which == 3 ? (i+1)%3 : i];
        const Vector& b = parent.vertices[which == 3 ? (i+2)%3 : which];
        for (int j = 0; j < 3; ++j) result.vertices[i][j] = a[j] + b[j];
    }
    return result;
}

struct Range { Int lo, hi; };
using Box = std::array<Region, 4>;

Range disagreement_range(const Triangle& a, const Triangle& b) {
    std::array<Vector, 6> normals{a[0], a[1], a[2], b[0], b[1], b[2]};
    Int agreement = 0, disagreement = 0;
    for (int f = 0; f < 6; ++f) agreement += face_area_lower(normals, f);
    for (int j = 3; j < 6; ++j) for (Int& x : normals[j]) x = -x;
    for (int f = 0; f < 6; ++f) disagreement += face_area_lower(normals, f);
    return {disagreement, unit-agreement};
}

std::array<std::unordered_map<std::uint64_t, Range>, 6> pair_cache;

Range pair_range(const Box& box, int edge) {
    const auto [i, j] = edges[edge];
    const auto key = (std::uint64_t(box[i].id) << 32) | box[j].id;
    auto& cache = pair_cache[edge];
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second;

    Range enclosure{unit, 0};
    for (int k = 0; k < 4; ++k) {
        const auto a = child(box[i], k);
        for (int l = 0; l < 4; ++l) {
            const auto b = child(box[j], l);
            const auto bounds = disagreement_range(a.vertices, b.vertices);
            enclosure.lo = std::min(enclosure.lo, bounds.lo);
            enclosure.hi = std::max(enclosure.hi, bounds.hi);
        }
    }
    assert(0 <= enclosure.lo && enclosure.lo <= enclosure.hi && enclosure.hi <= unit);
    cache.emplace(key, enclosure);
    return enclosure;
}

struct Constraint {
    int cells = 0, constant = 0;
    std::array<int, 6> coefficients{};
};

std::vector<Constraint> make_constraints() {
    std::array<Constraint, 7> cells{};
    for (int i = 0; i < 7; ++i) {
        // Each row represents four times one cell volume.
        const bool singleton = i < 4;
        const int pair_mask = singleton ? 0 : 1 | (1 << (i-3));
        cells[i].constant = singleton ? 2 : -2;
        for (int e = 0; e < 6; ++e) {
            const auto [a, b] = edges[e];
            cells[i].coefficients[e] = singleton ? -int(i != a && i != b)
                                                 : ((pair_mask >> a) ^ (pair_mask >> b)) & 1;
        }
    }
    std::vector<Constraint> result;
    for (int k = 1; k < 7; ++k) {
        for (unsigned selected = 1; selected < 127; ++selected) {
            if (std::popcount(selected) != k) continue;
            Constraint joint;
            joint.cells = k;
            for (int i = 0; i < 7; ++i) if ((selected >> i) & 1) {
                joint.constant += cells[i].constant;
                for (int e = 0; e < 6; ++e) joint.coefficients[e] += cells[i].coefficients[e];
            }
            result.push_back(joint);
        }
    }
    return result;
}

const auto constraints = make_constraints();

// Orient the four normals so the all-positive cell has zero volume.
bool possible(int flips, const std::array<Range, 6>& original) {
    std::array<Range, 6> d;
    constexpr Int lower = unit*(cap_den-6*cap_num)/cap_den;
    constexpr Int upper = (unit*8*cap_num+cap_den-1)/cap_den;
    for (int e = 0; e < 6; ++e) {
        const auto [a, b] = edges[e];
        auto [lo, hi] = original[e];
        if (((flips >> a) ^ (flips >> b)) & 1) d[e] = {unit-hi, unit-lo};
        else d[e] = {lo, hi};
        d[e].lo = std::max(d[e].lo, lower);
        d[e].hi = std::min(d[e].hi, upper);
        if (d[e].lo > d[e].hi) return false;
    }
    for (const auto& test : constraints) {
        Int minimum = test.constant*unit;
        for (int e = 0; e < 6; ++e) {
            const int c = test.coefficients[e];
            minimum += c * (c >= 0 ? d[e].lo : d[e].hi);
        }
        if (minimum*cap_den > 4*test.cells*cap_num*unit) return false;
    }
    return true;
}

std::uint64_t nodes = 0, leaves = 0;
int deepest = 0;

bool cover(Box& box, int depth, unsigned orientations) {
    ++nodes;
    deepest = std::max(deepest, depth);
    std::array<Range, 6> bounds;
    for (int e = 0; e < 6; ++e) bounds[e] = pair_range(box, e);
    for (int m = 0; m < 8; ++m)
        if ((orientations & (1u << m)) && !possible(m, bounds)) orientations &= ~(1u << m);
    if (!orientations) {
        ++leaves;
        return true;
    }
    if (depth == max_depth) {
        std::cerr << "Unresolved at depth " << depth << ", orientations " << orientations << '\n';
        return false;
    }
    const int i = depth % 4;
    const Region parent = box[i];
    for (int k = 0; k < 4; ++k) {
        box[i] = child(parent, k);
        if (!cover(box, depth+1, orientations)) return false;
    }
    box[i] = parent;
    return true;
}
// A, B, C, D have star centers (1,+,+), (1,+,-), (1,-,+), (1,-,-).
Region direction_region(char label, int sign) {
    if (label == 'X' || label == 'Y' || label == 'Z') {
        const int axis = label == 'X' ? 0 : label == 'Y' ? 1 : 2;
        Triangle face{};
        for (auto& vertex : face) vertex[axis] = 1;
        face[0][(axis+1)%3] = 1;
        face[1][(axis+1)%3] = -1;
        face[2][(axis+2)%3] = sign;
        return {face, sign < 0 ? 1u : 2u};
    }
    Triangle star{{{1,1,0},{1,0,1},{0,1,1}}};
    for (auto& v : star) {
        if (label == 'C' || label == 'D') v[1] = -v[1];
        if (label == 'B' || label == 'D') v[2] = -v[2];
    }
    return {star, 1};
}
} // namespace

int main(int argc, char** argv) {
    const std::string_view labels = argc == 1 ? "XBCD" : argv[1];
    if (argc > 2 || labels.size() != 4 || labels.find_first_not_of("XYZABCD") != labels.npos) {
        std::cerr << "Usage: centered_face_three_stars [four labels from XYZABCD]\n";
        return 2;
    }
    const int roots = 1 << std::count_if(labels.begin(), labels.end(), [](char c) { return c >= 'X'; });
    for (int root = 0; root < roots; ++root) {
        Box box;
        int bits = root;
        for (int i = 0; i < 4; ++i) {
            box[i] = direction_region(labels[i], bits & 1 ? 1 : -1);
            if (labels[i] >= 'X') bits >>= 1;
        }
        if (!cover(box, 0, 255)) return 1;
    }
    assert(3*nodes + roots == 4*leaves);
    std::cout << labels << ": cap " << cap_num << '/' << cap_den << " excluded. Roots " << roots
              << ", nodes " << nodes << ", leaves " << leaves << ", maximum depth " << deepest
              << ", unresolved 0\n";
    return 0;
}
