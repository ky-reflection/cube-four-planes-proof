#ifndef CUBE_SECTION_SIMPLEX_HPP
#define CUBE_SECTION_SIMPLEX_HPP
#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>
namespace cube_section {
using Scalar = std::int64_t;
using Triple = std::array<Scalar, 3>;
constexpr int normals[11][3] = {{1, 0, 0}, {0, 1, 0},    {0, 0, 1},   {-1, 0, 0},  {0, -1, 0}, {0, 0, -1},
                                {1, 1, 1}, {-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
struct Basis {
    int ids[3], inv[3][3];
};
inline const auto bases = [] {
    std::vector<Basis> out;
    for (int i = 0; i < 11; ++i)
        for (int j = 0; j < i; ++j)
            for (int k = 0; k < j; ++k) {
                Basis b{{k, j, i}, {}};
                int cof[3][3]{};
                for (int a = 0; a < 3; ++a)
                    for (int n = 0; n < 3; ++n) {
                        auto x = normals[b.ids[(a + 1) % 3]], y = normals[b.ids[(a + 2) % 3]];
                        cof[a][n] = x[(n + 1) % 3] * y[(n + 2) % 3] - x[(n + 2) % 3] * y[(n + 1) % 3];
                    }
                int d = 0;
                for (int h = 0; h < 3; ++h)
                    d += cof[0][h] * normals[k][h];
                if (!d)
                    continue;
                if (4 % d)
                    throw std::logic_error("unexpected section determinant");
                for (int a = 0; a < 3; ++a)
                    for (int n = 0; n < 3; ++n)
                        b.inv[a][n] = 4 * cof[n][a] / d;
                out.push_back(b);
            }
    return out;
}();
struct Column {
    std::vector<Triple> vertices;
    void init(Scalar unit, Triple lo, Triple hi, Scalar rlo, Scalar rhi) {
        Scalar v[11] = {lo[0],          lo[1],           lo[2], -hi[0], -hi[1], -hi[2],
                        unit + 2 * rlo, -unit - 2 * rhi, -unit, -unit,  -unit};
        for (auto &b : bases) {
            Triple p{};
            for (int a = 0; a < 3; ++a)
                for (int n = 0; n < 3; ++n)
                    p[a] += b.inv[a][n] * v[b.ids[n]];
            bool valid = true;
            for (int k = 0; k < 3; ++k)
                valid &= 4 * lo[k] <= p[k] && p[k] <= 4 * hi[k];
            if (!valid)
                continue;
            Scalar sum = p[0] + p[1] + p[2];
            if (sum < 4 * v[6] || sum > -4 * v[7])
                continue;
            for (int k = 0; k < 3; ++k)
                valid &= 2 * p[k] >= sum - 4 * unit;
            if (valid && std::find(vertices.begin(), vertices.end(), p) == vertices.end())
                vertices.push_back(p);
        }
    }
    Triple lowest(std::array<int, 3> cost) const {
        if (vertices.empty())
            throw std::logic_error("empty section polytope");
        Scalar best = 0;
        const Triple *at = nullptr;
        for (auto &v : vertices) {
            Scalar f = 0;
            for (int i = 0; i < 3; ++i)
                f += cost[i] * v[i];
            if (!at || f < best) {
                best = f;
                at = &v;
            }
        }
        Triple z = *at;
        for (int i = 0; i < 3; ++i)
            z[i] = (z[i] + (cost[i] < 0 ? 3 : 0)) / 4;
        return z;
    }
};
} // namespace cube_section
#endif
