#include <numeric>
// Certified cover for affine vertex types from vertex_partition_census.py.
// Optional argument: per-type node limit, zero for unlimited.
#include "affine_face_integral.hpp"
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string_view>
#include <unordered_map>
#include <vector>

using namespace seven_integer;
using std::array;
namespace {
constexpr Int parameter_scale = 2048, cap_numerator = 79759, cap_denominator = 1000000;
constexpr Bounds cube{{{-space_grid / 2, space_grid / 2},
                       {-space_grid / 2, space_grid / 2},
                       {-space_grid / 2, space_grid / 2}}};
struct Plane {
    Vector n{};
    Int h = 0;
};

// The edge-length and Gram certificate gives t <= 8972/10^6 at this target.
constexpr Int tetra_cap = 8972, unit_grid = 32768;
constexpr Int half_min = cap_denominator - 7 * cap_numerator - tetra_cap;
struct Interval {
    Int lo = 0, hi = 0;
};
Interval negate(Interval a) { return {-a.hi, -a.lo}; }
Interval product(Interval a, Interval b) {
    auto v = std::minmax({a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi});
    return {v.first, v.second};
}
Int norm_floor(Int coord, Int square) {
    if (!square || !coord)
        return 0;
    const Int numerator = coord * coord * unit_grid * unit_grid;
    Int r = std::sqrt(double(numerator) / square);
    while ((r + 1) * (r + 1) * square <= numerator)
        ++r;
    while (r * r * square > numerator)
        --r;
    return r;
}

constexpr Int graph_scale = 4096;
struct Projection {
    array<Interval, 3> coefficient{};
    int sign = 0;
};
struct Row {
    array<Projection, 3> projection{};
    array<Interval, 3> unit;
    array<int, 4> lo{}, hi{};
    array<Plane, 2> side{};
    array<std::uint64_t, 2> boundary{};
    int width = 0, axis = 0, child = -1;
};
int population(int v) { return std::popcount(unsigned(v)); }
int first(int v) { return std::countr_zero(unsigned(v)); }
bool bit(int v, int k) { return (v >> k) & 1; }

enum class Method { two_edges, graph, marginals, tetrahedron, correlations };
struct Domain {
    enum class Kind { edge, graph, star };
    int bits, count;
    Kind kind;
    Method method;
    int k = 0, i = 0, j = 0, corner = 0;
    std::vector<Row> rows;
    explicit Domain(int split, Method m = Method::two_edges) : method(m) {
        bits = population(split) > 4 ? 255 ^ split : split;
        count = population(bits);
        assert(count >= 2 && count <= 4);
        int vertices[8], n = 0;
        for (int v = 0; v < 8; ++v)
            if (bit(bits, v))
                vertices[n++] = v;
        Row r;
        for (int t = 0; t < 4; ++t)
            r.hi[t] = parameter_scale;
        if (count == 2) {
            assert(std::has_single_bit(unsigned(vertices[0] ^ vertices[1])));
            kind = Kind::edge;
            corner = vertices[0];
            k = first(vertices[0] ^ vertices[1]);
            i = (k + 1) % 3;
            j = (k + 2) % 3;
            // Normalize one fixed transverse coefficient to one; w may lie on either
            // side of one.
            r.lo[2] = parameter_scale / 2;
            r.hi[2] = 2 * parameter_scale;
        } else {
            int same = 7;
            for (int t = 1; t < n; ++t)
                same &= ~(vertices[t] ^ vertices[0]);
            if (same) {
                assert(std::has_single_bit(unsigned(same)));
                kind = Kind::graph;
                k = first(same);
                i = (k + 1) % 3;
                j = (k + 2) % 3;
                corner = vertices[0];
                for (int t = 0; t < 4; ++t) {
                    int v = (bit(corner, k) << k) | ((t & 1) << i) | ((t >> 1) << j);
                    if (!bit(bits, v)) {
                        r.lo[t] = -parameter_scale;
                        r.hi[t] = 0;
                    }
                }
            } else {
                kind = Kind::star;
                for (int v = 0; v < 8; ++v)
                    if (bit(bits, v) && bit(bits, v ^ 1) && bit(bits, v ^ 2) && bit(bits, v ^ 4))
                        corner = v;
                assert(count == 4 && bit(bits, corner) && bit(bits, corner ^ 1) && bit(bits, corner ^ 2) &&
                       bit(bits, corner ^ 4));
            }
        }
        bool ok = prepare(r);
        if (!ok)
            std::abort();
        rows.push_back(r);
    }
    // Exact CDF of a weighted sum of independent unit uniforms.
    static bool star_cap_exceeds(const array<int, 4> &a, Int num, Int den) {
        Int divisor = 1;
        unsigned nonzero = 0;
        int degree = 0;
        for (int j = 0; j < 3; ++j)
            if (a[j]) {
                nonzero |= 1u << j;
                divisor *= a[j] * (++degree);
            }
        Int sum = 0;
        for (unsigned mask = nonzero;; mask = (mask - 1) & nonzero) {
            Int z = parameter_scale;
            for (int j = 0; j < 3; ++j)
                if ((mask >> j) & 1)
                    z -= a[j];
            if (z > 0) {
                Int v = 1;
                for (int j = 0; j < degree; ++j)
                    v *= z;
                sum += (population(mask) & 1) ? -v : v;
            }
            if (!mask)
                break;
        }
        return sum * den > num * divisor;
    }
    bool prepare(Row &r) const {
        if (kind == Kind::edge) {
            const Int den = method == Method::tetrahedron ? cap_denominator : 25;
            const Int required = method == Method::tetrahedron ? 6 * half_min : 54;
            auto num = [den](Int a, Int b) { return den * (a * a + a * b + b * b); };
            for (int pass = 0; pass < 2; ++pass) {
                for (int t = 0; t < 2; ++t)
                    r.hi[t] = std::min(r.hi[t], r.hi[2]);
                Int best = num(r.hi[0], r.hi[1]), need = required * parameter_scale * r.lo[2];
                if (best < need)
                    return false;
                r.hi[2] = std::min<Int>(r.hi[2], (best + required * parameter_scale - 1) /
                                                     (required * parameter_scale));
                // w >= 18/25 follows from p0,p1 <= w and the cap volume bound.
                r.lo[2] = std::max<Int>(r.lo[2], required * parameter_scale / (3 * den));
                for (int t = 0; t < 2; ++t) {
                    if (num(r.lo[t], r.hi[1 - t]) >= need)
                        continue;
                    int a = r.lo[t], b = r.hi[t];
                    while (b - a > 1) {
                        int m = (a + b) / 2;
                        if (num(m, r.hi[1 - t]) < need)
                            a = m;
                        else
                            b = m;
                    }
                    r.lo[t] = a;
                }
            }
        } else if (kind == Kind::graph) {
            // The four heights of an affine graph have equal diagonal sums.
            if (method != Method::two_edges) {
                const bool tight = method == Method::tetrahedron;
                const int minsum =
                    count == 4 ? (tight ? 2 * half_min * parameter_scale / cap_denominator
                                        : 18 * parameter_scale / 25)
                               : (tight ? 17 * parameter_scale / 20
                                        : (method == Method::marginals ? 24 * parameter_scale / 35 : 0));
                const int maxsum = count == 3 ? parameter_scale
                                              : (tight ? (2 * (cap_denominator - half_min) * parameter_scale +
                                                          cap_denominator - 1) /
                                                             cap_denominator
                                                       : (32 * parameter_scale + 24) / 25);
                for (int j = 0; j < 4; ++j) {
                    r.lo[j] = std::max(r.lo[j], minsum - r.hi[j ^ 3]);
                    r.hi[j] = std::min(r.hi[j], maxsum - r.lo[j ^ 3]);
                }
            }
            for (int pass = 0; pass < 2; ++pass)
                for (int t = 0; t < 4; ++t) {
                    r.lo[t] = std::max(r.lo[t], r.lo[t ^ 1] + r.lo[t ^ 2] - r.hi[t ^ 3]);
                    r.hi[t] = std::min(r.hi[t], r.hi[t ^ 1] + r.hi[t ^ 2] - r.lo[t ^ 3]);
                }
        } else {
            // Each two star coefficients sum to at least one.
            for (int t = 0; t < 3; ++t)
                for (int h = 0; h < 3; ++h)
                    if (t != h)
                        r.lo[t] = std::max<Int>(r.lo[t], parameter_scale - r.hi[h]);
        }
        if (method >= Method::marginals && kind == Kind::star) {
            const bool tight = method == Method::tetrahedron;
            const Int den = tight ? cap_denominator : 25, cap_lo = tight ? half_min - 1 : 9,
                      cap_hi = tight ? cap_denominator - half_min : 16;
            for (int j = 0; j < 3; ++j) {
                if (r.lo[j] > r.hi[j])
                    return false;
                auto small = r.hi;
                small[j] = r.lo[j];
                if (star_cap_exceeds(small, cap_hi, den)) {
                    int a = r.lo[j], b = r.hi[j];
                    if (star_cap_exceeds(r.hi, cap_hi, den))
                        return false;
                    while (b - a > 1) {
                        int mid = (a + b) / 2;
                        small[j] = mid;
                        if (star_cap_exceeds(small, cap_hi, den))
                            a = mid;
                        else
                            b = mid;
                    }
                    r.lo[j] = a;
                }
                auto large = r.lo;
                large[j] = r.hi[j];
                if (!star_cap_exceeds(large, cap_lo, den)) {
                    int a = r.lo[j], b = r.hi[j];
                    // Test just below the closed lower cap bound to retain equality.
                    if (!star_cap_exceeds(r.lo, cap_lo, den))
                        return false;
                    while (b - a > 1) {
                        int mid = (a + b) / 2;
                        large[j] = mid;
                        if (star_cap_exceeds(large, cap_lo, den))
                            a = mid;
                        else
                            b = mid;
                    }
                    r.hi[j] = b;
                }
            }
        }
        int dim = kind == Kind::graph ? 4 : 3;
        r.width = 0;
        for (int t = 0; t < dim; ++t) {
            if (r.lo[t] > r.hi[t])
                return false;
            if (r.hi[t] - r.lo[t] > r.width) {
                r.width = r.hi[t] - r.lo[t];
                r.axis = t;
            }
        }
        for (int b = 0; b < 2; ++b) {
            Plane p;
            int s = b ? 1 : -1;
            if (kind == Kind::edge) {
                int a = b ? r.lo[0] : r.hi[0], z = b ? r.lo[1] : r.hi[1], w = b ? r.hi[2] : r.lo[2];
                p.n[k] = 2 * s * (z - a);
                p.n[i] = (bit(corner, i) ? 2 : -2) * s * parameter_scale;
                p.n[j] = (bit(corner, j) ? 2 : -2) * s * w;
                p.h = s * (parameter_scale + w - a - z);
            } else if (kind == Kind::graph) {
                // Affine minorant of the four signed corner heights. The smaller
                // diagonal sum fixes the intercept at the center of the face.
                array<Int, 4> a;
                for (int t = 0; t < 4; ++t)
                    a[t] = 2 * s * (b ? r.lo[t] : r.hi[t]);
                Int diff = (a[0] + a[3] - a[1] - a[2]) / 2;
                if (diff > 0) {
                    a[0] -= diff;
                    a[3] -= diff;
                } else {
                    a[1] += diff;
                    a[2] += diff;
                }
                p.n[k] = (bit(corner, k) ? 4 : -4) * s * parameter_scale;
                p.n[i] = 2 * (a[1] - a[0]);
                p.n[j] = 2 * (a[2] - a[0]);
                p.h = 2 * s * parameter_scale - a[0] - a[3];
            } else {
                const auto &a = b ? r.hi : r.lo;
                for (int t = 0; t < 3; ++t)
                    p.n[t] = (bit(corner, t) ? 2 : -2) * s * a[t];
                p.h = s * (a[0] + a[1] + a[2] - 2 * parameter_scale);
            }
            r.side[b] = p;
            // 56 boundary points with each coordinate in {-3,-1,1,3}/6.
            std::uint64_t witness = 0;
            for (int t = 0; t < 64; ++t) {
                int x = t % 4, y = t / 4 % 4, z = t / 16;
                if (x != 0 && x != 3 && y != 0 && y != 3 && z != 0 && z != 3)
                    continue;
                if (p.n[0] * (2 * x - 3) + p.n[1] * (2 * y - 3) + p.n[2] * (2 * z - 3) > 6 * p.h)
                    witness |= std::uint64_t(1) << t;
            }
            r.boundary[b] = witness;
        }
        if (method == Method::tetrahedron) {
            array<Interval, 3> n;
            if (kind == Kind::edge) {
                n[k] = {r.lo[1] - r.hi[0], r.hi[1] - r.lo[0]};
                n[i] = {parameter_scale, parameter_scale};
                n[j] = {r.lo[2], r.hi[2]};
                if (!bit(corner, i))
                    n[i] = negate(n[i]);
                if (!bit(corner, j))
                    n[j] = negate(n[j]);
            } else if (kind == Kind::graph) {
                n[k] = {parameter_scale, parameter_scale};
                if (!bit(corner, k))
                    n[k] = negate(n[k]);
                n[i] = {std::max(r.lo[1] - r.hi[0], r.lo[3] - r.hi[2]),
                        std::min(r.hi[1] - r.lo[0], r.hi[3] - r.lo[2])};
                n[j] = {std::max(r.lo[2] - r.hi[0], r.lo[3] - r.hi[1]),
                        std::min(r.hi[2] - r.lo[0], r.hi[3] - r.lo[1])};
            } else {
                for (int t = 0; t < 3; ++t) {
                    n[t] = {r.lo[t], r.hi[t]};
                    if (!bit(corner, t))
                        n[t] = negate(n[t]);
                }
            }
            Interval twice_height;
            if (kind == Kind::edge)
                twice_height = {std::max<Int>(0, parameter_scale + r.lo[2] - r.hi[0] - r.hi[1]),
                                parameter_scale + r.hi[2] - r.lo[0] - r.lo[1]};
            else if (kind == Kind::graph) {
                Int L = std::max(r.lo[0] + r.lo[3], r.lo[1] + r.lo[2]),
                    U = std::min(r.hi[0] + r.hi[3], r.hi[1] + r.hi[2]);
                if (count == 3) {
                    L = std::max(L, 17 * parameter_scale / 20);
                    U = std::min(U, parameter_scale);
                } else {
                    L = std::max(L, 2 * half_min * parameter_scale / cap_denominator);
                    U = std::min(U,
                                 (2 * (cap_denominator - half_min) * parameter_scale + cap_denominator - 1) /
                                     cap_denominator);
                }
                twice_height = {parameter_scale - U, parameter_scale - L};
            } else
                twice_height = {r.lo[0] + r.lo[1] + r.lo[2] - 2 * parameter_scale,
                                r.hi[0] + r.hi[1] + r.hi[2] - 2 * parameter_scale};
            for (int axis = 0; axis < 3; ++axis) {
                auto d = n[axis];
                auto &out = r.projection[axis];
                out.sign = d.lo > 0 ? 1 : (d.hi < 0 ? -1 : 0);
                if (!out.sign)
                    continue;
                if (out.sign < 0)
                    d = negate(d);
                if (d.lo < parameter_scale / 8) {
                    out.sign = 0;
                    continue;
                }
                for (int t = 0; t < 3; ++t) {
                    Interval value = t == 2 ? twice_height : n[(axis + t + 1) % 3];
                    if (out.sign < 0)
                        value = negate(value);
                    auto den = d;
                    if (t == 2) {
                        den.lo *= 2;
                        den.hi *= 2;
                    }
                    out.coefficient[t] = {
                        floor_division(value.lo * graph_scale, value.lo < 0 ? den.lo : den.hi),
                        -floor_division(-value.hi * graph_scale, value.hi < 0 ? den.hi : den.lo)};
                }
            }
            array<Int, 3> small{}, large{};
            for (int t = 0; t < 3; ++t) {
                if (n[t].hi < n[t].lo)
                    return false;
                small[t] = std::max({Int{0}, n[t].lo, -n[t].hi});
                large[t] = std::max(-n[t].lo, n[t].hi);
            }
            for (int t = 0; t < 3; ++t) {
                int u = (t + 1) % 3, v = (t + 2) % 3;
                const Int low = small[t] * small[t] + large[u] * large[u] + large[v] * large[v];
                const Int high = large[t] * large[t] + small[u] * small[u] + small[v] * small[v];
                const Int a = norm_floor(small[t], low),
                          b = high ? std::min(unit_grid, norm_floor(large[t], high) + 1) : unit_grid;
                r.unit[t] = {n[t].lo >= 0 ? a : -b, n[t].hi <= 0 ? -a : b};
            }
        }
        return true;
    }
    int children(int id) {
        if (rows[id].child >= 0)
            return rows[id].child;
        Row a = rows[id], b = a;
        int j = a.axis, mid = (a.lo[j] + a.hi[j]) / 2;
        a.hi[j] = b.lo[j] = mid;
        int child = int(rows.size());
        rows[id].child = child;
        a.child = b.child = -1;
        if (!prepare(a))
            a.width = -1;
        if (!prepare(b))
            b.width = -1;
        rows.push_back(a);
        rows.push_back(b);
        return child;
    }
};

#include "affine_diagonal_geometry.hpp"

using Ids = array<int, 4>;
using Sectors = array<Int, 16>;
struct SubsetKeyHash {
    std::size_t operator()(const Ids &a) const {
        std::size_t h = 0;
        for (int x : a)
            h = h * 1000003 + unsigned(x);
        return h;
    }
};
#include "affine_cell_correlations.hpp"
template <class RowDomain> struct Solver {
    CorrChecks corr;
    bool use_correlations;
    array<RowDomain, 4> dom;
    array<std::unordered_map<Ids, Sectors, SubsetKeyHash>, 16> cache;
    array<std::unordered_map<Ids, Sectors, SubsetKeyHash>, 16> upper_cache;
    Method method;
    bool use_faces = false, use_axis = false;
    Int budget;
    Int nodes = 0, leaves = 0;
    int deepest = 0;
    explicit Solver(Ids codes, Method mode = Method::two_edges, Int limit = 0)
        : use_correlations(mode == Method::correlations),
          dom{RowDomain(codes[0], mode), RowDomain(codes[1], mode), RowDomain(codes[2], mode),
              RowDomain(codes[3], mode)},
          method(use_correlations ? Method::tetrahedron : mode), budget(limit) {}
    template <std::size_t N, Bound bound = Bound::Lower> Int lower(const Ids &ids, int subset, int m) const {
        if constexpr (std::is_same_v<RowDomain, DiagonalDomain>)
            return diagonal_volume<N, bound>(dom, ids, subset, m);
        array<Vector, N> ns;
        array<Int, N> hs;
        int j = 0;
        for (int i = 0; i < 4; ++i)
            if (bit(subset, i)) {
                const auto &p = dom[i].rows[ids[i]].side[bit(m, i) ^ (bound == Bound::Upper)];
                ns[j] = p.n;
                if constexpr (bound == Bound::Upper) {
                    for (auto &x : ns[j])
                        x = -x;
                    hs[j++] = -p.h;
                } else
                    hs[j++] = p.h;
            }
        if (use_faces)
            return cube_face_integral::volume<bound>(ns, hs);
        return volume<bound>(ns, hs, cube);
    }
    const Sectors &sectors(const Ids &ids, int subset) {
        Ids key = ids;
        for (int i = 0; i < 4; ++i)
            if (!bit(subset, i))
                key[i] = 0;
        auto &tab = cache[subset];
        auto it = tab.find(key);
        if (it != tab.end())
            return it->second;
        Sectors v{};
        int n = population(subset);
        for (int m = subset;; m = (m - 1) & subset) {
            if (n == 1)
                v[m] = lower<1>(ids, subset, m);
            else if (n == 2)
                v[m] = lower<2>(ids, subset, m);
            else
                v[m] = lower<3>(ids, subset, m);
            if (m == 0)
                break;
        }
        // Clearing a full cache only causes exact recomputation.
        if (tab.size() >= 100000)
            tab.clear();
        return tab.emplace(key, v).first->second;
    }
    const Sectors &upper_triples(const Ids &ids, int omitted) {
        int sub = 15 ^ (1 << omitted);
        Ids key = ids;
        key[omitted] = 0;
        auto &tab = upper_cache[sub];
        auto it = tab.find(key);
        if (it != tab.end())
            return it->second;
        Sectors v{};
        for (int m = sub;; m = (m - 1) & sub) {
            v[m] = lower<3, Bound::Upper>(ids, sub, m);
            if (!m)
                break;
        }
        if (tab.size() >= 100000)
            tab.clear();
        return tab.emplace(key, v).first->second;
    }
    bool marginal_feasible(int g, const array<const Sectors *, 4> &L, const array<const Sectors *, 4> &U) {
        constexpr Int ceiling = (cap_numerator * volume_grid + cap_denominator - 1) / cap_denominator;
        Sectors lo{}, hi;
        hi.fill(ceiling);
        hi[g] = 0;
        if (method == Method::tetrahedron)
            hi[15 ^ g] = (tetra_cap * volume_grid + cap_denominator - 1) / cap_denominator;
        for (int sweep = 0; sweep < 3; ++sweep) {
            for (int i = 0; i < 4; ++i) {
                int d = 1 << i, sub = 15 ^ d;
                for (int a = sub;; a = (a - 1) & sub) {
                    int b = a | d;
                    const Int l = (*L[i])[a], u = (*U[i])[a];
                    lo[a] = std::max(lo[a], l - hi[b]);
                    lo[b] = std::max(lo[b], l - hi[a]);
                    hi[a] = std::min(hi[a], u - lo[b]);
                    hi[b] = std::min(hi[b], u - lo[a]);
                    if (lo[a] > hi[a] || lo[b] > hi[b])
                        return false;
                    if (!a)
                        break;
                }
            }
            Int l = 0, u = 0;
            for (int i = 0; i < 16; ++i) {
                l += lo[i];
                u += hi[i];
            }
            if (l > volume_grid || u < volume_grid)
                return false;
        }
        return true;
    }
    // Uniform upper bound for the difference of two clipped affine graphs.
    static Int graph_distance(const Projection &p, const Projection &q) {
        array<Int, 3> d{};
        for (int k = 0; k < 3; ++k) {
            auto a = p.coefficient[k], b = q.coefficient[k];
            d[k] = std::max(a.hi - b.lo, b.hi - a.lo);
        }
        Int a = std::max(d[0], d[1]), b = std::min(d[0], d[1]), c = d[2];
        Wide den = graph_scale, num = c;
        if (2 * c < a + b) {
            if (2 * c + b <= a) {
                num = 3 * a * a + 12 * c * c + b * b;
                den = Wide(12) * a * graph_scale;
            } else {
                Int tail = a + b - 2 * c;
                num = Wide(24) * a * b * c + Wide(tail) * tail * tail;
                den = Wide(24) * a * b * graph_scale;
            }
        }
        return -floor_division(-num * volume_grid, den);
    }
    bool cover(Ids ids, unsigned possible, int depth = 0) {
        ++nodes;
        if (budget && nodes > budget)
            return false;
        deepest = std::max(deepest, depth);
        int best = -1, width = 0;
        for (int i = 0; i < 4; ++i) {
            int w = dom[i].rows[ids[i]].width;
            if (w < 0) {
                ++leaves;
                return true;
            }
            if (w > width) {
                width = w;
                best = i;
            }
        }
        if (method == Method::tetrahedron) {
            for (int i = 0; i < 4 && possible; ++i)
                for (int j = 0; j < i && possible; ++j) {
                    Interval dot{};
                    for (int k = 0; k < 3; ++k) {
                        auto v = product(dom[i].rows[ids[i]].unit[k], dom[j].rows[ids[j]].unit[k]);
                        dot.lo += v.lo;
                        dot.hi += v.hi;
                    }
                    if (use_axis) {
                        for (int k = 0; k < 3; ++k) {
                            const auto &p = dom[i].rows[ids[i]].projection[k],
                                       &q = dom[j].rows[ids[j]].projection[k];
                            if (!p.sign || !q.sign)
                                continue;
                            Int upper = graph_distance(p, q) * cap_denominator;
                            for (int opposite = 0; opposite < 2; ++opposite) {
                                bool same = (p.sign == q.sign) ^ bool(opposite);
                                Int required = same ? cap_denominator - 6 * cap_numerator - 8452
                                                    : cap_denominator - 8 * cap_numerator;
                                if (upper < required * volume_grid)
                                    for (int g = 0; g < 16; ++g)
                                        if ((bit(g, i) ^ bit(g, j)) == opposite)
                                            possible &= ~(1u << g);
                            }
                        }
                    }
                    for (int opposite = 0; opposite < 2; ++opposite) {
                        auto v = opposite ? negate(dot) : dot;
                        if (v.lo * cap_denominator > (use_axis ? 143814 : 329912) * unit_grid * unit_grid ||
                            v.hi * cap_denominator < -(use_axis ? 573370 : 574310) * unit_grid * unit_grid)
                            for (int g = 0; g < 16; ++g)
                                if ((bit(g, i) ^ bit(g, j)) == opposite)
                                    possible &= ~(1u << g);
                    }
                }
        }
        // The empty exterior cell and its bounded opposite have no strict
        // cube-boundary point.
        for (int m = 0; m < 16; ++m)
            if (possible & (1u << m)) {
                for (int d : {m, 15 ^ m}) {
                    std::uint64_t mask = ~std::uint64_t(0);
                    for (int i = 0; i < 4; ++i)
                        mask &= dom[i].rows[ids[i]].boundary[bit(d, i)];
                    if (mask) {
                        possible &= ~(1u << m);
                        break;
                    }
                }
            }
        for (int sub : {1, 2, 4, 8, 3, 5, 6, 9, 10, 12, 7, 11, 13, 14}) {
            if (!possible)
                break;
            const auto &lower_by_sign = sectors(ids, sub);
            // A partial sign cell has this many completions; one is absent if it
            // matches the exterior word.
            int cap = 1 << (4 - population(sub));
            for (int m = sub;; m = (m - 1) & sub) {
                if (lower_by_sign[m] * cap_denominator > cap * cap_numerator * volume_grid) {
                    possible = 0;
                    break;
                }
                if (lower_by_sign[m] * cap_denominator > (cap - 1) * cap_numerator * volume_grid) {
                    for (int t = 15 ^ sub;; t = (t - 1) & (15 ^ sub)) {
                        possible &= ~(1u << (m | t));
                        if (method == Method::tetrahedron &&
                            lower_by_sign[m] * cap_denominator >
                                ((cap - 1) * cap_numerator + tetra_cap) * volume_grid)
                            possible &= ~(1u << (15 ^ (m | t)));
                        if (!t)
                            break;
                    }
                }
                if (!m)
                    break;
            }
        }
        if constexpr (std::is_same_v<RowDomain, DiagonalDomain>) {
            if (use_correlations && possible)
                possible = corr.feasible(*this, ids, possible);
        }
        if (method >= Method::marginals && possible) {
            array<const Sectors *, 4> L, U;
            for (int i = 0; i < 4; ++i) {
                L[i] = &sectors(ids, 15 ^ (1 << i));
                U[i] = &upper_triples(ids, i);
            }
            for (int g = 0; g < 16; ++g)
                if ((possible & (1u << g)) && !marginal_feasible(g, L, U))
                    possible &= ~(1u << g);
        }
        if (possible) {
            for (int m = 0; m < 16; ++m)
                if (lower<4>(ids, 15, m) * cap_denominator > cap_numerator * volume_grid) {
                    possible = 0;
                    break;
                }
        }
        if (!possible) {
            ++leaves;
            return true;
        }
        if (width <= 1 || depth == 90)
            return false;
        int a = dom[best].children(ids[best]);
        ids[best] = a;
        if (!cover(ids, possible, depth + 1))
            return false;
        ids[best] = a + 1;
        return cover(ids, possible, depth + 1);
    }
    bool run() {
        unsigned p = 65535;
        for (int v = 0; v < 8; ++v) {
            int m = 0;
            for (int i = 0; i < 4; ++i)
                m |= bit(dom[i].bits, v) << i;
            p &= ~((1u << m) | (1u << (15 ^ m)));
        }
        if (use_axis) {
            for (int g = 0; g < 16; ++g) {
                unsigned axes = 0;
                for (int i = 0; i < 4; ++i)
                    if (dom[i].kind == Domain::Kind::graph) {
                        int k = dom[i].k, side = bit(dom[i].corner, k) ^ bit(g, i);
                        unsigned flag = 1u << (2 * k + side);
                        if (axes & flag)
                            p &= ~(1u << g);
                        axes |= flag;
                    }
            }
        }
        bool ok = cover({}, p);
        if (ok)
            assert(2 * leaves == nodes + 1);
        return ok;
    }
};
} // namespace
int main(int argc, char **argv) {
    bool use_faces = false, use_axis = false, diagonal = false;
    Method requested = Method::tetrahedron;
    int arg = 1;
    if (arg < argc && std::string_view(argv[arg]).starts_with("--")) {
        std::string_view option = argv[arg++];
        if (option == "--graph")
            requested = Method::graph;
        else if (option == "--marginals")
            requested = Method::marginals;
        else if (option == "--correlations") {
            diagonal = true;
            requested = Method::correlations;
        } else if (option == "--diagonal")
            diagonal = true;
        else if (option == "--axis") {
            use_faces = true;
            use_axis = true;
        } else if (option == "--faces")
            use_faces = true;
        else if (option != "--tetrahedron")
            return 2;
    }
    if (argc > arg + 1)
        return 2;
    Int budget = arg < argc ? std::atoll(argv[arg]) : 0;
    Ids c;
    int fields, cases = 0;
    bool complete = true;
    while ((fields = std::scanf("%d %d %d %d", &c[0], &c[1], &c[2], &c[3])) == 4) {
        if (std::any_of(c.begin(), c.end(), [](int x) { return x < 1 || x > 127; }))
            return 2;
        ++cases;
        bool extended = std::count_if(c.begin(), c.end(), [](int x) {
                            int n = population(x);
                            return n == 2 || n == 6;
                        }) < 2;
        Int nodes = 0, leaves = 0;
        int depth = 0;
        bool ok = true;
        if (diagonal) {
            if (!extended)
                return 2;
            auto charts = diagonal_charts(c);
            // Represent the finite union of direction charts as a balanced binary
            // tree.
            auto run_charts = [&](auto &&visit, int lo, int hi, int level) -> bool {
                if (hi - lo > 1) {
                    ++nodes;
                    int mid = (lo + hi) / 2;
                    return visit(visit, lo, mid, level + 1) && visit(visit, mid, hi, level + 1);
                }
                Solver<DiagonalDomain> s(charts[lo], requested,
                                         budget ? std::max<Int>(1, budget - nodes) : 0);
                s.use_axis = true;
                bool closed = s.run();
                nodes += s.nodes;
                leaves += s.leaves;
                depth = std::max(depth, level + s.deepest);
                return closed;
            };
            ok = run_charts(run_charts, 0, int(charts.size()), 0);
        } else {
            Solver<Domain> s(c, extended ? requested : Method::two_edges, budget);
            s.use_faces = use_faces;
            s.use_axis = use_axis;
            ok = s.run();
            nodes = s.nodes;
            leaves = s.leaves;
            depth = s.deepest;
        }
        if (ok)
            assert(2 * leaves == nodes + 1);
        complete = complete && ok;
        std::printf("%d %d %d %d  %s %lld %lld %d\n", c[0], c[1], c[2], c[3], ok ? "proved" : "OPEN",
                    (long long)nodes, (long long)leaves, depth);
        std::fflush(stdout);
    }
    if (fields != EOF || !cases)
        return 2;
    return complete ? 0 : 1;
}
