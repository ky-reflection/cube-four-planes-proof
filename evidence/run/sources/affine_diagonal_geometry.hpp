// Included after the common row data types inside the cover program.
// Diagonal coordinates on a face and the 24 center-to-edge tetrahedra of the
// cube.
constexpr int graph_grid = 2048;
constexpr int mesh_tets = 24;
using Mesh = std::array<std::array<int, 2>, 15>;
static const auto spatial_points = [] {
    std::array<std::array<int, 3>, 15> list{};
    for (int a = 0; a < 3; ++a)
        for (int s = 0; s < 2; ++s)
            list[1 + 2 * a + s][a] = 2 * s - 1;
    for (int v = 0; v < 8; ++v)
        for (int a = 0; a < 3; ++a)
            list[7 + v][a] = 2 * ((v >> a) & 1) - 1;
    return list;
}();
static const auto tetrahedra = [] {
    std::array<std::array<int, 4>, mesh_tets> list{};
    int t = 0;
    for (int a = 0; a < 3; ++a)
        for (int s = 0; s < 2; ++s)
            for (int b = 0; b < 3; ++b)
                if (a != b)
                    for (int r = 0; r < 2; ++r) {
                        int v = (s << a) | (r << b);
                        list[t++] = {0, 1 + 2 * a + s, 7 + v, 7 + v + (1 << (3 - a - b))};
                    }
    return list;
}();
constexpr Int slice_grid = Int{1} << 29;
template <bool upper> Int round_slice(Int num, Int den) {
    if constexpr (upper)
        return num / den + (num % den > 0);
    else
        return num / den - (num % den < 0);
}
using SliceScores = array<array<int, 4>, 4>;
template <bool upper> Int slice(SliceScores &a, int n, Int w) {
    int chosen = -1, pos = -1, neg = -1;
    for (int j = 0; j < n; ++j) {
        int p = 0, m = 0, ip = -1, in = -1;
        for (int v = 0; v < 4; ++v) {
            if (a[j][v] > 0) {
                ++p;
                ip = v;
            } else if (a[j][v] < 0) {
                ++m;
                in = v;
            }
        }
        if (m && !p)
            return 0;
        if (m && chosen < 0) {
            chosen = j;
            pos = ip;
            neg = in;
        }
    }
    if (chosen < 0)
        return w;
    if (w == 0)
        return 0;
    Int A = a[chosen][pos], B = -a[chosen][neg], d = A + B;
    int cp[4], cn[4], val[4];
    for (int j = 0; j < n; ++j) {
        cp[j] = a[j][pos];
        cn[j] = a[j][neg];
        val[j] = int(round_slice<upper>(B * cp[j] + A * cn[j], d));
    }
    for (int j = 0; j < n; ++j)
        a[j][neg] = val[j];
    Int left = slice<upper>(a, n, round_slice<upper>(w * A, d));
    for (int j = 0; j < n; ++j) {
        a[j][neg] = cn[j];
        a[j][pos] = val[j];
    }
    Int right = slice<upper>(a, n, round_slice<upper>(w * B, d));
    for (int j = 0; j < n; ++j)
        a[j][pos] = cp[j];
    return left + right;
}
struct DiagonalRow : Row {
    Mesh centered{};
    Interval cap{};
    Mesh mesh{};
};
Int cap_at(int e, std::array<Int, 3> weights, bool upper) {
    Int z = 2 * e, den = 1;
    int d = 0, active = 0;
    for (int i = 0; i < 3; ++i)
        if (weights[i]) {
            z += weights[i];
            den *= 2 * weights[i] * (++d);
            active |= 1 << i;
        }
    Int sum = 0;
    for (int m = active;; m = (m - 1) & active) {
        Int v = z;
        for (int i = 0; i < 3; ++i)
            if (m >> i & 1)
                v -= 2 * weights[i];
        if (v > 0) {
            Int t = 1;
            for (int k = 0; k < d; ++k)
                t *= v;
            sum += (population(m) & 1) ? -t : t;
        }
        if (!m)
            break;
    }
    return upper ? -floor_division(-Wide(sum) * volume_grid, den)
                 : floor_division(Wide(sum) * volume_grid, den);
}

struct DiagonalDomain {
    using Kind = Domain::Kind;
    int k, i, j, sign, bits, corner, count;
    Kind kind;
    bool track_caps;
    std::vector<DiagonalRow> rows;
    explicit DiagonalDomain(int encoded, Method mode) : track_caps(mode == Method::correlations) {
        k = encoded / 256;
        int split = encoded % 256;
        bits = population(split) > 4 ? 255 ^ split : split;
        count = population(bits);
        corner = first(bits);
        bool face = false;
        for (int d = 0; d < 3; ++d) {
            unsigned sides = 0;
            for (int v = 0; v < 8; ++v)
                if (bit(bits, v))
                    sides |= 1u << bit(v, d);
            face |= sides != 3;
        }
        kind = count == 2 ? Kind::edge : (face ? Kind::graph : Kind::star);
        i = (k + 1) % 3;
        j = (k + 2) % 3;
        sign = 0;
        for (int v = 0; v < 8; ++v)
            if (!bit(v, k) && bit(bits, v) != bit(bits, v ^ (1 << k))) {
                int wanted = bit(bits, v) ? -1 : 1;
                assert(!sign || sign == wanted);
                sign = wanted;
            }
        assert(sign);
        constexpr int G = graph_grid;
        DiagonalRow r;
        r.lo = {-int((232 * G + 999) / 1000), -2 * G, -2 * G, 0};
        r.hi = {-r.lo[0], 2 * G, 2 * G, 0};
        if (count < 4)
            r.hi[0] = 0;
        if (kind == Kind::graph) {
            int cap = (count == 4 ? (133530 * G + 999999) / 1000000 : (15 * G + 99) / 100);
            r.lo[0] = -cap;
            r.hi[0] = count == 4 ? cap : 0;
        }
        if (!prepare(r))
            r.width = -1;
        rows.push_back(r);
    }
    bool prepare(DiagonalRow &r) const {
        constexpr int G = graph_grid;
        for (int sweep = 0; sweep < 3; ++sweep) {
            for (int v = 0; v < 8; ++v) {
                int x = 2 * bit(v, i) - 1, y = 2 * bit(v, j) - 1, z = 2 * bit(v, k) - 1;
                int d = x == y ? 1 : 2, t = y, s = bit(bits, v) ? 1 : -1;
                // s(e+t*A_d+sign*z) >= 0.
                int C = -s * sign * z * G;
                auto clip = [&](int a, int da, int b, int db) {
                    int lo = C - (db > 0 ? r.hi[b] : -r.lo[b]);
                    if (da > 0)
                        r.lo[a] = std::max(r.lo[a], lo);
                    else
                        r.hi[a] = std::min(r.hi[a], -lo);
                };
                clip(0, s, d, s * t);
                clip(d, s * t, 0, s);
            }
            for (int d : {1, 2}) {
                int v = 3 - d, u = 2 * G - std::max({0, r.lo[v], -r.hi[v]});
                r.lo[d] = std::max(r.lo[d], -u);
                r.hi[d] = std::min(r.hi[d], u);
            }
        }
        r.width = 0;
        for (int d = 0; d < 3; ++d) {
            int w = r.hi[d] - r.lo[d];
            if (w < 0)
                return false;
            if (w <= 1)
                continue;
            w *= d == 0 ? 3 : 1;
            if (w > r.width) {
                r.width = w;
                r.axis = d;
            }
        }
        auto interval = [&](int a, int b, int c) {
            int coeff[] = {a, b, c};
            Interval h;
            for (int v = 0; v < 3; ++v) {
                int lo = coeff[v] * r.lo[v], hi = coeff[v] * r.hi[v];
                h.lo += std::min(lo, hi);
                h.hi += std::max(lo, hi);
            }
            return h;
        };
        array<Interval, 15> h;
        for (int p = 0; p < 15; ++p) {
            auto xyz = spatial_points[p];
            auto f = interval(2, xyz[i] + xyz[j], xyz[j] - xyz[i]);
            f.lo += 2 * sign * G * xyz[k];
            f.hi += 2 * sign * G * xyz[k];
            if (p >= 7) {
                if (bit(bits, p - 7))
                    f.lo = std::max<Int>(0, f.lo);
                else
                    f.hi = std::min<Int>(0, f.hi);
            }
            h[p] = f;
        }
        // At a face center, each of its opposite corner pairs averages to that
        // center.
        for (int a = 0; a < 3; ++a)
            for (int s = 0; s < 2; ++s) {
                auto &f = h[1 + 2 * a + s];
                for (int v = 0; v < 8; ++v)
                    if (bit(v, a) == s) {
                        auto l = h[7 + v], u = h[7 + (v ^ (7 ^ (1 << a)))];
                        f.lo = std::max(f.lo, floor_division(l.lo + u.lo, 2));
                        f.hi = std::min(f.hi, -floor_division(-l.hi - u.hi, 2));
                    }
            }
        for (int p = 0; p < 15; ++p) {
            if (h[p].lo > h[p].hi)
                return false;
            r.mesh[p] = {int(-h[p].hi * 64), int(h[p].lo * 64)};
        }
        array<Interval, 3> n;
        n[k] = {2 * sign * G, 2 * sign * G};
        n[i] = interval(0, 1, -1);
        n[j] = interval(0, 1, 1);
        for (auto &t : n) {
            t.lo = std::max<Int>(-2 * G, t.lo);
            t.hi = std::min<Int>(2 * G, t.hi);
            if (t.lo > t.hi)
                return false;
        }
        if (track_caps) {
            std::array<Int, 3> l{}, u{};
            for (int t = 0; t < 3; ++t) {
                l[t] = std::max({Int(0), n[t].lo, -n[t].hi});
                u[t] = std::max(n[t].hi, -n[t].lo);
            }
            r.cap = {cap_at(r.lo[0], r.lo[0] < 0 ? l : u, false),
                     cap_at(r.hi[0], r.hi[0] < 0 ? u : l, true)};
            if (r.cap.lo * cap_denominator > (cap_denominator - 433235) * volume_grid ||
                r.cap.hi * cap_denominator < 433235 * volume_grid)
                return false;
            for (int p = 0; p < 15; ++p) {
                auto xyz = spatial_points[p];
                auto f = interval(0, xyz[i] + xyz[j], xyz[j] - xyz[i]);
                f.lo = std::max(f.lo + 2 * sign * G * xyz[k], h[p].lo - 2 * r.hi[0]);
                f.hi = std::min(f.hi + 2 * sign * G * xyz[k], h[p].hi - 2 * r.lo[0]);
                if (f.lo > f.hi)
                    return false;
                r.centered[p] = {int(-f.hi * 64), int(f.lo * 64)};
            }
        }
        for (int t = 0; t < 3; ++t) {
            auto a = n[t];
            int u = (t + 1) % 3, v = (t + 2) % 3;
            auto minabs = [](Interval x) { return std::max({Int{0}, x.lo, -x.hi}); };
            auto maxabs = [](Interval x) { return std::max(-x.lo, x.hi); };
            Int m = minabs(a), M = maxabs(a),
                l = m * m + maxabs(n[u]) * maxabs(n[u]) + maxabs(n[v]) * maxabs(n[v]);
            Int h = M * M + minabs(n[u]) * minabs(n[u]) + minabs(n[v]) * minabs(n[v]);
            Int L = norm_floor(m, l), H = h ? std::min(unit_grid, norm_floor(M, h) + 1) : unit_grid;
            r.unit[t] = {a.lo >= 0 ? L : -H, a.hi <= 0 ? -L : H};
        }
        // The proof only needs vertex witnesses: any strict sign word there
        // excludes both exterior choices.
        return true;
    }
    int children(int id) {
        if (rows[id].child >= 0)
            return rows[id].child;
        DiagonalRow a = rows[id], b = a;
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

template <std::size_t N, Bound bound, bool centered = false>
Int diagonal_volume(const std::array<DiagonalDomain, 4> &dom, const std::array<int, 4> &ids, int subset,
                    int m) {
    constexpr bool upper = bound == Bound::Upper;
    array<const Mesh *, N> mesh{};
    array<int, N> side{};
    int k = 0;
    for (int i = 0; i < 4; ++i)
        if (bit(subset, i)) {
            mesh[k] = centered ? &dom[i].rows[ids[i]].centered : &dom[i].rows[ids[i]].mesh;
            side[k++] = bit(m, i) ^ upper;
        }
    Int sum = 0;
    for (auto &tet : tetrahedra) {
        SliceScores a{};
        int active = 0;
        bool empty = false;
        for (size_t i = 0; i < N; ++i) {
            int p = 0, neg = 0;
            for (int v = 0; v < 4; ++v) {
                int h = (*mesh[i])[tet[v]][side[i]] * (upper ? -1 : 1);
                a[active][v] = h;
                p += h > 0;
                neg += h < 0;
            }
            if (neg && !p) {
                empty = true;
                break;
            }
            if (neg)
                ++active;
        }
        if (!empty)
            sum += active ? slice<upper>(a, active, slice_grid) : slice_grid;
    }
    constexpr Int den = mesh_tets * (slice_grid / volume_grid);
    return round_slice<upper>(sum, den);
}
std::vector<std::array<int, 4>> diagonal_charts(const std::array<int, 4> &codes) {
    std::array<std::vector<int>, 4> axes;
    for (int a = 0; a < 4; ++a) {
        int bits = population(codes[a]) > 4 ? 255 ^ codes[a] : codes[a], v = first(bits), same = 7;
        for (int u = 0; u < 8; ++u)
            if (bit(bits, u))
                same &= ~(u ^ v);
        for (int k = 0; k < 3; ++k)
            if (!same || bit(same, k))
                axes[a].push_back(k);
    }
    std::vector<std::array<int, 4>> charts;
    for (int i : axes[0])
        for (int j : axes[1])
            for (int k : axes[2])
                for (int l : axes[3]) {
                    charts.push_back(
                        {codes[0] + 256 * i, codes[1] + 256 * j, codes[2] + 256 * k, codes[3] + 256 * l});
                }
    return charts;
}
