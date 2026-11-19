#ifndef CUBE_AFFINE_FACE_INTEGRAL_HPP
#define CUBE_AFFINE_FACE_INTEGRAL_HPP
// Integrate the signed pyramids from the origin to each supporting facet.
#include "seven_integer_volume.hpp"
namespace cube_face_integral {
using namespace seven_integer;
struct Constraint {
    Vector n;
    Int h;
};
constexpr Int coord_grid = Int{1} << 23;
struct L {
    Int a, b, c;
}; // a*x+b*y >= c
// Each polygon vertex is at the intersection of lines left and right.
// These two indices also identify the two incident polygon edges.
struct V {
    Wide x, y, d;
    int left, right;
};
inline V vertex(const L *line, int p, int q) {
    auto a = line[p], b = line[q];
    Wide d = Wide(a.a) * b.b - Wide(a.b) * b.a, x = Wide(a.c) * b.b - Wide(a.b) * b.c,
         y = Wide(a.a) * b.c - Wide(a.c) * b.a;
    assert(d);
    if (d < 0) {
        d = -d;
        x = -x;
        y = -y;
    }
    return {x, y, d, p, q};
}
inline std::pair<Int, Int> projected_area(const L *input, int count) {
    L l[16] = {{2, 0, -1}, {0, 2, -1}, {-2, 0, -1}, {0, -2, -1}};
    V a[24] = {{-1, -1, 2, 0, 1}, {1, -1, 2, 1, 2}, {1, 1, 2, 2, 3}, {-1, 1, 2, 3, 0}}, b[24];
    int n = 4;
    for (int h = 4; h < 4 + count; ++h) {
        l[h] = input[h - 4];
        int m = 0;
        int s[24];
        bool neg = false, pos = false;
        for (int i = 0; i < n; ++i) {
            auto v = a[i];
            Wide f = Wide(l[h].a) * v.x + Wide(l[h].b) * v.y - Wide(l[h].c) * v.d;
            s[i] = (f > 0) - (f < 0);
            neg |= s[i] < 0;
            pos |= s[i] > 0;
        }
        if (!neg)
            continue;
        if (!pos)
            return {0, 0};
        for (int i = 0; i < n; ++i) {
            int j = (i + 1) % n, z = (i + n - 1) % n;
            auto v = a[i];
            if (s[i] >= 0) {
                if (!s[i]) {
                    if (s[z] < 0)
                        v.left = h;
                    if (s[j] < 0)
                        v.right = h;
                }
                b[m++] = v;
            }
            if (s[i] * s[j] < 0) {
                int same = a[i].right;
                assert(same == a[j].left);
                b[m++] = s[i] > 0 ? vertex(l, same, h) : vertex(l, h, same);
            }
        }
        n = m;
        std::copy_n(b, n, a);
    }
    if (n < 3)
        return {0, 0};
    Point p[24];
    for (int i = 0; i < n; ++i)
        p[i] = {floor_division(Wide(a[i].x) * coord_grid, a[i].d),
                floor_division(Wide(a[i].y) * coord_grid, a[i].d)};
    Int twice = 0;
    for (int i = 0; i < n; ++i) {
        auto x = p[i], y = p[(i + 1) % n];
        twice += x.x * y.y - x.y * y.x;
    }
    // Each coordinate error is <1 on the spatial grid, and all |x| <= grid/2.
    Int error = n * (2 * coord_grid + 2);
    return {std::max<Int>(0, twice - error), std::min(2 * coord_grid * coord_grid, twice + error)};
}
inline bool zero(const Vector &v) { return !v[0] && !v[1] && !v[2]; }
inline Vector cross3(const Vector &a, const Vector &b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
template <Bound which, std::size_t N>
Int volume(const std::array<Vector, N> &normals, const std::array<Int, N> &offsets) {
    static_assert(N >= 1 && N <= 4);
    for (auto row : normals)
        for (Int v : row)
            assert(std::abs(v) <= parameter_grid);
    for (Int h : offsets)
        assert(std::abs(h) <= parameter_grid);
    Constraint cuts[6 + N];
    int count = 0;
    for (int i = 0; i < 3; ++i)
        for (int s : {-2, 2}) {
            Vector n{};
            n[i] = s;
            cuts[count++] = {n, -1};
        }
    for (size_t j = 0; j < N; ++j) {
        const auto &n = normals[j];
        Int h = offsets[j];
        if (zero(n)) {
            if (h > 0)
                return 0;
            continue;
        }
        bool used = false;
        for (int k = 0; k < count; ++k) {
            auto v = cuts[k];
            if (!zero(cross3(n, v.n)))
                continue;
            int t = 0;
            while (!n[t])
                ++t;
            Int a = std::abs(n[t]), b = std::abs(v.n[t]);
            if ((n[t] > 0) == (v.n[t] > 0)) {
                if (h * b > v.h * a)
                    cuts[k] = {n, h};
                used = true;
                break;
            }
            if (h * b >= -v.h * a)
                return 0;
        }
        if (!used)
            cuts[count++] = {n, h};
    }
    Int answer = 0;
    for (int f = 0; f < count; ++f) {
        auto p = cuts[f];
        if (!p.h)
            continue;
        int k = 0;
        for (int i = 1; i < 3; ++i)
            if (std::abs(p.n[i]) > std::abs(p.n[k]))
                k = i;
        int i = (k + 1) % 3, j = (k + 2) % 3, s = p.n[k] > 0 ? 1 : -1;
        L lines[12];
        int m = 0;
        bool empty = false;
        for (int g = 0; g < count; ++g)
            if (g != f) {
                auto q = cuts[g];
                if (!q.n[k] && (!q.n[i] || !q.n[j]) && q.h * 2 == -std::abs(q.n[i] + q.n[j]))
                    continue;
                L line{s * (q.n[i] * p.n[k] - q.n[k] * p.n[i]), s * (q.n[j] * p.n[k] - q.n[k] * p.n[j]),
                       s * (q.h * p.n[k] - q.n[k] * p.h)};
                if (line.a || line.b)
                    lines[m++] = line;
                else if (line.c > 0) {
                    empty = true;
                    break;
                }
            }
        if (empty)
            continue;
        auto area = projected_area(lines, m);
        Int value = ((-p.h > 0) == (which == Bound::Lower)) ? area.first : area.second;
        Wide num = Wide(-p.h) * value * volume_grid,
             den = Wide(6) * coord_grid * coord_grid * std::abs(p.n[k]);
        answer += which == Bound::Lower ? floor_division(num, den) : -floor_division(-num, den);
    }
    return std::clamp(answer, Int{0}, volume_grid);
}
} // namespace cube_face_integral

#endif
