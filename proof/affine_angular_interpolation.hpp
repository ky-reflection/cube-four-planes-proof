// Pair disagreements, angular derivatives and multilinear interpolation.
struct PhysicalPair {
    bool has[2][2]{};
    I at{}, deriv[2][2]{}, hessian[2][2]{}, corner[16]{};
};
I signed_div(I a, I b) { return b.h < 0 ? divi(neg(a), neg(b)) : divi(a, b); }
I dmean_at(Int a, Int b) {
    if (a < 0)
        return neg(dmean_at(-a, b));
    if (a == 0)
        return {};
    return a < b ? divi({a, a}, {6 * b, 6 * b})
                 : sub({mass / 4, mass / 4}, divide(sq(divi({b, b}, {a, a})), 12));
}
I dmean(I a, I b) {
    I q{mass, -mass};
    for (Int x : {a.l, a.h})
        for (Int y : {absmin(b), absmax(b)}) {
            auto v = dmean_at(x, y);
            q.l = std::min(q.l, v.l);
            q.h = std::max(q.h, v.h);
        }
    return q;
}
I hull(I a, I b) { return {std::min(a.l, b.l), std::max(a.h, b.h)}; }
struct DirectionPairs {
    using S = std::array<I, 8>;
    using Norm = std::array<I, 3>;
    struct CacheKey {
        std::array<int, 8> p;
        bool operator==(const CacheKey &) const = default;
    };
    struct Hash {
        size_t operator()(const CacheKey &k) const {
            size_t h = 0;
            for (int x : k.p)
                h = 131 * h + unsigned(x);
            return h;
        }
    };
    std::unordered_map<CacheKey, PhysicalPair, Hash> store[6];
    static int sg(int i, int k) { return i == 0 || i == k + 1 ? 1 : -1; }
    static I scores(Coordinate c) { return {c.lo * (mass / grid), c.hi * (mass / grid)}; }
    static int coeff(int r, int h) {
        auto [i, j] = edges[h];
        return r < 4 ? -(r != i && r != j) : int(((reps[r] >> i) ^ (reps[r] >> j)) & 1);
    }
    std::unordered_map<CacheKey, I, Hash> center_cache[6];
    I corner_value(const Matrix &c, const Pivots &p, int h, int vertex) {
        auto [i, j] = edges[h];
        CacheKey key{};
        int iz = 0;
        std::array<seven_integer::Vector, 2> n{};
        for (int w = 0; w < 2; ++w) {
            int r = w ? j : i;
            for (int col = 0; col < 4; ++col)
                if (col != r) {
                    int t = grid;
                    if (col != p[r]) {
                        t = vertex >> iz & 1 ? c[r][col].hi : c[r][col].lo;
                        key.p[iz++] = t;
                    }
                    for (int k = 0; k < 3; ++k)
                        if (col != (r ^ (k + 1)))
                            n[w][k] += sg(r, k) * t;
                }
        }
        auto &cache = center_cache[h];
        auto it = cache.find(key);
        if (it != cache.end())
            return it->second;
        Wide lo = 0, hi = 0;
        for (int k = 0; k < 3; ++k)
            for (int face : {-1, 1}) {
                cube_face_integral::L lines[2];
                for (int w = 0; w < 2; ++w) {
                    int t = w ? -1 : 1;
                    lines[w] = {2 * t * n[w][(k + 1) % 3], 2 * t * n[w][(k + 2) % 3], -t * face * n[w][k]};
                }
                auto area = cube_face_integral::projected_area(lines, 2);
                lo += area.first;
                hi += area.second;
            }
        constexpr Int gc = cube_face_integral::coord_grid;
        I q{floor_ratio(lo * mass, 6 * gc * gc), ceil_ratio(hi * mass, 6 * gc * gc)};
        if (cache.size() > 200000)
            cache.clear();
        cache.emplace(key, q);
        return q;
    }
    PhysicalPair compute(const Matrix &c, const Pivots &piv, int h) {
        auto [i, j] = edges[h];
        CacheKey key;
        int iz = 0;
        for (int r : {i, j})
            for (int col = 0; col < 4; ++col)
                if (r != col && piv[r] != col) {
                    key.p[iz++] = c[r][col].lo;
                    key.p[iz++] = c[r][col].hi;
                }
        auto &table = store[h];
        auto it = table.find(key);
        if (it != table.end())
            return it->second;
        Norm n[2]{};
        std::array<seven_integer::Vector, 2> center{};
        for (int w = 0; w < 2; ++w) {
            int r = w ? j : i;
            for (int k = 0; k < 3; ++k) {
                for (int col = 0; col < 4; ++col)
                    if (col != r && col != (r ^ (k + 1))) {
                        n[w][k] = add(n[w][k], scores(c[r][col]));
                        center[w][k] += sg(r, k) * (c[r][col].lo + c[r][col].hi);
                    }
                n[w][k] = scalar(divide(n[w][k], 2), sg(r, k));
            }
        }
        PhysicalPair out;
        Wide lo = 0, hi = 0;
        for (int k = 0; k < 3; ++k)
            for (int face : {-1, 1}) {
                cube_face_integral::L lines[2];
                for (int w = 0; w < 2; ++w) {
                    int t = w ? -1 : 1;
                    lines[w] = {2 * t * center[w][(k + 1) % 3], 2 * t * center[w][(k + 2) % 3],
                                -t * face * center[w][k]};
                }
                auto area = cube_face_integral::projected_area(lines, 2);
                lo += area.first;
                hi += area.second;
            }
        constexpr Int gc = cube_face_integral::coord_grid;
        out.at = {floor_ratio(lo * mass, 6 * gc * gc), ceil_ratio(hi * mass, 6 * gc * gc)};
        for (int w = 0; w < 2; ++w) {
            int r = w ? j : i;
            for (auto &g : out.deriv[w])
                g = {-8 * mass, 8 * mass};
            for (auto &g : out.hessian[w])
                g = {-1000000 * mass, 1000000 * mass};
            int best = 0;
            for (int k = 1; k < 3; ++k)
                if (absmin(n[w][k]) > absmin(n[w][best]))
                    best = k;
            for (int k = 0; k < 3; ++k) {
                if (absmin(n[w][k]) < mass / 4 || 2 * absmin(n[w][k]) < absmin(n[w][best]))
                    continue;
                int x = (k + 1) % 3, y = (k + 2) % 3, K = r ^ (k + 1), X = r ^ (x + 1), Y = r ^ (y + 1);
                I den = scalar(n[w][k], sg(r, k));
                auto rat = [&](int keep, int drop) {
                    auto t = linked_ratio(c[r][keep], c[r][drop], c[r][K]);
                    return I{t.lo, t.hi};
                };
                I u = scalar(rat(Y, X), sg(r, x) * sg(r, k)), v = scalar(rat(X, Y), sg(r, y) * sg(r, k));
                I A = sub(n[1 - w][x], mul(n[1 - w][k], u)), B = sub(n[1 - w][y], mul(n[1 - w][k], v));
                I gx = dmean(A, B), gy = dmean(B, A);
                auto pu = share(c[r][K], c[r][Y]), qu = share(c[r][K], c[r][X]);
                I P{pu.lo, pu.hi}, R{qu.lo, qu.hi};
                I a = scalar(A, sg(r, x) * sg(r, k)), b = scalar(B, sg(r, y) * sg(r, k));
                I f0 = divide(add(a, b), 2), f1 = sub(f0, mul(a, P)), f2 = sub(f0, mul(b, R)),
                  area = mul(P, R);
                // Integrate the clipped corner triangle; split crossed triangles in four.
                struct Tri {
                    I f[3], x[3], y[3];
                };
                Tri src{{f0, f1, f2},
                        {{mass / 2, mass / 2}, sub({mass / 2, mass / 2}, P), {mass / 2, mass / 2}},
                        {{mass / 2, mass / 2}, {mass / 2, mass / 2}, sub({mass / 2, mass / 2}, R)}};
                I mx{}, my{};
                auto integrate = [&](auto &&self, const Tri &t, I weight, int depth) -> void {
                    Int l = std::min({t.f[0].l, t.f[1].l, t.f[2].l}),
                        u = std::max({t.f[0].h, t.f[1].h, t.f[2].h});
                    if (l > 0 || u < 0) {
                        int s = l > 0 ? 1 : -1;
                        mx = add(mx, mul(weight, scalar(divide(add(add(t.x[0], t.x[1]), t.x[2]), 3), s)));
                        my = add(my, mul(weight, scalar(divide(add(add(t.y[0], t.y[1]), t.y[2]), 3), s)));
                        return;
                    }
                    if (!depth) {
                        Int X = std::max({absmax(t.x[0]), absmax(t.x[1]), absmax(t.x[2])}),
                            Y = std::max({absmax(t.y[0]), absmax(t.y[1]), absmax(t.y[2])});
                        mx = add(mx, mul(weight, {-X, X}));
                        my = add(my, mul(weight, {-Y, Y}));
                        return;
                    }
                    I f[6], xx[6], yy[6];
                    for (int j = 0; j < 3; ++j) {
                        f[j] = t.f[j];
                        xx[j] = t.x[j];
                        yy[j] = t.y[j];
                        f[3 + j] = divide(add(t.f[j], t.f[(j + 1) % 3]), 2);
                        xx[3 + j] = divide(add(t.x[j], t.x[(j + 1) % 3]), 2);
                        yy[3 + j] = divide(add(t.y[j], t.y[(j + 1) % 3]), 2);
                    }
                    constexpr int subtri[4][3] = {{0, 3, 5}, {1, 4, 3}, {2, 5, 4}, {3, 4, 5}};
                    for (auto const &xyz : subtri) {
                        Tri q;
                        for (int j = 0; j < 3; ++j) {
                            int l = xyz[j];
                            q.f[j] = f[l];
                            q.x[j] = xx[l];
                            q.y[j] = yy[l];
                        }
                        self(self, q, divide(weight, 4), depth - 1);
                    }
                };
                integrate(integrate, src, area, 1);
                gx = sub(gx, scalar(mx, sg(r, x) * sg(r, k)));
                gy = sub(gy, scalar(my, sg(r, y) * sg(r, k)));
                Norm cross{};
                for (int l = 0; l < 3; ++l) {
                    int y = (l + 1) % 3, z = (l + 2) % 3;
                    cross[l] = sub(mul(n[w][y], n[1 - w][z]), mul(n[w][z], n[1 - w][y]));
                }
                I length{};
                for (auto e : cross) {
                    length.l = std::max(length.l, absmin(e));
                    length.h = std::max(length.h, absmax(e));
                }
                I divisor = mul(add(scores(c[r][K]), scores(c[r][X])), add(scores(c[r][K]), scores(c[r][Y])));
                I density = divisor.l > mass / 512
                                ? divi(mul(scores(c[r][K]), add(scores(c[r][X]), scores(c[r][Y]))),
                                       mul(add(scores(c[r][K]), scores(c[r][X])),
                                           add(scores(c[r][K]), scores(c[r][Y]))))
                                : I{0, 2 * mass};
                int cur = 0;
                for (int col = 0; col < 4; ++col)
                    if (col != r && col != piv[r]) {
                        I ax{(col != X) * sg(r, x) * mass, (col != X) * sg(r, x) * mass},
                            ay{(col != Y) * sg(r, y) * mass, (col != Y) * sg(r, y) * mass};
                        if (col != K) {
                            ax = sub(ax, scalar(u, sg(r, k)));
                            ay = sub(ay, scalar(v, sg(r, k)));
                        }
                        I g = neg(divi(divide(add(mul(ax, gx), mul(ay, gy)), 2), den));
                        out.deriv[w][cur] = intersection(out.deriv[w][cur], g);
                        if (out.deriv[w][cur].l > out.deriv[w][cur].h)
                            throw std::logic_error("derivative intersection");
                        if (length.l > mass / 24) {
                            Norm vv{};
                            I dprod{};
                            for (int vaxis = 0; vaxis < 3; ++vaxis) {
                                Int a = col == (r ^ (vaxis + 1)) ? 0 : sg(r, vaxis) * mass / 2;
                                vv[vaxis] = {a, a};
                                dprod = add(dprod, mul(vv[vaxis], cross[vaxis]));
                            }
                            auto line = divi(sq(dprod), scalar(mul(sq(length), length), 12));
                            auto second = add(scalar(mul(divi(scalar(vv[k], sg(r, k)), den), g), -2),
                                              scalar(mul(divi(scalar(n[1 - w][k], sg(r, k)), den), line), 2));
                            auto left = add(mul(scalar(ax, sg(r, x) * sg(r, k)), src.x[1]),
                                            mul(scalar(ay, sg(r, y) * sg(r, k)), src.y[1]));
                            auto right = add(mul(scalar(ax, sg(r, x) * sg(r, k)), src.x[2]),
                                             mul(scalar(ay, sg(r, y) * sg(r, k)), src.y[2]));
                            // ax,ay are twice the projected variation. The opposite face supplies a second
                            // copy.
                            I squareMean =
                                positive(divide(add(add(sq(left), sq(right)), mul(left, right)), 12));
                            I lim = divi(scalar(mul(density, squareMean), 2), sq(den));
                            I signed_value{-lim.h, lim.h};
                            if (f1.l > 0 && f2.l > 0)
                                signed_value = lim;
                            if (f1.h < 0 && f2.h < 0)
                                signed_value = neg(lim);
                            second = add(second, scalar(signed_value, sg(r, k)));
                            auto &H = out.hessian[w][cur];
                            H = out.has[w][cur] ? intersection(H, second) : second;
                            out.has[w][cur] = true;
                            if (H.l > H.h)
                                throw std::logic_error("Hessian intersection");
                        }
                        ++cur;
                    }
            }
        }
        for (int j = 0; j < 16; ++j)
            out.corner[j] = corner_value(c, piv, h, j);
        if (table.size() > 120000)
            table.clear();
        table.emplace(key, out);
        return out;
    }
    struct State {
        PhysicalPair pairs[6];
        S dev{};
    };
    State prepare(const Matrix &c, const Pivots &piv) {
        State out;
        for (int h = 0; h < 6; ++h)
            out.pairs[h] = compute(c, piv, h);
        int a = 0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (j != i && j != piv[i]) {
                    auto s = scores(c[i][j]);
                    Int t = (s.h - s.l) / 2;
                    out.dev[a++] = {-t, t};
                }
        return out;
    }
    using Vec = std::array<Int, 7>;
    struct Interpolation {
        std::array<Vec, 256> corners{};
        I err[6][4]{};
        Interpolation(const State &s) {
            for (int h = 0; h < 6; ++h)
                for (int k = 0; k < 4; ++k) {
                    auto &p = s.pairs[h];
                    int a = 2 * edges[h][k / 2] + k % 2;
                    auto d = p.deriv[k / 2][k % 2], H = p.hessian[k / 2][k % 2];
                    Int radius = s.dev[a].h;
                    Int lip = ceil_ratio(Wide(d.h - d.l) * radius, 2 * mass);
                    I errh = mul(divide(sq({radius, radius}), 2), H);
                    err[h][k] = !p.has[k / 2][k % 2] ? I{-lip, lip}
                                                     : I{std::max(-lip, -std::max<Int>(0, errh.h)),
                                                         std::min(lip, -std::min<Int>(0, errh.l))};
                }
            for (int m = 0; m < 256; ++m) {
                I d[6];
                for (int h = 0; h < 6; ++h) {
                    auto [i, j] = edges[h];
                    d[h] = s.pairs[h].corner[((m >> (2 * i)) & 3) + 4 * ((m >> (2 * j)) & 3)];
                }
                for (int r = 0; r < 7; ++r) {
                    Int base = (r < 4 ? 6 : -6) * mass;
                    for (int h = 0; h < 6; ++h) {
                        int c = 3 * coeff(r, h);
                        base += c * (c >= 0 ? d[h].l : d[h].h);
                    }
                    corners[m][r] = base;
                }
            }
        }
        Vec worst(const std::array<int, 14> &w) const {
            std::array<int, 7> v;
            for (int r = 0; r < 7; ++r)
                v[r] = w[2 * r] + w[2 * r + 1];
            int best = -1;
            Wide low = 0;
            for (int m = 0; m < 256; ++m) {
                Wide z = 0;
                for (int r = 0; r < 7; ++r)
                    z += Wide(v[r]) * corners[m][r];
                if (best < 0 || z < low) {
                    low = z;
                    best = m;
                }
            }
            Vec z = corners[best];
            for (int h = 0; h < 6; ++h) {
                Int c = 0;
                for (int r = 0; r < 7; ++r)
                    c += v[r] * coeff(r, h);
                Int delta = 0;
                for (int k = 0; k < 4; ++k)
                    delta += c >= 0 ? err[h][k].l : err[h][k].h;
                for (int r = 0; r < 7; ++r)
                    z[r] += 3 * coeff(r, h) * delta;
            }
            return z;
        }
    };
    struct Atom {
        I value{};
        S g{};
    };
    static Atom addatom(const Atom &a, const Atom &b) {
        Atom z;
        z.value = add(a.value, b.value);
        for (int j = 0; j < 8; ++j)
            z.g[j] = add(a.g[j], b.g[j]);
        return z;
    }
    static Int eval(const Atom &a, const S &dev) {
        Int z = a.value.l;
        for (int j = 0; j < 8; ++j)
            z += mul(a.g[j], dev[j]).l;
        return z;
    }
    static std::array<Int, 128> bounds(const State &state) {
        static const auto multi = [] {
            std::array<std::array<int, 7>, 128> z{};
            for (int m = 1; m < 128; ++m) {
                auto &a = z[m];
                a = z[m & (m - 1)];
                int r = std::countr_zero(unsigned(m));
                a[6] += r < 4 ? 1 : -1;
                for (int h = 0; h < 6; ++h)
                    a[h] += coeff(r, h);
            }
            return z;
        }();
        std::array<Int, 128> lower{};
        for (int m = 1; m < 128; ++m) {
            Atom s;
            s.value = {6 * multi[m][6] * mass, 6 * multi[m][6] * mass};
            for (int h = 0; h < 6; ++h) {
                int c = 3 * multi[m][h];
                if (!c)
                    continue;
                auto &pp = state.pairs[h];
                s.value = add(s.value, scalar(pp.at, c));
                for (int i = 0; i < 2; ++i)
                    for (int j = 0; j < 2; ++j) {
                        int a = 2 * edges[h][i] + j;
                        s.g[a] = add(s.g[a], scalar(pp.deriv[i][j], c));
                    }
            }
            lower[m] = eval(s, state.dev);
        }
        return lower;
    }
};
