// Weighted centers and sharp pairwise transfer bounds shared by the affine covers.
struct AngularAtoms {
    RowControls rows[4];
    AffineOverlap atoms[18];
    int coefficients[7][18]{};
    Int constants[7]{};
    bool exclude(const std::array<I, 4> &e, I t, const I (&corr)[6], const std::array<Linear4, 7> &odd,
                 const std::array<I, 7> &off) const {
        using A = cube_affine_box::Linear<12>;
        std::array<A, 14> forms{};
        static const auto delta = [] {
            std::array<std::array<int, 6>, 7> z{};
            for (int r = 0; r < 7; ++r)
                for (int j = 0; j < 6; ++j) {
                    auto [a, b] = edges[j];
                    z[r][j] = r < 4 ? -(r != a && r != b) : int(((reps[r] >> a) ^ (reps[r] >> b)) & 1);
                }
            return z;
        }();
        auto mid = [](I v) { return (v.l + v.h) / 2; };
        for (int r = 0; r < 7; ++r)
            for (int s = 0; s < 2; ++s) {
                auto &f = forms[2 * r + s];
                int sign = 2 * s - 1;
                f.center = constants[r] - 6 * (r < 4 ? 1 : -1) * mid(t) + 6 * sign * mid(off[r]);
                for (int j = 0; j < 6; ++j)
                    f.center -= 3 * delta[r][j] * mid(corr[j]);
                for (int j = 0; j < 18; ++j) {
                    int c = coefficients[r][j];
                    f.center += c * atoms[j].center;
                    for (int k = 0; k < 8; ++k)
                        f.slope[k] += c * atoms[j].slope[k];
                }
                for (int k = 0; k < 4; ++k) {
                    auto a = 6 * sign * mid(odd[r][k]);
                    f.center += floor_ratio(Wide(a) * mid(e[k]), mass);
                    f.slope[8 + k] = floor_ratio(Wide(a) * (e[k].h - e[k].l), 2 * mass);
                }
            }
        auto w = cube_affine_box::propose_linear_weights<12, 14>(forms, 12 * volume_target);
        std::array<Int, 7> sum{}, diff{};
        Int wt = 0;
        for (int r = 0; r < 7; ++r) {
            sum[r] = w[2 * r] + w[2 * r + 1];
            diff[r] = w[2 * r + 1] - w[2 * r];
            wt += sum[r];
        }
        if (!wt)
            return false;
        Wide base = 0;
        Int coef[18]{};
        for (int r = 0; r < 7; ++r) {
            base += Wide(sum[r]) * (constants[r] - 12 * volume_target);
            for (int k = 0; k < 18; ++k)
                coef[k] += sum[r] * coefficients[r][k];
        }
        Wide rowLower = 0, affine = 0;
        std::array<Wide, 8> slope{};
        for (int i = 0; i < 4; ++i) {
            Int best = std::numeric_limits<Int>::max();
            for (int k = 0; k < 9; ++k)
                if (rows[i].denominator[k]) {
                    Wide num = 0;
                    for (int a = 0; a < 3; ++a)
                        num += Wide(coef[3 * i + a]) * rows[i].numerator[a][k];
                    best = std::min(best, floor_ratio(num * mass, rows[i].denominator[k]));
                }
            rowLower += best;
        }
        for (int k = 0; k < 18; ++k) {
            Int c = coef[k];
            auto &a = atoms[k];
            Int v = a.center + (c < 0 ? a.error.hi : a.error.lo);
            affine += Wide(c) * v;
            for (int j = 0; j < 8; ++j)
                slope[j] += Wide(c) * a.slope[j];
            if (k >= 12) {
                Wide l = Wide(c) * v;
                for (auto x : a.slope)
                    l -= std::abs(c * x);
                rowLower += l;
            }
        }
        for (auto x : slope)
            affine -= (x < 0 ? -x : x);
        base += std::max(rowLower, affine);
        Int ct = 0;
        for (int r = 0; r < 7; ++r)
            ct += sum[r] * (r < 4 ? 1 : -1);
        base -= Wide(6) * ct * (ct > 0 ? t.h : t.l);
        for (int h = 0; h < 6; ++h) {
            Int a = 0;
            for (int r = 0; r < 7; ++r)
                a += sum[r] * delta[r][h];
            base -= Wide(3) * a * (a > 0 ? corr[h].h : corr[h].l);
        }
        for (int j = 0; j < 4; ++j) {
            Wide l = 0, h = 0;
            for (int r = 0; r < 7; ++r) {
                auto v = odd[r][j];
                Int c = diff[r];
                l += Wide(c) * (c < 0 ? v.h : v.l);
                h += Wide(c) * (c < 0 ? v.l : v.h);
            }
            auto bounds = std::min({l * e[j].l, l * e[j].h, h * e[j].l, h * e[j].h});
            base += 6 * floor_ratio(bounds, mass);
        }
        for (int r = 0; r < 7; ++r) {
            Int c = diff[r];
            base += 6 * Wide(c) * (c < 0 ? off[r].h : off[r].l);
        }
        return base > 0;
    }
};

struct WeightedEndpoint {
    Int base, coeff[6];
};
struct CenterSupport {
    int h[6]{};
    Int base = 0;
    std::vector<WeightedEndpoint> verts;
};
static const auto convex_support = [] {
    std::array<CenterSupport, 128> output{};
    for (int m = 1; m < 127; ++m) {
        auto &o = output[m];
        int tcount = 0;
        for (int r = 0; r < 7; ++r)
            if (m >> r & 1) {
                tcount += r < 4 ? 1 : -1;
                for (int h = 0; h < 6; ++h) {
                    auto [i, j] = edges[h];
                    o.h[h] += r < 4 ? -(r != i && r != j) : int(((reps[r] >> i) ^ (reps[r] >> j)) & 1);
                }
            }
        std::vector<WeightedEndpoint> v;
        for (Int t : {0, 8452})
            for (int a = 1; a < 15; ++a)
                for (int b = 1; b < 15; ++b)
                    if (a != b) {
                        Int rest = 2 * 58313 + t - 79759, e[4];
                        for (int i = 0; i < 4; ++i)
                            e[i] = -58313 + 79759 * (a >> i & 1) + rest * (b >> i & 1);
                        WeightedEndpoint p{};
                        p.base = 12 * Int(std::popcount(unsigned(m))) * 79759 + 6 * t * tcount;
                        for (int r = 0; r < 7; ++r)
                            if (m >> r & 1) {
                                int s = reps[r];
                                p.base -= 6 * (79759 * ((s == a) || ((15 ^ s) == a)) +
                                               rest * ((s == b) || ((15 ^ s) == b)));
                            }
                        o.base = std::max(o.base, p.base);
                        for (int h = 0; h < 6; ++h) {
                            auto [i, j] = edges[h];
                            Int z = e[i] + (o.h[h] >= 0 ? e[j] : -e[j]);
                            p.coeff[h] = 3 * std::abs(o.h[h]) * z * z;
                        }
                        v.push_back(p);
                    }
        for (auto &p : v) {
            Int possible = p.base;
            for (auto b : p.coeff)
                possible += (b + micro - 1) / micro;
            if (possible < o.base)
                continue;
            auto dominates = [](auto &u, auto &w) {
                if (u.base < w.base)
                    return false;
                for (int h = 0; h < 6; ++h)
                    if (u.coeff[h] < w.coeff[h])
                        return false;
                return true;
            };
            bool dominated = false;
            for (auto &u : o.verts)
                if (dominates(u, p)) {
                    dominated = true;
                    break;
                }
            if (!dominated) {
                std::erase_if(o.verts, [&](auto &u) { return dominates(p, u); });
                o.verts.push_back(p);
            }
        }
    }
    return output;
}();
bool excludes(const AngularErrors &errors, int m, Int lower) {
    auto &support = convex_support[m];
    Int L = floor_ratio(Wide(lower) * micro, mass);
    if (L <= support.base)
        return false;
    for (auto &v : support.verts) {
        Wide bound = Wide(v.base) * micro * micro;
        for (int h = 0; h < 6; ++h)
            bound += Wide(v.coeff[h]) * errors.rates[h][support.h[h] >= 0];
        if (Wide(L) * micro * micro <= bound)
            return false;
    }
    return true;
}
