// Common nonnegative weights for volume comparisons over each cap simplex.
struct OffsetGame {
    using Vec = JointCenter::Vec;
    using Pay = JointCenter::Pay;
    using W = JointCenter::Weights;
    struct Endpoint {
        I t, delta[6]{};
        I odd[4]{};
        std::array<Int, 4> e{};
    };
    JointCenter angular;
    const DirectionPairs::Interpolation *lattice;
    std::optional<ProbabilityColumns> columns;
    PureRows P{};
    std::array<Endpoint, 81> verts{};
    int coeff[7][6]{};
    W warm{};
    OffsetGame(const MomentDomain &d, const OddCaps &box) : angular(*d.dual), lattice(d.lattice) {
        I H[4];
        for (int j = 0; j < 4; ++j)
            H[j] = intersection(d.max_h[j], divi(box[j], d.ratios[j]));
        P = section_coefficients(*d.odd, d.ratios, H);
        if (use_columns())
            columns.emplace(d, box, P);
        for (int r = 0; r < 7; ++r)
            for (int h = 0; h < 6; ++h) {
                auto [i, j] = edges[h];
                coeff[r][h] = r < 4 ? -(r != i && r != j) : int(((reps[r] >> i) ^ (reps[r] >> j)) & 1);
            }
        for (int v = 0; v < 81; ++v) {
            auto &p = verts[v];
            int vv = v;
            I B{};
            for (int i = 0; i < 4; ++i) {
                int t = vv % 3;
                vv /= 3;
                p.e[i] = t == 0 ? std::min<Int>(0, box[i].l) : t == 1 ? 0 : std::max<Int>(0, box[i].h);
                B = add(B, mul(divi(d.tree[i], d.ratios[i]), {p.e[i], p.e[i]}));
            }
            if (d.den.l > mass / 100000) {
                I z = positive(B);
                p.t = {0, divi(mul(sq(z), z), d.den).h};
            } else
                p.t = {0, 8 * mass};
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    p.odd[i] = add(p.odd[i], mul(P[i][j], {p.e[j], p.e[j]}));
            for (int h = 0; h < 6; ++h) {
                auto [i, j] = edges[h];
                Int e = p.e[i], f = p.e[j];
                Int a = ceil_ratio(Wide(d.errors->rates[h][0]) * sq({e - f, e - f}).h, micro);
                Int b = ceil_ratio(Wide(d.errors->rates[h][1]) * sq({e + f, e + f}).h, micro);
                p.delta[h] = {-a, b};
            }
        }
    }
    Pay score_old(W w, const std::vector<int> &ids) const {
        auto z = angular.center(w);
        if (lattice) {
            auto alt = lattice->worst(w);
            for (auto &x : alt)
                x -= 12 * volume_target;
            if (JointCenter::dot(alt, w) > JointCenter::dot(z, w))
                z = alt;
        }
        int sum[7]{}, diff[7]{}, c[6]{}, ct = 0, odd[4]{};
        for (int r = 0; r < 7; ++r) {
            sum[r] = w[2 * r] + w[2 * r + 1];
            diff[r] = w[2 * r + 1] - w[2 * r];
            for (int h = 0; h < 6; ++h)
                c[h] += coeff[r][h] * sum[r];
            ct += (r < 4 ? 1 : -1) * sum[r] + (r < 4 ? diff[r] : 0);
            for (int j = 0; j < 4; ++j)
                odd[j] += r < 4 ? 2 * (r == j) * diff[r] : ((reps[r] >> j & 1) ? -1 : 1) * diff[r];
        }
        std::array<int, 4> oddcost;
        std::copy(odd, odd + 4, oddcost.begin());
        ProbabilityColumns::BySign sect{};
        if (columns)
            sect = columns->supports(oddcost);
        int best = -1;
        Wide least = 0;
        Pay out{};
        for (auto id : ids) {
            auto &v = verts[id];
            Int t = ct > 0 ? v.t.h : 0;
            Int O[7], delta[6];
            for (int h = 0; h < 6; ++h)
                delta[h] = c[h] > 0 ? v.delta[h].h : v.delta[h].l;
            for (int i = 0; i < 4; ++i) {
                if (!columns) {
                    O[i] = odd[i] >= 0 ? v.odd[i].l : v.odd[i].h;
                    continue;
                }
                I calc{};
                for (int j = 0; j < 4; ++j) {
                    Int p = sect[j][v.e[j] >= 0][i];
                    calc = add(calc, mul({p, p}, {v.e[j], v.e[j]}));
                }
                O[i] = odd[i] >= 0 ? calc.l : calc.h;
            }
            // Express all signs in scale 12V; pair O appears with factor 3, so no half-grid division is
            // needed.
            for (int r = 4; r < 7; ++r) {
                O[r] = 0;
                for (int j = 0; j < 4; ++j)
                    O[r] += ((reps[r] >> j & 1) ? -1 : 1) * (v.e[j] + O[j]);
            }
            Pay val;
            Wide tmp = 0;
            for (int r = 0; r < 7; ++r) {
                Int y = z[r] - 6 * (r < 4 ? 1 : -1) * t;
                for (int h = 0; h < 6; ++h)
                    y -= 3 * coeff[r][h] * delta[h];
                Int diff = r < 4 ? 6 * (O[r] - t) : 3 * O[r];
                val[2 * r] = y - diff;
                val[2 * r + 1] = y + diff;
                tmp += Wide(val[2 * r]) * w[2 * r] + Wide(val[2 * r + 1]) * w[2 * r + 1];
            }
            if (best < 0 || tmp < least) {
                best = id;
                least = tmp;
                out = val;
            }
        }
        return out;
    }
    Pay score_flow(W w, const std::vector<int> &ids) const {
        auto z = angular.center(w);
        if (lattice) {
            auto alt = lattice->worst(w);
            for (auto &x : alt)
                x -= 12 * volume_target;
            if (JointCenter::dot(alt, w) > JointCenter::dot(z, w))
                z = alt;
        }
        int sum[7]{}, diff[7]{}, c[6]{}, ct = 0, odd[4]{};
        for (int r = 0; r < 7; ++r) {
            sum[r] = w[2 * r] + w[2 * r + 1];
            diff[r] = w[2 * r + 1] - w[2 * r];
            for (int h = 0; h < 6; ++h)
                c[h] += coeff[r][h] * sum[r];
            ct += (r < 4 ? 1 : -1) * sum[r] + (r < 4 ? diff[r] : 0);
            for (int j = 0; j < 4; ++j)
                odd[j] += r < 4 ? 2 * (r == j) * diff[r] : ((reps[r] >> j & 1) ? -1 : 1) * diff[r];
        }
        std::array<int, 4> oddcost;
        for (int i = 0; i < 4; ++i)
            oddcost[i] = odd[i] - ct;
        ProbabilityColumns::BySign sect{};
        if (columns)
            sect = columns->supports(oddcost);
        int best = -1;
        Wide least = 0;
        Pay out{};
        for (auto id : ids) {
            auto &v = verts[id];
            Int twiceT = 0, O[7], delta[6];
            for (int h = 0; h < 6; ++h)
                delta[h] = c[h] > 0 ? v.delta[h].h : v.delta[h].l;
            for (int i = 0; i < 4; ++i) {
                if (!columns) {
                    O[i] = odd[i] >= 0 ? v.odd[i].l : v.odd[i].h;
                    continue;
                }
                I calc{};
                for (int j = 0; j < 4; ++j) {
                    Int p = sect[j][v.e[j] >= 0][i];
                    calc = add(calc, mul({p, p}, {v.e[j], v.e[j]}));
                }
                O[i] = oddcost[i] >= 0 ? calc.l : calc.h;
                twiceT += O[i] - v.e[i];
            }
            // Express all signs in scale 12V; pair O appears with factor 3, so no half-grid division is
            // needed.
            for (int r = 4; r < 7; ++r) {
                O[r] = 0;
                for (int j = 0; j < 4; ++j)
                    O[r] += ((reps[r] >> j & 1) ? -1 : 1) * (v.e[j] + O[j]);
            }
            Pay val;
            Wide tmp = 0;
            for (int r = 0; r < 7; ++r) {
                Int y = z[r] - 3 * (r < 4 ? 1 : -1) * twiceT;
                for (int h = 0; h < 6; ++h)
                    y -= 3 * coeff[r][h] * delta[h];
                Int diff = r < 4 ? 6 * O[r] - 3 * twiceT : 3 * O[r];
                val[2 * r] = y - diff;
                val[2 * r + 1] = y + diff;
                tmp += Wide(val[2 * r]) * w[2 * r] + Wide(val[2 * r + 1]) * w[2 * r + 1];
            }
            if (best < 0 || tmp < least) {
                best = id;
                least = tmp;
                out = val;
            }
        }
        return out;
    }
    Pay scenario(W w, const std::vector<int> &ids) const {
        if (mode != ProjectionMode::flow)
            return score_old(w, ids);
        auto a = score_flow(w, ids);
        Wide va = 0;
        for (int i = 0; i < 14; ++i)
            va += Wide(w[i]) * a[i];
        if (va > 0)
            return a;
        auto b = score_old(w, ids);
        Wide vb = 0;
        for (int i = 0; i < 14; ++i)
            vb += Wide(w[i]) * b[i];
        return va >= vb ? a : b;
    }
    bool prove(const std::vector<int> &ids, int rounds) {
        W w = warm;
        if (*std::max_element(w.begin(), w.end()) == 0)
            w.fill(1);
        std::vector<Pay> constraints;
        for (int j = 0; j < rounds; ++j) {
            auto z = scenario(w, ids);
            Wide value = 0;
            for (int r = 0; r < 14; ++r)
                value += Wide(z[r]) * w[r];
            if (value > 0) {
                warm = w;
                return true;
            }
            constraints.push_back(z);
            w = proposal<14>(constraints);
            if (*std::max_element(w.begin(), w.end()) == 0)
                return false;
        }
        return false;
    }
    bool all() {
        std::vector<int> ids(81);
        for (int i = 0; i < 81; ++i)
            ids[i] = i;
        if (prove(ids, 9))
            return true;
        constexpr int p3[]{1, 3, 9, 27};
        for (int orth = 0; orth < 16; ++orth) {
            std::array<int, 4> s{};
            for (int i = 0; i < 4; ++i)
                s[i] = p3[i] * ((orth >> i & 1) ? 1 : -1);
            ids.clear();
            for (int m = 0; m < 16; ++m) {
                int ix = 40;
                for (int j = 0; j < 4; ++j)
                    if (m >> j & 1)
                        ix += s[j];
                ids.push_back(ix);
            }
            if (prove(ids, 17))
                continue;
            int order[]{0, 1, 2, 3};
            do {
                ids = {40};
                for (int j : order)
                    ids.push_back(ids.back() + s[j]);
                if (!prove(ids, 20))
                    return false;
            } while (std::next_permutation(order, order + 4));
        }
        return true;
    }
};
