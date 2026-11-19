// Necessary halfspaces from the 364 extreme proper-cell deficit vectors.
struct DeficitCut {
    constexpr static int K = 23; // seven centered pairs, eight cap inequalities and eight triple identities
    using Action = std::array<int, K>;
    using Values = std::array<Int, K>;
    struct End {
        int e[4], t, delta[7], odd[4];
    };
    static const std::vector<End> &extremes() {
        static const auto vertices = [] {
            std::vector<End> out;
            for (int t : {0, 8452})
                for (int a = 1; a < 15; ++a)
                    for (int b = 1; b < 15; ++b)
                        if (a != b) {
                            End d{};
                            d.t = t;
                            int def[16]{};
                            def[a] = 79759;
                            def[b] = 2 * 58313 + t - 79759;
                            for (int j = 0; j < 4; ++j) {
                                int m = 1 << j;
                                d.e[j] = -58313 + (a & m ? def[a] : 0) + (b & m ? def[b] : 0);
                                d.odd[j] = def[15 ^ m] - def[m];
                            }
                            for (int r = 0; r < 7; ++r)
                                d.delta[r] = def[reps[r]] + def[15 ^ reps[r]];
                            out.push_back(d);
                        }
            return out;
        }();
        return vertices;
    }
    struct Scenario {
        Int body[K]{};
        I D[6]{};
        int e[4]{};
    };
    const MomentDomain &model;
    JointCenter angular;
    PureRows pure{};
    std::optional<ProbabilityColumns> columns;
    std::vector<Scenario> samples;
    static I microvalue(int x) {
        return {floor_ratio(Wide(x) * mass, micro), ceil_ratio(Wide(x) * mass, micro)};
    }
    static Int huber(I v, I range, Int rate) {
        Int hi = sq(v).h;
        if (v.h < range.l) {
            auto a = I{range.l, range.l};
            hi = mul(a, sub(scalar(v, 2), a)).h;
        } else if (v.l > range.h) {
            auto a = I{range.h, range.h};
            hi = mul(a, sub(scalar(v, 2), a)).h;
        }
        return ceil_ratio(Wide(rate) * hi, micro);
    }
    DeficitCut(const MomentDomain &m, const OddCaps &e) : model(m), angular(*m.dual) {
        I H[4];
        for (int j = 0; j < 4; ++j)
            H[j] = intersection(m.max_h[j], divi(e[j], m.ratios[j]));
        pure = section_coefficients(*m.odd, m.ratios, H);
        if (use_columns())
            columns.emplace(m, e, pure);
        for (auto &x : extremes()) {
            Scenario v;
            I cap[4], t = microvalue(x.t);
            for (int r = 0; r < 7; ++r)
                v.body[r] = 12 * (volume_target - microvalue(79759).h) +
                            6 * sub(microvalue(x.delta[r]), scalar(t, r < 4 ? 1 : -1)).l;
            for (int j = 0; j < 4; ++j) {
                cap[j] = microvalue(x.e[j]);
                v.e[j] = x.e[j];
                v.body[7 + 2 * j] = 6 * (cap[j].l - e[j].h);
                v.body[8 + 2 * j] = 6 * (e[j].l - cap[j].h);
            }
            for (int k = 0; k < 4; ++k) {
                I residual = add(microvalue(x.odd[k]), t);
                if (!columns)
                    for (int j = 0; j < 4; ++j)
                        residual = sub(residual, mul(pure[k][j], cap[j]));

                v.body[15 + 2 * k] = 6 * residual.l;
                v.body[16 + 2 * k] = -6 * residual.h;
            }
            for (int h = 0; h < 6; ++h) {
                auto [i, j] = edges[h];
                Int neg = huber(sub(cap[i], cap[j]), sub(e[i], e[j]), m.errors->rates[h][0]);
                Int pos = huber(add(cap[i], cap[j]), add(e[i], e[j]), m.errors->rates[h][1]);
                v.D[h] = {-neg, pos};
            }
            samples.push_back(v);
        }
    }
    int seeking = -1;
    Values score(const Action &w) {
        std::array<int, 14> wt{};
        for (int r = 0; r < 7; ++r)
            wt[2 * r] = w[r];
        auto z = angular.center(wt);
        Wide lowFace[16]{};
        int faceIndex[16]{};
        auto &lat = *model.lattice;
        int closest = -1;
        Wide angleLow = 0;
        JointCenter::Vec err{};
        for (int h = 0; h < 6; ++h) {
            int f = 0;
            for (int r = 0; r < 7; ++r)
                f += w[r] * DirectionPairs::coeff(r, h);
            Int delta = 0;
            for (auto x : lat.err[h])
                delta += f >= 0 ? x.l : x.h;
            for (int r = 0; r < 7; ++r)
                err[r] += 3 * DirectionPairs::coeff(r, h) * delta;
        }
        std::array<Wide, 256> work;
        std::array<int, 256> which;
        for (int id = 0; id < 256; ++id) {
            auto &extreme = lat.corners[id];
            Wide t = 0;
            for (int r = 0; r < 7; ++r)
                t += Wide(w[r]) * (extreme[r] + err[r] - 12 * volume_target);
            which[id] = id;
            work[id] = t;
            if (closest < 0 || t < angleLow) {
                closest = id;
                angleLow = t;
            }
        }
        for (int j = 0; j < 8; ++j) {
            int num = 128 >> j;
            Wide vals[2]{work[0], work[1]};
            int args[2]{which[0], which[1]};
            for (int k = 0; k < num; ++k) {
                for (int x = 0; x < 2; ++x) {
                    int v = 2 * k + x;
                    if (work[v] < vals[x]) {
                        vals[x] = work[v];
                        args[x] = which[v];
                    }
                }
                int v = 2 * k + (work[2 * k + 1] < work[2 * k]);
                work[k] = work[v];
                which[k] = which[v];
            }
            for (int x = 0; x < 2; ++x) {
                lowFace[2 * j + x] = vals[x];
                faceIndex[2 * j + x] = args[x];
            }
        }
        if (seeking >= 0) {
            closest = faceIndex[seeking];
            angleLow = lowFace[seeking];
        }
        if (seeking >= 0 || angleLow > JointCenter::dot(z, wt)) {
            z = lat.corners[closest];
            for (int r = 0; r < 7; ++r)
                z[r] += err[r] - 12 * volume_target;
        }
        int c[6]{};
        for (int h = 0; h < 6; ++h)
            for (int r = 0; r < 7; ++r)
                c[h] += w[r] * DirectionPairs::coeff(r, h);
        std::array<int, 4> oddCost;
        for (int i = 0; i < 4; ++i)
            oddCost[i] = 6 * (w[16 + 2 * i] - w[15 + 2 * i]);
        ProbabilityColumns::BySign sect{};
        if (columns)
            sect = columns->supports(oddCost);
        Int price[4][2]{};
        for (int j = 0; j < 4; ++j)
            for (int sign = 0; sign < 2; ++sign)
                for (int i = 0; i < 4; ++i)
                    price[j][sign] += Int(oddCost[i]) * sect[j][sign][i];
        Wide low = 0;
        int idx = -1;
        for (int j = 0; j < (int)samples.size(); ++j) {
            auto &v = samples[j];
            Wide a = 0;
            for (int k = 0; k < K; ++k)
                a += Wide(w[k]) * v.body[k];
            for (int h = 0; h < 6; ++h)
                a -= 3 * Wide(c[h]) * (c[h] > 0 ? v.D[h].h : v.D[h].l);
            Wide flow = 0;
            for (int k = 0; k < 4; ++k)
                flow += Wide(v.e[k]) * price[k][v.e[k] >= 0];
            a += floor_ratio(flow, micro);
            if (idx < 0 || a < low) {
                idx = j;
                low = a;
            }
        }
        auto &v = samples[idx];
        Values a;
        for (int j = 0; j < K; ++j)
            a[j] = v.body[j];
        for (int i = 0; i < 4; ++i) {
            Wide cost = 0;
            for (int j = 0; j < 4; ++j)
                cost += Wide(6) * v.e[j] * sect[j][v.e[j] >= 0][i];
            a[15 + 2 * i] -= ceil_ratio(cost, micro);
            a[16 + 2 * i] += floor_ratio(cost, micro);
        }
        for (int r = 0; r < 7; ++r) {
            a[r] += z[r];
            for (int h = 0; h < 6; ++h)
                a[r] -= 3 * DirectionPairs::coeff(r, h) * (c[h] > 0 ? v.D[h].h : v.D[h].l);
        }
        for (int k = 0; k < 8; ++k) {
            auto L = low + lowFace[2 * k], R = low + lowFace[2 * k + 1];
            conditional[2 * k] = L;
            conditional[2 * k + 1] = R;
            if (L > 0 && R > 0) {
                projection.lo[k] = mass;
                projection.hi[k] = 0;
            } else if (L > 0)
                projection.lo[k] = std::max(projection.lo[k], Int(L * mass / (L - R)));
            else if (R > 0)
                projection.hi[k] = std::min(projection.hi[k], mass - Int(R * mass / (R - L)));
        }
        return a;
    }
    struct Projection {
        std::array<Int, 8> lo{}, hi{};
        Projection() { hi.fill(mass); }
    };
    Projection projection;
    std::array<Wide, 16> conditional{};
    bool ok() { return attempt(-1, 15); }
    bool attempt(int face, int rounds) {
        seeking = face;
        Action w{};
        for (int j = 0; j < 7; ++j)
            w[j] = 1;
        std::vector<Values> batch;
        for (int i = 0; i < rounds; ++i) {
            auto a = score(w);
            Wide z = 0;
            for (int j = 0; j < K; ++j)
                z += Wide(w[j]) * a[j];
            if (z > 0)
                return true;
            batch.push_back(a);
            w = proposal<K>(batch);
            if (*std::max_element(w.begin(), w.end()) == 0)
                return false;
        }
        return false;
    }
};
