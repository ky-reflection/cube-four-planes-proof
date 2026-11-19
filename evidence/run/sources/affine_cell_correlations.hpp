// Necessary pair and marginal conditions for the four affine sign cuts.
struct PairInfo {
    Interval same, different;
};
struct CorrChecks {
    array<std::unordered_map<Ids, PairInfo, SubsetKeyHash>, 16> cache;
    static constexpr Int W = volume_grid, Q = (cap_numerator * W + cap_denominator - 1) / cap_denominator,
                         T = (8452 * W + cap_denominator - 1) / cap_denominator;
    template <class Parent> PairInfo pair(Parent &s, const Ids &ids, int i, int j) {
        int subset = (1 << i) | (1 << j);
        Ids key{};
        key[i] = ids[i];
        key[j] = ids[j];
        auto &tab = cache[subset];
        auto it = tab.find(key);
        if (it != tab.end())
            return it->second;
        auto &m = s.sectors(ids, subset);
        Interval D = {m[1 << i] + m[1 << j], W - m[0] - m[subset]};
        // At zero offsets opposite signs have the same volume. The centered mesh retains row direction
        // uncertainty.
        Int cl = 2 * diagonal_volume<2, Bound::Lower, true>(s.dom, ids, subset, 1 << i);
        Int cu = W - 2 * diagonal_volume<2, Bound::Lower, true>(s.dom, ids, subset, 0);
        auto p = s.dom[i].rows[ids[i]].cap, q = s.dom[j].rows[ids[j]].cap;
        Int sum_error = std::max(std::abs(p.lo + q.lo - W), std::abs(p.hi + q.hi - W));
        Int diff_error = std::max(std::abs(p.lo - q.hi), std::abs(p.hi - q.lo));
        Int lower_error = -floor_division(-Wide(sum_error) * sum_error, 3 * (345462 * W / 1000000));
        Int upper_error = -floor_division(-Wide(diff_error) * diff_error, 3 * (345462 * W / 1000000));
        D.lo = std::max(D.lo, cl - lower_error);
        D.hi = std::min(D.hi, cu + upper_error);
        PairInfo val{{W - D.hi, W - D.lo}, D};
        if (tab.size() > 100000)
            tab.clear();
        tab.emplace(key, val);
        return val;
    }
    static int rep(int k) {
        constexpr int masks[] = {1, 2, 4, 8, 3, 5, 9};
        return masks[k];
    }
    static auto supporting_directions() {
        std::vector<array<int, 4>> rows;
        for (int a = -2; a <= 2; ++a)
            for (int b = -2; b <= 2; ++b)
                for (int c = -2; c <= 2; ++c)
                    for (int d = -2; d <= 2; ++d) {
                        array<int, 4> r{a, b, c, d};
                        int ortho = 0, g = 0;
                        for (int v : r)
                            g = std::gcd(g, std::abs(v));
                        if (g != 1)
                            continue;
                        for (int k = 0; k < 7; ++k) {
                            int dot = 0;
                            for (int i = 0; i < 4; ++i)
                                dot += r[i] * (bit(rep(k), i) ? 1 : -1);
                            ortho += dot == 0;
                        }
                        if (ortho >= 3)
                            rows.push_back(r);
                    }
        return rows;
    }
    template <class Parent> unsigned feasible(Parent &s, const Ids &ids, unsigned possibles) {
        static const auto dirs = supporting_directions();
        array<array<PairInfo, 4>, 4> info{};
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < i; ++j)
                info[i][j] = pair(s, ids, i, j);
        std::array<Interval, 4> minors{};
        Interval normals[4][3], heights[4];
        constexpr Int G = graph_grid;
        for (int i = 0; i < 4; ++i) {
            auto &dom = s.dom[i];
            auto &r = dom.rows[ids[i]];
            normals[i][dom.k] = {2 * G * dom.sign, 2 * G * dom.sign};
            normals[i][dom.i] = {r.lo[1] - r.hi[2], r.hi[1] - r.lo[2]};
            normals[i][dom.j] = {r.lo[1] + r.lo[2], r.hi[1] + r.hi[2]};
            heights[i] = {-r.hi[0], -r.lo[0]};
        }
        for (int o = 0; o < 4; ++o) {
            int a = (o + 1) % 4, b = (o + 2) % 4, c = (o + 3) % 4;
            Interval m{};
            for (int x = 0; x < 3; ++x) {
                int y = (x + 1) % 3, z = (x + 2) % 3;
                auto u = product(normals[b][y], normals[c][z]), v = product(normals[b][z], normals[c][y]);
                auto d = product(normals[a][x], {u.lo - v.hi, u.hi - v.lo});
                m.lo += d.lo;
                m.hi += d.hi;
            }
            // Removing row o leaves a cyclic permutation of the three remaining rows.
            minors[o] = o % 2 ? negate(m) : m;
        }
        std::array<Int, 8> center{};
        for (int m = 0; m < 8; ++m)
            center[m] = diagonal_volume<4, Bound::Lower, true>(s.dom, ids, 15, m);
        for (int g = 0; g < 16; ++g)
            if (possibles >> g & 1) {
                array<Interval, 4> P;
                for (int i = 0; i < 4; ++i) {
                    auto r = s.dom[i].rows[ids[i]].cap;
                    P[i] = bit(g, i) ? r : Interval{W - r.hi, W - r.lo};
                }
                bool pos = true, neg = true;
                Int Bmax = 0;
                for (int i = 0; i < 4; ++i) {
                    auto m = bit(g, i) ? minors[i] : negate(minors[i]);
                    pos &= m.hi > 0;
                    neg &= m.lo < 0;
                    Interval mag{std::max({Int(0), m.lo, -m.hi}), std::max(m.hi, -m.lo)};
                    auto b = product(mag, bit(g, i) ? heights[i] : negate(heights[i]));
                    Bmax += b.hi;
                }
                if ((!pos && !neg) || Bmax < 0) {
                    possibles &= ~(1u << g);
                    continue;
                }
                Int tmin = 0, tmax = T;
                array<Interval, 7> Y{};
                for (int pass = 0; pass < 2 && tmin <= tmax; ++pass) {
                    for (int k = 0; k < 7; ++k) {
                        int r = rep(k);
                        Int l = 0, h = 0;
                        for (int i = 0; i < 4; ++i)
                            for (int j = 0; j < i; ++j) {
                                if (k < 4 ? (k == i || k == j) : bit(r, i) == bit(r, j))
                                    continue;
                                auto d = bit(g, i) == bit(g, j) ? info[i][j].different : info[i][j].same;
                                l += d.lo;
                                h += d.hi;
                            }
                        if (k < 4) {
                            Y[k] = {W - tmax - (h + 1) / 2, W - tmin - l / 2};
                            tmin = std::max(tmin, W - (h + 1) / 2 - 2 * Q);
                            tmax = std::min(tmax, W - l / 2);
                        } else {
                            Y[k] = {tmin - W + l / 2, tmax - W + (h + 1) / 2};
                            tmin = std::max(tmin, W - (h + 1) / 2);
                            tmax = std::min(tmax, W - l / 2 + 2 * Q);
                        }
                        int mask = (15 ^ g) ^ r;
                        Int shift = (k < 4 ? -tmax : 0);
                        Int delta = 0;
                        for (int i = 0; i < 4; ++i)
                            for (int j = 0; j < i; ++j) {
                                if (k < 4 ? (k == i || k == j) : bit(r, i) == bit(r, j))
                                    continue;
                                bool sum = (k >= 4);
                                Int a = sum ? P[i].lo + P[j].lo - W : P[i].lo - P[j].hi;
                                Int b = sum ? P[i].hi + P[j].hi - W : P[i].hi - P[j].lo;
                                Int e = std::max(std::abs(a), std::abs(b));
                                delta += -floor_division(-Wide(e) * e,
                                                         3 * ((sum ? 345462 : 379464) * W / 1000000));
                            }
                        Y[k].lo = std::max(Y[k].lo, 2 * center[mask > 7 ? 15 ^ mask : mask] + shift -
                                                        (delta + 1) / 2);
                        if (Y[k].lo > Y[k].hi || Y[k].lo > 2 * Q || Y[k].hi < 0) {
                            tmax = -1;
                            break;
                        }
                    }
                }
                if (tmin > tmax) {
                    possibles &= ~(1u << g);
                    continue;
                }
                Int room[7];
                for (int k = 0; k < 7; ++k)
                    room[k] = std::min({Q, 2 * Q - Y[k].lo, Y[k].hi});
                for (auto dir : dirs) {
                    int mean = 0;
                    Int left = 0, right = 0;
                    for (int i = 0; i < 4; ++i) {
                        mean += dir[i];
                        left += dir[i] * (2 * (dir[i] > 0 ? P[i].lo : P[i].hi) - W);
                    }
                    left += mean * (mean > 0 ? tmin : tmax);
                    for (int k = 0; k < 7; ++k) {
                        int dot = 0;
                        for (int i = 0; i < 4; ++i)
                            dot += dir[i] * (bit(rep(k), i) ? 1 : -1);
                        right += std::abs(dot) * room[k];
                    }
                    if (left > right) {
                        possibles &= ~(1u << g);
                        break;
                    }
                }
            }
        return possibles;
    }
};
