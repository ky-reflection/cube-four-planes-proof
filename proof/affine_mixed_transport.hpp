// Model methods for the signed cell volumes, vertex inequalities and cap
// simplices.
bool simplex_cover(const Caps &cap) const {
    if (!haveLattice || !haveOdd)
        return false;
    using Pay = std::array<Int, 14>;
    using Weight = std::array<int, 14>;
    R B{};
    for (int j = 0; j < 4; ++j)
        B = B + lambda[j] * intersect(hlimit[j], cap[j] / density[j]);
    if (B.hi < 0)
        return true;
    B = nonnegative(B);
    cube_section::Column cols[4];
    for (int j = 0; j < 4; ++j) {
        R rr{0, ceil_ratio(Wide(224) * U, 1000)};
        if (den.lo > U / 32768)
            rr = intersect(rr, (square(B) / den) * (lambda[j] / density[j]));
        cube_section::Triple lo{}, hi{};
        int z = 0;
        for (int i = 0; i < 4; ++i)
            if (i != j) {
                lo[z] = prob[i][j].lo;
                hi[z++] = prob[i][j].hi;
            }
        cols[j].init(U, lo, hi, rr.lo, rr.hi);
        if (cols[j].vertices.empty())
            return true;
    }
    struct End {
        Int e[4]{};
        R delta[6]{};
    };
    std::array<End, 81> ends{};
    Int rate[6][2];
    for (int h = 0; h < 6; ++h) {
        Int bot = std::max(Int(499) * U / 1000, D0[h].lo),
            top = std::min(ceil_ratio(Wide(654538) * U, 1000000), D0[h].hi);
        if (bot > top)
            return true;
        rate[h][0] = ceil_ratio(Wide(U) * 1000000, 3 * bot);
        rate[h][1] = ceil_ratio(Wide(U) * 1000000, 3 * (U - top));
    }
    for (int m = 0; m < 81; ++m) {
        auto &p = ends[m];
        int t = m;
        for (int i = 0; i < 4; ++i) {
            int d = t % 3;
            t /= 3;
            p.e[i] = d == 0 ? std::min(Int(0), cap[i].lo) : d == 1 ? 0 : std::max(Int(0), cap[i].hi);
        }
        for (int h = 0; h < 6; ++h) {
            auto [i, j] = pairid[h];
            p.delta[h] = {-ceil_ratio(Wide(rate[h][0]) * square(pt(p.e[i] - p.e[j])).hi, 1000000),
                          ceil_ratio(Wide(rate[h][1]) * square(pt(p.e[i] + p.e[j])).hi, 1000000)};
        }
    }
    auto score = [&](const Weight &w, const std::vector<int> &ids) {
        int sums[7]{}, sgn[7]{}, coeff[6]{}, odds[4]{}, ct = 0;
        Int err[7]{};
        for (int r = 0; r < 7; ++r) {
            sums[r] = w[2 * r + 1] + w[2 * r];
            sgn[r] = w[2 * r + 1] - w[2 * r];
            for (int h = 0; h < 6; ++h)
                coeff[h] += co[1 << r][h] * sums[r];
            ct += (r < 4 ? 1 : -1) * sums[r] + (r < 4 ? sgn[r] : 0);
            for (int i = 0; i < 4; ++i)
                odds[i] += (r < 4 ? 2 * (r == i) : ((rep[r] >> i & 1) ? -1 : 1)) * sgn[r];
        }
        for (int i = 0; i < 4; ++i)
            odds[i] -= ct;
        for (int h = 0; h < 6; ++h) {
            Int error = coeff[h] >= 0 ? latticeErr[h].lo : latticeErr[h].hi;
            for (int r = 0; r < 7; ++r)
                err[r] += 3 * co[1 << r][h] * error;
        }
        int at = 0;
        Wide low = 0;
        for (int mask = 0; mask < 256; ++mask) {
            Wide z = 0;
            for (int r = 0; r < 7; ++r)
                z += Wide(sums[r]) * latticeVertex[mask][r];
            if (!mask || z < low) {
                low = z;
                at = mask;
            }
        }
        Int col[4][2][4]{};
        for (int j = 0; j < 4; ++j)
            for (int pos = 0; pos < 2; ++pos) {
                std::array<int, 3> cost;
                int n = 0;
                for (int i = 0; i < 4; ++i)
                    if (i != j)
                        cost[n++] = (pos ? 1 : -1) * odds[i];
                auto v = cols[j].lowest(cost);
                n = 0;
                for (int i = 0; i < 4; ++i)
                    if (i != j)
                        col[j][pos][i] = v[n++];
            }
        Pay answer{};
        Wide least = 0;
        bool first = true;
        for (int id : ids) {
            auto &p = ends[id];
            Int O[7]{}, twiceT = 0;
            for (int i = 0; i < 4; ++i) {
                Wide v = 0;
                for (int j = 0; j < 4; ++j)
                    v += Wide(p.e[j]) * col[j][p.e[j] >= 0][i];
                O[i] = odds[i] >= 0 ? floor_ratio(v, U) : ceil_ratio(v, U);
                twiceT += O[i] - p.e[i];
            }
            for (int r = 4; r < 7; ++r)
                for (int i = 0; i < 4; ++i)
                    O[r] += ((rep[r] >> i & 1) ? -1 : 1) * (p.e[i] + O[i]);
            Pay out;
            Wide obj = 0;
            for (int r = 0; r < 7; ++r) {
                Int z = latticeVertex[at][r] + err[r] - 12 * Q - 3 * (r < 4 ? 1 : -1) * twiceT;
                for (int h = 0; h < 6; ++h)
                    z -= 3 * co[1 << r][h] * (coeff[h] > 0 ? p.delta[h].hi : p.delta[h].lo);
                Int odd = r < 4 ? 6 * O[r] - 3 * twiceT : 3 * O[r];
                out[2 * r] = z - odd;
                out[2 * r + 1] = z + odd;
                obj += Wide(out[2 * r]) * w[2 * r] + Wide(out[2 * r + 1]) * w[2 * r + 1];
            }
            if (first || obj < least) {
                least = obj;
                answer = out;
                first = false;
            }
        }
        return answer;
    };
    Weight warm{};
    warm.fill(1);
    auto prove = [&](const std::vector<int> &ids, int rounds) {
        auto w = warm;
        std::vector<Pay> active;
        for (int j = 0; j < rounds; ++j) {
            auto z = score(w, ids);
            Wide value = 0;
            for (int i = 0; i < 14; ++i)
                value += Wide(w[i]) * z[i];
            if (value > 0) {
                warm = w;
                return true;
            }
            active.push_back(z);
            w = proposal<14>(active);
            if (*std::max_element(w.begin(), w.end()) == 0)
                break;
        }
        return false;
    };
    std::vector<int> ids;
    constexpr int pow3[]{1, 3, 9, 27};
    for (int m = 0; m < 16; ++m) {
        int v = 0;
        for (int j = 0; j < 4; ++j)
            if (m >> j & 1)
                v += 2 * pow3[j];
        ids.push_back(v);
    }
    if (prove(ids, 10))
        return true;
    for (int o = 0; o < 16; ++o) {
        int axes[4], order[]{0, 1, 2, 3};
        for (int j = 0; j < 4; ++j)
            axes[j] = pow3[j] * ((o >> j & 1) ? 1 : -1);
        ids.clear();
        for (int m = 0; m < 16; ++m) {
            int v = 40;
            for (int j = 0; j < 4; ++j)
                if (m >> j & 1)
                    v += axes[j];
            ids.push_back(v);
        }
        if (prove(ids, 14))
            continue;
        do {
            ids = {40};
            for (int j : order)
                ids.push_back(ids.back() + axes[j]);
            if (!prove(ids, 18))
                return false;
        } while (std::next_permutation(order, order + 4));
    }
    return true;
}

bool vertex_projection(const Caps &caps, Box &b, bool &shorter) const {
    constexpr int K = 57;
    using Action = std::array<int, K>;
    using Values = std::array<Int, K>;
    if (!haveOdd || !haveLattice)
        return false;
    struct Scenario {
        Values v{};
        int e[4]{};
        R delta[6]{};
    };
    constexpr Int micro = 1000000;
    auto mic = [&](Int x) { return R{floor_ratio(Wide(x) * U, micro), ceil_ratio(Wide(x) * U, micro)}; };
    R H[4], bias[4], sumB{};
    for (int j = 0; j < 4; ++j) {
        H[j] = intersect(hlimit[j], caps[j] / density[j]);
        bias[j] = lambda[j] / density[j];
        sumB = sumB + lambda[j] * H[j];
    }
    if (sumB.hi < 0)
        return true;
    sumB = nonnegative(sumB);
    R strength{0, U};
    if (den.lo > U / 32768)
        strength = square(sumB) / den;
    cube_section::Column columns[4];
    for (int j = 0; j < 4; ++j) {
        R rr{0, ceil_ratio(Wide(224) * U, 1000)};
        if (den.lo > U / 32768)
            rr = intersect(rr, strength * bias[j]);
        cube_section::Triple lo{}, hi{};
        int q = 0;
        for (int i = 0; i < 4; ++i)
            if (i != j) {
                lo[q] = prob[i][j].lo;
                hi[q++] = prob[i][j].hi;
            }
        columns[j].init(U, lo, hi, rr.lo, rr.hi);
        if (columns[j].vertices.empty())
            return true;
    }
    struct Facet {
        int j, k, sg;
        Int dir[2];
    };
    std::array<Facet, 32> facets{};
    for (int j = 0; j < 4; ++j) {
        int k = type.axis[j], x = (k + 1) % 3, y = (k + 2) % 3;
        for (int v = 0; v < 8; ++v) {
            int X = 2 * ((v >> x) & 1) - 1, Y = 2 * ((v >> y) & 1) - 1, Z = 2 * ((v >> k) & 1) - 1,
                sg = (type.bits[j] >> v & 1) ? 1 : -1;
            int axis = 2 * j + (X != Y);
            auto &f = facets[8 * j + v];
            f = {j, axis, sg, {}};
            for (int side = 0; side < 2; ++side) {
                Int x = side ? b[axis].hi : b[axis].lo;
                f.dir[side] = -3 * sg * (Y * x + Z * type.orient[j] * U);
            }
        }
    }
    Int rate[6][2];
    for (int h = 0; h < 6; ++h) {
        auto [i, j] = pairid[h];
        Int low = std::max(Int(499) * U / 1000, D0[h].lo),
            high = std::min(ceil_ratio(Wide(654538) * U, micro), D0[h].hi);
        if (low > high)
            return true;
        rate[h][0] = ceil_ratio(Wide(U) * micro, 3 * low);
        rate[h][1] = ceil_ratio(Wide(U) * micro, 3 * (U - high));
        if (type.face[i] && type.face[j]) {
            auto &p = pairs[h];
            if (type.axis[i] == type.axis[j]) {
                rate[h][0] = 0;
                if (p.ac[0].lo > U / 4096)
                    rate[h][1] = std::min(rate[h][1], ceil_ratio(Wide(U) * micro, p.ac[0].lo));
            } else if (p.ac[2].lo > U / 4096) {
                R f = p.ac[0] + p.ac[1];
                if (2 * U + f.lo > U / 4096)
                    rate[h][0] = std::min(rate[h][0], ceil_ratio(Wide(U) * micro, 2 * U + f.lo));
                if (2 * U - f.hi > U / 4096)
                    rate[h][1] = std::min(rate[h][1], ceil_ratio(Wide(U) * micro, 2 * U - f.hi));
            }
        }
    }
    for (int h = 0; h < 6; ++h)
        if (transfer[h].valid)
            for (int f = 0; f < 2; ++f)
                rate[h][f] = std::min(rate[h][f], transfer[h].rate[f]);
    auto huber = [&](R s, R interval, Int r) {
        Int val = square(s).hi;
        if (s.hi < interval.lo)
            val = (pt(interval.lo) * (times(s, 2) - pt(interval.lo))).hi;
        else if (s.lo > interval.hi)
            val = (pt(interval.hi) * (times(s, 2) - pt(interval.hi))).hi;
        return ceil_ratio(Wide(val) * r, micro);
    };
    std::vector<Scenario> data;
    data.reserve(364);
    for (int t : {0, 8452})
        for (int a = 1; a < 15; ++a)
            for (int b = 1; b < 15; ++b)
                if (a != b) {
                    Scenario v;
                    int deficits[16]{};
                    deficits[a] = 79759;
                    deficits[b] = 2 * 58313 + t - 79759;
                    R B{};
                    for (int j = 0; j < 4; ++j) {
                        int m = 1 << j;
                        v.e[j] = -58313 + (a & m ? deficits[a] : 0) + (b & m ? deficits[b] : 0);
                        R e = mic(v.e[j]);
                        v.v[7 + 2 * j] = 6 * (e.lo - caps[j].hi);
                        v.v[8 + 2 * j] = 6 * (caps[j].lo - e.hi);
                        R o = mic(deficits[15 ^ m] - deficits[m] + t);
                        v.v[15 + 2 * j] = 6 * o.lo;
                        v.v[16 + 2 * j] = -6 * o.hi;
                        B = B + bias[j] * e;
                    }
                    for (int r = 0; r < 7; ++r) {
                        int m = rep[r];
                        v.v[r] = 6 * mic(deficits[m] + deficits[15 ^ m] - (r < 4 ? 1 : -1) * t).lo;
                    }
                    v.v[23] = -B.hi;
                    v.v[24] = den.lo > U / 32768 ? 6 * (mic(t) - strength * B).lo : 0;
                    for (int h = 0; h < 6; ++h) {
                        auto [i, j] = pairid[h];
                        v.delta[h] = {-huber(mic(v.e[i] - v.e[j]), caps[i] - caps[j], rate[h][0]),
                                      huber(mic(v.e[i] + v.e[j]), caps[i] + caps[j], rate[h][1])};
                    }
                    for (int k = 0; k < 32; ++k) {
                        auto f = facets[k];
                        v.v[25 + k] = 6 * times(mic(v.e[f.j]) / density[f.j], f.sg).lo;
                    }
                    data.push_back(v);
                }
    const Box original = b;
    Action w{};
    for (int r = 0; r < 7; ++r)
        w[r] = 1;
    std::vector<Values> batch;
    for (int step = 0; step < 18; ++step) {
        int c[6]{}, cost[4]{};
        Int err[7]{};
        for (int h = 0; h < 6; ++h) {
            for (int r = 0; r < 7; ++r)
                c[h] += co[1 << r][h] * w[r];
            Int d = c[h] >= 0 ? latticeErr[h].lo : latticeErr[h].hi;
            for (int r = 0; r < 7; ++r)
                err[r] += 3 * co[1 << r][h] * d;
        }
        for (int i = 0; i < 4; ++i)
            cost[i] = 6 * (w[16 + 2 * i] - w[15 + 2 * i]);
        Int col[4][2][4]{}, price[4][2]{};
        for (int j = 0; j < 4; ++j)
            for (int sign = 0; sign < 2; ++sign) {
                std::array<int, 3> z{};
                int n = 0;
                for (int i = 0; i < 4; ++i)
                    if (i != j)
                        z[n++] = cost[i] * (sign ? 1 : -1);
                auto p = columns[j].lowest(z);
                n = 0;
                auto &v = col[j][sign];
                for (int i = 0; i < 4; ++i)
                    if (i != j)
                        v[i] = p[n++];
                for (int i = 0; i < 4; ++i)
                    price[j][sign] += Int(cost[i]) * v[i];
            }
        Wide least = 0;
        int at = -1;
        for (int id = 0; id < (int)data.size(); ++id) {
            auto &v = data[id];
            Wide value = 0;
            for (int k = 0; k < K; ++k)
                value += Wide(w[k]) * v.v[k];
            for (int h = 0; h < 6; ++h)
                value -= 3 * Wide(c[h]) * (c[h] > 0 ? v.delta[h].hi : v.delta[h].lo);
            Wide flow = 0;
            for (int j = 0; j < 4; ++j)
                flow += Wide(v.e[j]) * price[j][v.e[j] >= 0];
            value += floor_ratio(flow, micro);
            if (at < 0 || value < least) {
                least = value;
                at = id;
            }
        }
        auto &v = data[at];
        Values rows = v.v;
        Wide base = 0, increment[8]{};
        for (int k = 0; k < 32; ++k)
            if (w[25 + k]) {
                auto f = facets[k];
                base += Wide(w[25 + k]) * f.dir[0];
                increment[f.k] += Wide(w[25 + k]) * (f.dir[1] - f.dir[0]);
            }
        Wide affine[256];
        affine[0] = base;
        for (int m = 1; m < 256; ++m)
            affine[m] = affine[m & (m - 1)] + increment[std::countr_zero(unsigned(m))];
        Wide angles[256], minAngle = 0;
        int best = -1;
        for (int m = 0; m < 256; ++m) {
            Wide z = affine[m];
            for (int r = 0; r < 7; ++r)
                z += Wide(w[r]) * (latticeVertex[m][r] + err[r] - 12 * Q);
            angles[m] = z;
            if (best < 0 || z < minAngle) {
                best = m;
                minAngle = z;
            }
        }
        if (least + minAngle > 0)
            return true;
        for (int k = 0; k < 8; ++k) {
            int count = 128 >> k;
            Wide x = angles[0], y = angles[1];
            for (int u = 0; u < count; ++u) {
                x = std::min(x, angles[2 * u]);
                y = std::min(y, angles[2 * u + 1]);
                angles[u] = std::min(angles[2 * u], angles[2 * u + 1]);
            }
            x += least;
            y += least;
            Int width = original[k].hi - original[k].lo;
            if (x > 0 && y > 0)
                return true;
            constexpr Int gridstep = U / 4096;
            Int slice = std::max(gridstep, (width / 8 / gridstep) * gridstep);
            if (x > 0) {
                Int size = Int(x * width / (x - y));
                b[k].lo = std::max(b[k].lo, original[k].lo + (size / slice) * slice);
            }
            if (y > 0) {
                Int size = Int(y * width / (y - x));
                b[k].hi = std::min(b[k].hi, original[k].hi - (size / slice) * slice);
            }
            if (b[k].lo > b[k].hi)
                return true;
            shorter |= b[k].lo > original[k].lo || b[k].hi < original[k].hi;
        }
        for (int r = 0; r < 7; ++r) {
            rows[r] += latticeVertex[best][r] + err[r] - 12 * Q;
            for (int h = 0; h < 6; ++h)
                rows[r] -= 3 * co[1 << r][h] * (c[h] > 0 ? v.delta[h].hi : v.delta[h].lo);
        }
        for (int i = 0; i < 4; ++i) {
            Wide z = 0;
            for (int j = 0; j < 4; ++j)
                z += Wide(6) * v.e[j] * col[j][v.e[j] >= 0][i];
            rows[15 + 2 * i] -= ceil_ratio(z, micro);
            rows[16 + 2 * i] += floor_ratio(z, micro);
        }
        for (int k = 0; k < 32; ++k) {
            auto f = facets[k];
            rows[25 + k] += f.dir[(best >> f.k) & 1];
        }
        batch.push_back(rows);
        w = proposal<K>(batch);
        if (*std::max_element(w.begin(), w.end()) == 0)
            break;
    }
    return false;
}
