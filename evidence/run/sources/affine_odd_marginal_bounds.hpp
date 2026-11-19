// Bound each signed marginal by linear combinations of the seven opposite-cell differences.
bool MomentDomain::inverse_caps(Caps &cap, const I (&H)[4], const std::array<Int, 128> &slack) const {
    const auto P = section_coefficients(*odd, ratios, H);
    ProbabilityColumns columns(*this, cap, P);
    double L[7][4]{};
    for (int j = 0; j < 4; ++j) {
        double p[4]{}, twice = -1;
        for (int i = 0; i < 4; ++i)
            if (i != j) {
                p[i] = double(P[i][j].l + P[i][j].h) / (2 * mass);
                twice += p[i];
            }
        for (int i = 0; i < 4; ++i)
            L[i][j] = 2 * p[i] - twice;
        for (int r = 4; r < 7; ++r)
            for (int i = 0; i < 4; ++i)
                L[r][j] += ((reps[r] >> i & 1) ? -1 : 1) * ((i == j) + p[i]);
    }
    struct Form {
        std::array<int, 7> y{};
        double price = 1e100;
    };
    Form best[4][3];
    auto cost = [&](const std::array<int, 7> &y) {
        std::array<std::pair<int, int>, 7> vals;
        for (int r = 0; r < 7; ++r)
            vals[r] = {std::abs(y[r]), r};
        std::sort(vals.begin(), vals.end());
        Int value = 0;
        int mask = 127, prior = 0;
        for (auto [n, r] : vals) {
            value += Int(n - prior) * slack[mask];
            prior = n;
            mask ^= 1 << r;
        }
        return value;
    };
    for (int mask = 1; mask < 128; ++mask)
        if (std::popcount(unsigned(mask)) == 4) {
            int row[4], nr = 0;
            for (int r = 0; r < 7; ++r)
                if (mask >> r & 1)
                    row[nr++] = r;
            double a[4][8]{};
            for (int i = 0; i < 4; ++i) {
                a[i][4 + i] = 1;
                for (int j = 0; j < 4; ++j)
                    a[i][j] = L[row[j]][i];
            }
            bool ok = true;
            for (int j = 0; j < 4; ++j) {
                int k = j;
                for (int i = j + 1; i < 4; ++i)
                    if (std::abs(a[i][j]) > std::abs(a[k][j]))
                        k = i;
                if (std::abs(a[k][j]) < .02) {
                    ok = false;
                    break;
                }
                for (int l = 0; l < 8; ++l)
                    std::swap(a[k][l], a[j][l]);
                double t = a[j][j];
                for (double &v : a[j])
                    v /= t;
                for (int i = 0; i < 4; ++i)
                    if (i != j) {
                        double t = a[i][j];
                        for (int l = 0; l < 8; ++l)
                            a[i][l] -= t * a[j][l];
                    }
            }
            if (!ok)
                continue;
            for (int j = 0; j < 4; ++j) {
                double scale = 0;
                for (int i = 0; i < 4; ++i)
                    scale = std::max(scale, std::abs(a[i][4 + j]));
                if (!std::isfinite(scale) || scale <= 0)
                    continue;
                Form v;
                for (int i = 0; i < 4; ++i)
                    v.y[row[i]] = std::lround(1024 * a[i][4 + j] / scale);
                v.price = cost(v.y) * scale;
                for (auto &f : best[j])
                    if (v.price < f.price)
                        std::swap(f, v);
            }
        }
    for (int target = 0; target < 4; ++target)
        for (auto &f : best[target])
            if (f.price < 1e99) {
                auto &y = f.y;
                int S = y[0] + y[1] + y[2] + y[3], base[4]{};
                std::array<int, 4> b{};
                for (int i = 0; i < 4; ++i) {
                    int t = 0;
                    for (int r = 4; r < 7; ++r)
                        t += ((reps[r] >> i & 1) ? -1 : 1) * y[r];
                    base[i] = S + t;
                    b[i] = 2 * y[i] - S + t;
                }
                auto vals = columns.supports(b);
                I coeff[4];
                for (int j = 0; j < 4; ++j) {
                    I q{base[j] * mass, base[j] * mass};
                    for (int i = 0; i < 4; ++i) {
                        q.l += Int(b[i]) * vals[j][1][i];
                        q.h += Int(b[i]) * vals[j][0][i];
                    }
                    coeff[j] = q;
                }
                if (absmin(coeff[target]) < mass * 4)
                    continue;
                Int room = 2 * cost(y);
                I rhs{-room, room};
                for (int j = 0; j < 4; ++j)
                    if (j != target)
                        rhs = sub(rhs, mul(coeff[j], cap[j]));
                cap[target] = intersection(cap[target], signed_div(rhs, coeff[target]));
                if (cap[target].h < cap[target].l)
                    return false;
            }
    return true;
}
