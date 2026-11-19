// Joint density of two cuts and the boundary of their projected sections.
struct PairTransfer {
    bool valid = false;
    R A{}, B{}, C{};
    Int rate[2]{};
    R at(R e, R f) const { return A * square(e) + times(B * e * f, 2) + C * square(f); }
};
R signed_ratio(R a, R b) { return b.hi < 0 ? (-a) / (-b) : a / b; }
R surface_density(R u, R v, R z) {
    Int lo = minabs(z), hi = maxabs(z), amax = maxabs(u), bmax = maxabs(v), amin = minabs(u),
        bmin = minabs(v);
    if (amax + bmax < 2 * lo)
        return {};
    Int upper = lo > 0 ? ceil_ratio(Wide(U) * U, 2 * lo) : 8 * U;
    if (std::max(amin, bmin) > 0)
        upper = std::min(upper, ceil_ratio(Wide(U) * U, std::max(amin, bmin)));
    R prod = pt(amin) * pt(bmin);
    if (prod.lo > U / 16384)
        upper = std::min(upper, ceil_ratio(Wide(U) * (amax + bmax - 2 * lo), 2 * prod.lo));
    Int lower = 0;
    if (amin + bmin > 2 * hi && amax && bmax) {
        Int denom = (pt(amax) * pt(bmax)).hi;
        lower = std::min(floor_ratio(Wide(U) * U, std::max(amax, bmax)),
                         floor_ratio(Wide(U) * (amin + bmin - 2 * hi), 2 * denom));
    }
    return {lower, upper};
}
R edge_term(int i, int j, int sg, const std::array<V, 4> &normals, R h, R other) {
    int k = type.axis[i], x = (k + 1) % 3, y = (k + 2) % 3;
    R u = normals[i][x], v = normals[i][y], z = h + pt(sg * U / 2);
    auto d = surface_density(u, v, z);
    if (!d.hi)
        return {};
    R X{-U / 2, U / 2}, Y = X;
    if (minabs(u) > U / 64)
        X = intersect(X, signed_ratio(z - v * Y, u));
    if (minabs(v) > U / 64)
        Y = intersect(Y, signed_ratio(z - u * X, v));
    if (X.lo > X.hi || Y.lo > Y.hi)
        return {};
    R a = normals[j][x] * X + normals[j][y] * Y - times(normals[j][k], sg * type.orient[i]) / 2 - other;
    if (a.lo > 0)
        return d;
    if (a.hi < 0)
        return -d;
    return {-d.hi, d.hi};
}
PairTransfer transfer_coefficients(int i, int j, const std::array<V, 4> &n, const Caps &H, const Caps &rho) {
    PairTransfer out;
    V w{};
    Int low = 0, high = 0;
    int k = 0;
    for (int t = 0; t < 3; ++t) {
        int x = (t + 1) % 3, y = (t + 2) % 3;
        w[t] = n[i][x] * n[j][y] - n[i][y] * n[j][x];
        if (minabs(w[t]) > low) {
            low = minabs(w[t]);
            k = t;
        }
        high = std::max(high, maxabs(w[t]));
    }
    if (low < U / 8)
        return out;
    int x = (k + 1) % 3, y = (k + 2) % 3;
    R zero_coordinate[2] = {signed_ratio(H[i] * n[j][y] - H[j] * n[i][y], w[k]),
                            signed_ratio(n[i][x] * H[j] - n[j][x] * H[i], w[k])};
    Int full = floor_ratio(Wide(U) * U, high), end = floor_ratio(Wide(U) * U, maxabs(w[k]));
    for (int a = 0; a < 2; ++a) {
        int axis = a ? y : x;
        Int space = U - 2 * maxabs(zero_coordinate[a]);
        if (space <= 0) {
            end = 0;
            break;
        }
        if (maxabs(w[axis]) > 0)
            end = std::min(end, floor_ratio(Wide(space) * U, maxabs(w[axis])));
    }
    R ell{(2 * full + end) / 3, ceil_ratio(Wide(U) * U, low)};
    R diag[2];
    for (int t = 0; t < 2; ++t) {
        int r = t ? j : i, l = t ? i : j, axis = type.axis[r];
        R bend{};
        if (!type.face[r])
            for (int s = 0; s < 4; ++s) {
                R seg{U * s / 4, U * (s + 1) / 4};
                auto h = H[r] * seg, z = H[l] * seg;
                // Integral of (1-s) on this quarter is (7-2*s)/32.
                bend =
                    bend + times(edge_term(r, l, 1, n, h, z) - edge_term(r, l, -1, n, h, z), 7 - 2 * s) / 32;
            }
        diag[t] = -(times(n[l][axis], type.orient[r]) * ell + bend) / square(rho[r]);
    }
    out.valid = true;
    out.A = diag[0];
    out.C = diag[1];
    out.B = ell / (rho[i] * rho[j]);
    for (int s = 0; s < 2; ++s) {
        R a = s ? out.A : -out.A, c = s ? out.C : -out.C;
        R d = times(out.B, 2) - a - c;
        out.rate[s] = 3000000;
        if (d.lo > U / 32768) {
            R h = (square(out.B) - a * c) / d;
            Int K = std::max({Int(0), a.hi, c.hi, h.hi});
            out.rate[s] = ceil_ratio(Wide(1000000) * K, U);
        }
    }
    return out;
}
