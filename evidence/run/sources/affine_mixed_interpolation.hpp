// Face-pyramid interpolation for a pair of normalized physical directions.
struct PairLattice {
    R corner[16]{}, curvature[2][2]{};
    bool smooth = false;
};
std::unordered_map<PairKey, PairLattice, KeyHash> lattice_cache[6];
std::unordered_map<PairKey, R, KeyHash> polygon_cache[6];
R oriented_volume(const std::array<seven_integer::Vector, 2> &n) {
    Wide lo = 0, hi = 0;
    for (int k = 0; k < 3; ++k)
        for (int face : {-1, 1}) {
            cube_face_integral::L lines[2];
            for (int w = 0; w < 2; ++w) {
                int t = w ? -1 : 1;
                lines[w] = {2 * t * n[w][(k + 1) % 3], 2 * t * n[w][(k + 2) % 3], -t * face * n[w][k]};
            }
            auto a = cube_face_integral::projected_area(lines, 2);
            lo += a.first;
            hi += a.second;
        }
    constexpr Int gc = cube_face_integral::coord_grid;
    return {floor_ratio(lo * U, 6 * gc * gc), ceil_ratio(hi * U, 6 * gc * gc)};
}
PairLattice face_lattice(int h, const Box &box, const std::array<V, 4> &n) {
    auto [i, j] = pairid[h];
    PairKey key{};
    int at = 0;
    for (int r : {i, j})
        for (int k = 0; k < 2; ++k) {
            auto b = box[2 * r + k];
            key.vals[at++] = b.lo;
            key.vals[at++] = b.hi;
        }
    auto &cache = lattice_cache[h];
    auto iter = cache.find(key);
    if (iter != cache.end())
        return iter->second;
    PairLattice out{};
    for (int m = 0; m < 16; ++m) {
        std::array<seven_integer::Vector, 2> normals{};
        PairKey vert{};
        int pos = 0;
        for (int w = 0; w < 2; ++w) {
            int r = w ? j : i, k = type.axis[r];
            Int v[2];
            for (int z = 0; z < 2; ++z) {
                R value = box[2 * r + z];
                Int c = m >> pos & 1 ? value.hi : value.lo;
                vert.vals[pos++] = c;
                assert(c % (U / 8192) == 0);
                v[z] = c / (U / 8192);
            }
            normals[w][k] = 16384 * type.orient[r];
            normals[w][(k + 1) % 3] = v[0] - v[1];
            normals[w][(k + 2) % 3] = v[0] + v[1];
        }
        auto &c = polygon_cache[h];
        auto it = c.find(vert);
        if (it != c.end())
            out.corner[m] = it->second;
        else {
            out.corner[m] = oriented_volume(normals);
            if (c.size() > 200000)
                c.clear();
            c.emplace(vert, out.corner[m]);
        }
    }
    V t;
    R length{};
    for (int k = 0; k < 3; ++k) {
        int x = (k + 1) % 3, y = (k + 2) % 3;
        t[k] = n[i][x] * n[j][y] - n[j][x] * n[i][y];
        length = {std::max(length.lo, minabs(t[k])), std::max(length.hi, maxabs(t[k]))};
    }
    if (length.lo > U / 24) {
        out.smooth = true;
        for (int w = 0; w < 2; ++w) {
            int r = w ? j : i, l = w ? i : j, k = type.axis[r], x = (k + 1) % 3, y = (k + 2) % 3;
            R u = n[r][x], v = n[r][y];
            Int aa = maxabs(u), bb = maxabs(v);
            Int omit = 0;
            if (!type.face[r] && aa + bb > U) {
                // A symmetric unimodal density is at most one at distance 1/2.
                Int pdf = U, den = (pt(minabs(u)) * pt(minabs(v))).lo;
                if (den > U / 512)
                    pdf = std::min(pdf, ceil_ratio(Wide(aa + bb - U) * U, 2 * den));
                omit = (pdf + 1) / 2;
            }
            for (int z = 0; z < 2; ++z) {
                auto prod = (times(t[x], z ? -1 : 1) + t[y]) / 2;
                auto line = square(prod) / (times(square(length) * length, 12));
                out.curvature[w][z] = times(times(n[l][k], type.orient[r]) * line, 2) + R{-omit, omit};
            }
        }
    }
    if (cache.size() > 70000)
        cache.clear();
    cache.emplace(key, out);
    return out;
}
