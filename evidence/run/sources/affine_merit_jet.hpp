// Common third-order interval arithmetic. D is the number of local coordinates.
constexpr int H = D * (D + 1) / 2;
constexpr int T = H * (D + 2) / 3;
constexpr int hh(int i, int j) { return i * (i + 1) / 2 + j; }
using V3 = array<struct Jet, 3>;
Range ratio(Int a, Int b) { return {floor_ratio(Wide(a) * scale, b), ceil_ratio(Wide(a) * scale, b)}; }
constexpr Range exact(Int a) { return {a * scale, a * scale}; }
Range reciprocal(Range x) {
    if (x.hi < 0)
        return -(one / (-x));
    if (x.lo <= 0)
        throw std::logic_error("singular");
    return one / x;
}
struct Jet {
    Range v{};
    array<Range, D> d{};
    array<Range, H> h{};
    array<Range, T> t{};
    Jet() = default;
    explicit Jet(Range a) : v(a) {}
    static Jet variable(Range v, int i) {
        Jet x(v);
        x.d[i] = one;
        return x;
    }
};
Jet operator+(const Jet &a, const Jet &b) {
    Jet r(a.v + b.v);
    for (int i = 0; i < D; ++i)
        r.d[i] = a.d[i] + b.d[i];
    for (int i = 0; i < H; ++i)
        r.h[i] = a.h[i] + b.h[i];
    for (int i = 0; i < T; ++i)
        r.t[i] = a.t[i] + b.t[i];
    return r;
}
Jet operator-(const Jet &a) {
    Jet r(-a.v);
    for (int i = 0; i < D; ++i)
        r.d[i] = -a.d[i];
    for (int i = 0; i < H; ++i)
        r.h[i] = -a.h[i];
    for (int i = 0; i < T; ++i)
        r.t[i] = -a.t[i];
    return r;
}
Jet operator-(const Jet &a, const Jet &b) { return a + -b; }
Jet operator*(const Jet &a, const Jet &b) {
    Jet r(a.v * b.v);
    for (int i = 0; i < D; ++i)
        r.d[i] = a.d[i] * b.v + a.v * b.d[i];
    int k = 0;
    for (int i = 0; i < D; ++i)
        for (int j = 0; j <= i; ++j, ++k) {
            r.h[k] = a.h[k] * b.v + b.h[k] * a.v + a.d[i] * b.d[j] + a.d[j] * b.d[i];
        }
    int index = 0;
    for (int i = 0; i < D; ++i)
        for (int j = 0; j <= i; ++j)
            for (int k = 0; k <= j; ++k, ++index)
                r.t[index] = a.t[index] * b.v + b.t[index] * a.v + a.h[hh(i, j)] * b.d[k] +
                             a.h[hh(i, k)] * b.d[j] + a.h[hh(j, k)] * b.d[i] + b.h[hh(i, j)] * a.d[k] +
                             b.h[hh(i, k)] * a.d[j] + b.h[hh(j, k)] * a.d[i];
    return r;
}
Jet operator*(const Jet &a, Range b) {
    Jet r(a.v * b);
    for (int i = 0; i < D; ++i)
        r.d[i] = a.d[i] * b;
    for (int i = 0; i < H; ++i)
        r.h[i] = a.h[i] * b;
    for (int i = 0; i < T; ++i)
        r.t[i] = a.t[i] * b;
    return r;
}
Jet operator/(const Jet &a, int b) { return a * ratio(1, b); }
Jet unary(const Jet &a, Range v, Range dp, Range d2, Range d3) {
    Jet r(v);
    int k = 0;
    for (int i = 0; i < D; ++i) {
        r.d[i] = dp * a.d[i];
        for (int j = 0; j <= i; ++j, ++k)
            r.h[k] = d2 * a.d[i] * a.d[j] + dp * a.h[k];
    }
    int index = 0;
    for (int i = 0; i < D; ++i)
        for (int j = 0; j <= i; ++j)
            for (int k = 0; k <= j; ++k, ++index)
                r.t[index] =
                    d3 * a.d[i] * a.d[j] * a.d[k] + dp * a.t[index] +
                    d2 * (a.h[hh(i, j)] * a.d[k] + a.h[hh(j, k)] * a.d[i] + a.h[hh(i, k)] * a.d[j]);
    return r;
}
Jet inverse(const Jet &a) {
    Range b = reciprocal(a.v), b2 = b * b;
    return unary(a, b, -b2, exact(2) * b2 * b, exact(-6) * b2 * b2);
}
Jet operator/(const Jet &a, const Jet &b) { return a * inverse(b); }
Jet sq(const Jet &a) {
    Range v = a.v * a.v;
    v.lo = std::max<Int>(v.lo, 0);
    return unary(a, v, exact(2) * a.v, exact(2), exact(0));
}
Jet abs3(const Jet &a) {
    Range abs{std::max({Int{0}, a.v.lo, -a.v.hi}), std::max(-a.v.lo, a.v.hi)};
    return unary(a, abs * abs * abs, exact(3) * a.v * abs, exact(6) * abs,
                 Range{a.v.lo > 0 ? 6 * scale : -6 * scale, a.v.hi < 0 ? -6 * scale : 6 * scale});
}
V3 cross(const V3 &a, const V3 &b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Jet dot(const V3 &a, const V3 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
