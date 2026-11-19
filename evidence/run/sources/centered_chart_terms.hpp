// Centered triangle controls, edge overlaps and their direction bounds.
#ifndef CUBE_CENTERED_CHART_TERMS_HPP
#define CUBE_CENTERED_CHART_TERMS_HPP
Bound edge_overlap(const Matrix& c, int i, int j) {
    int k = -1, l = -1;
    for (int a = 0; a < 4; ++a) if (a != i && a != j) {
        if (k < 0) k = a;
        else l = a;
    }
    const auto p = share(c[k][i], c[k][j]);
    const auto q = share(c[l][i], c[l][j]);
    Bound delta{p.lo-q.hi, p.hi-q.lo};
    Bound result{0, 0};
    for (int sign = 0; sign < 2; ++sign) {
        if (delta.hi > 0) {
            const auto first = linked_ratio(c[k][i], c[k][j], c[k][l]);
            const auto second = linked_ratio(c[l][j], c[l][i], c[l][k]);
            const Int den_low = first.lo + second.lo;
            const Wide high_square = Wide{delta.hi} * delta.hi;
            // The transverse triangle legs are at most one, so the
            // denominator on this branch is also at least 1+delta.
            const Int cap = ceil_ratio(high_square, mass + delta.hi);
            const Int upper = den_low ? std::min(cap, ceil_ratio(high_square, den_low)) : cap;
            result.hi = std::max(result.hi, upper);
            if (delta.lo > 0 && first.hi != infinity && second.hi != infinity) {
                result.lo = std::max(result.lo, floor_ratio(Wide{delta.lo} * delta.lo, first.hi + second.hi));
            }
        }
        std::swap(k, l);
        delta = {-delta.hi, -delta.lo};
    }
    return result;
}

int pattern_id(Pattern p) {
    return (p[0]+1)*9 + (p[1]+1)*3 + p[2]+1;
}

struct Combination {
    Int constant;
    unsigned mask;
    int cells;
    std::array<int, 4> row_patterns;
    std::array<int, 6> overlap_coefficients;
};

struct AffineTriangle {
    Int center;
    std::array<Int, 2> slope;
    Bound error;
};

struct RowControls {
    std::array<Int, 9> denominator;
    std::array<std::array<Int, 9>, 3> numerator;
};

RowControls row_controls(Coordinate x, Coordinate y) {
    const std::array<Int,3> X{2*x.lo, x.lo+x.hi, 2*x.hi};
    const std::array<Int,3> Y{2*y.lo, y.lo+y.hi, 2*y.hi};
    const std::array<Int,3> XX{Int{x.lo}*x.lo, Int{x.lo}*x.hi, Int{x.hi}*x.hi};
    const std::array<Int,3> YY{Int{y.lo}*y.lo, Int{y.lo}*y.hi, Int{y.hi}*y.hi};
    RowControls result{};
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        const int index = 3*i+j;
        const Int n0 = (X[i]+Y[j])*grid*grid;
        const Int n1 = XX[i]*(2*grid+Y[j]);
        const Int n2 = YY[j]*(2*grid+X[i]);
        result.numerator[0][index] = n0;
        result.numerator[1][index] = n1;
        result.numerator[2][index] = n2;
        result.denominator[index] = n0+n1+n2+X[i]*Y[j]*grid;
    }
    return result;
}

std::array<AffineTriangle, 3> affine_row(Coordinate x, Coordinate y) {
    // Floating-point midpoint Taylor coefficients are only proposals.
    // Their uniform errors are certified below by integer Bernstein bounds.
    const std::array<double,3> midpoint{1., (x.lo+x.hi)/(2.*grid), (y.lo+y.hi)/(2.*grid)};
    const std::array<double,2> halfwidth{(x.hi-x.lo)/(2.*grid), (y.hi-y.lo)/(2.*grid)};
    const auto controls = row_controls(x,y);
    std::array<AffineTriangle, 3> result{};
    for (int term = 0; term < 3; ++term) {
        const double a = midpoint[term];
        const double b = midpoint[(term+1)%3];
        const double c = midpoint[(term+2)%3];
        const double value = a*a/((a+b)*(a+c));
        auto& atom = result[term];
        atom.center = std::llround(mass*value);
        for (int k = 1; k < 3; ++k) {
            const double derivative = term == k ? (2./a-1./(a+b)-1./(a+c)) : -1./(a+midpoint[k]);
            atom.slope[k-1] = std::llround(mass*halfwidth[k-1]*value*derivative);
        }
        atom.error = {infinity,-infinity};
        for (int i = 0; i <= 3; ++i) for (int j = 0; j <= 3; ++j) {
            Wide numerator = 0;
            Int denominator = 0;
            for (int u = 0; u <= 1; ++u) for (int v = 0; v <= 1; ++v) {
                const int weight = (u ? i : 3-i)*(v ? j : 3-j);
                if (!weight) continue;
                const int index = 3*(i-u)+(j-v);
                const Int local_den = controls.denominator[index];
                const Int local_linear = atom.center+(2*u-1)*atom.slope[0]+(2*v-1)*atom.slope[1];
                numerator += weight*(Wide{mass}*controls.numerator[term][index]-Wide{local_den}*local_linear);
                denominator += weight*local_den;
            }
            if (!denominator) continue;
            atom.error.lo = std::min(atom.error.lo, floor_ratio(numerator,denominator));
            atom.error.hi = std::max(atom.error.hi, ceil_ratio(numerator,denominator));
        }
        if (atom.error.lo > atom.error.hi) throw std::logic_error("empty affine row control net");
    }
    return result;
}

using Weights = std::array<Int, 7>;
using LinearSlopes = std::array<Int, 8>;

struct AffineOverlap {
    Int center;
    LinearSlopes slope;
    Bound error;
};

Bound add_interval(Bound a, Bound b) { return {a.lo+b.lo,a.hi+b.hi}; }
Bound negate_interval(Bound a) { return {-a.hi,-a.lo}; }

Bound multiply_interval(Bound a, Bound b) {
    if ((a.lo == 0 && a.hi == 0) || (b.lo == 0 && b.hi == 0)) return {};
    const std::array<Wide,4> corners{Wide{a.lo}*b.lo,Wide{a.lo}*b.hi,Wide{a.hi}*b.lo,Wide{a.hi}*b.hi};
    const auto [lo,hi] = std::minmax_element(corners.begin(),corners.end());
    return {floor_ratio(*lo,mass),ceil_ratio(*hi,mass)};
}

Bound divide_interval(Bound a, Bound b) {
    if (a.lo == 0 && a.hi == 0) return {};
    return {floor_ratio(Wide{mass}*a.lo,a.lo < 0 ? b.lo : b.hi),
            ceil_ratio(Wide{mass}*a.hi,a.hi < 0 ? b.hi : b.lo)};
}

struct FirstDerivative {
    Bound value;
    std::array<Bound,8> slope;
};

FirstDerivative sum_jet(const FirstDerivative& a, const FirstDerivative& b) {
    FirstDerivative result{add_interval(a.value,b.value),{}};
    for (int i = 0; i < 8; ++i) result.slope[i] = add_interval(a.slope[i],b.slope[i]);
    return result;
}

FirstDerivative negate_jet(const FirstDerivative& a) {
    FirstDerivative result{negate_interval(a.value),{}};
    for (int i = 0; i < 8; ++i) result.slope[i] = negate_interval(a.slope[i]);
    return result;
}

FirstDerivative ratio_jet(const FirstDerivative& a, const FirstDerivative& b) {
    FirstDerivative result{divide_interval(a.value,b.value),{}};
    for (int i = 0; i < 8; ++i) {
        const Bound product = multiply_interval(result.value,b.slope[i]);
        result.slope[i] = divide_interval(add_interval(a.slope[i],negate_interval(product)),b.value);
    }
    return result;
}

FirstDerivative overlap_branch(const FirstDerivative& delta, FirstDerivative denominator) {
    if (delta.value.hi <= 0) return {};
    const Bound positive{std::max(Int{0},delta.value.lo),delta.value.hi};
    // This branch is zero when delta <= 0, with matching first derivative
    // at zero. On its positive support the denominator is at least one.
    denominator.value.lo = std::max(denominator.value.lo,mass);
    denominator.value.hi = std::max(denominator.value.hi,mass);
    FirstDerivative square{multiply_interval(positive,positive),{}};
    for (int i = 0; i < 8; ++i) {
        square.slope[i] = multiply_interval({2*positive.lo,2*positive.hi},delta.slope[i]);
    }
    auto result = ratio_jet(square,denominator);
    if (delta.value.lo <= 0) {
        for (auto& component : result.slope) {
            component = {std::min(Int{0},component.lo),std::max(Int{0},component.hi)};
        }
    }
    return result;
}

FirstDerivative overlap_jet(const std::array<std::array<FirstDerivative,4>,4>& c, int i, int j) {
    int k = -1, l = -1;
    for (int a = 0; a < 4; ++a) if (a != i && a != j) {
        if (k < 0) k = a;
        else l = a;
    }
    const auto ab = sum_jet(c[k][i],c[k][j]);
    const auto cd = sum_jet(c[l][i],c[l][j]);
    const auto delta = sum_jet(ratio_jet(c[k][i],ab),negate_jet(ratio_jet(c[l][i],cd)));
    FirstDerivative forward{}, backward{};
    if (delta.value.hi > 0) {
        const auto denominator = sum_jet(ratio_jet(sum_jet(c[k][i],c[k][l]),ab),
                                         ratio_jet(sum_jet(c[l][j],c[l][k]),cd));
        forward = overlap_branch(delta,denominator);
    }
    if (delta.value.lo < 0) {
        const auto denominator = sum_jet(ratio_jet(sum_jet(c[l][i],c[l][k]),cd),
                                         ratio_jet(sum_jet(c[k][j],c[k][l]),ab));
        backward = overlap_branch(negate_jet(delta),denominator);
    }
    return sum_jet(forward,backward);
}

#endif
