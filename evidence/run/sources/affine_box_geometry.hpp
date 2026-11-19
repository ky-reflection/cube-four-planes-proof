#ifndef CUBE_AFFINE_BOX_GEOMETRY_HPP
#define CUBE_AFFINE_BOX_GEOMETRY_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <stdexcept>

namespace cube_affine_box {

using Int = std::int64_t;
using Wide = __int128_t;
constexpr Int scale = Int{1} << 40;
struct Range { Int lo, hi; };
constexpr Range one{scale,scale};

inline Int floor_ratio(Wide a, Int b) {
    if (b <= 0) throw std::logic_error("nonpositive denominator");
    return a < 0 ? -((-a+b-1)/b) : a/b;
}
inline Int ceil_ratio(Wide a, Int b) { return -floor_ratio(-a,b); }
inline Range operator+(Range a, Range b) { return {a.lo+b.lo,a.hi+b.hi}; }
inline Range operator-(Range a) { return {-a.hi,-a.lo}; }
inline Range operator-(Range a, Range b) { return a+-b; }
inline Range operator*(Range a, Range b) {
    const std::array<Wide,4> corners{Wide{a.lo}*b.lo,Wide{a.lo}*b.hi,Wide{a.hi}*b.lo,Wide{a.hi}*b.hi};
    const auto [lo,hi] = std::minmax_element(corners.begin(),corners.end());
    return {floor_ratio(*lo,scale),ceil_ratio(*hi,scale)};
}
inline Range operator/(Range a, Range b) {
    return {floor_ratio(Wide{scale}*a.lo,a.lo<0?b.lo:b.hi),ceil_ratio(Wide{scale}*a.hi,a.hi<0?b.hi:b.lo)};
}
inline Range operator/(Range a, Int b) { return {floor_ratio(a.lo,b),ceil_ratio(a.hi,b)}; }
inline Range nonnegative(Range a) { return {std::max(Int{0},a.lo),std::max(Int{0},a.hi)}; }
inline Range fraction(Range a) { return {std::clamp(a.lo,Int{0},scale),std::clamp(a.hi,Int{0},scale)}; }
inline Range intersect(Range a, Range b) { return {std::max(a.lo,b.lo),std::min(a.hi,b.hi)}; }

inline Int candidate_upper() {
    Int lo = scale/15, hi = scale/12+1;
    while (hi-lo > 1) {
        const Int mid = (lo+hi)/2;
        const Wide value = 48*mid-scale;
        if (value > 0 && value*value > Wide{8}*scale*scale) hi = mid;
        else lo = mid;
    }
    return hi;
}

constexpr std::array<int,6> faces{3,5,9,6,10,12};
using Areas = std::array<std::array<Range,16>,6>;
inline int face_shift(int face) {
    const int first=std::countr_zero(static_cast<unsigned>(face));
    return first^std::countr_zero(static_cast<unsigned>(face^(1<<first)));
}

// Each face has a default two-positive sign set. A corner triangle flips
// one sign, and an overlap of opposite types flips two signs.
template<class Number,class Triangle,class Overlap>
std::array<std::array<Number,16>,6> face_expressions(Number unit,Triangle triangle,Overlap overlap) {
    std::array<std::array<Number,16>,6> areas{};
    for (unsigned f=0; f<faces.size(); ++f) {
        const int face=faces[f];
        std::array<Number,4> triangles{}, overlap_sums{};
        Number total_triangles{}, total_overlap{};
        for (int i=0; i<4; ++i) {
            triangles[i]=triangle(face,i);
            total_triangles=total_triangles+triangles[i];
        }
        for (int i=0; i<4; ++i) if (face&(1<<i)) {
            for (int j=0; j<4; ++j) if (!(face&(1<<j))) {
                const auto o=overlap(face,i,j);
                areas[f][face^(1<<i)^(1<<j)]=o;
                overlap_sums[i]=overlap_sums[i]+o;
                overlap_sums[j]=overlap_sums[j]+o;
                total_overlap=total_overlap+o;
            }
        }
        areas[f][face]=fraction(unit-total_triangles+total_overlap);
        for (int i=0; i<4; ++i) areas[f][face^(1<<i)]=fraction(triangles[i]-overlap_sums[i]);
    }
    return areas;
}

template<class Triangle,class Overlap>
Areas face_areas(Triangle triangle,Overlap overlap) {
    return face_expressions(one,triangle,overlap);
}

// Expand the pyramid recurrence with every choice of the omitted labels.
// Height is 1-nu_F·p_i. Difference is nu_F·(p_j-p_i).
template<class Height, class Difference>
std::array<Int,16> cell_lower(const Areas& areas, Height height, Difference difference) {
    std::array<Int,16> lower{};
    for (int mask=1; mask<15; ++mask) {
        std::array<int,3> omitted{};
        int count=0;
        for (int i=0; i<4; ++i) if (!(mask&(1<<i))) omitted[count++]=i;
        for (int first=0; first<count; ++first) for (int second=0; second<count; ++second) {
            if (count>1 && first==second) continue;
            const int i=omitted[first], j=omitted[second];
            const int k=count==3?omitted[3-first-second]:-1;
            Range result{};
            for (unsigned f=0; f<faces.size(); ++f) {
                result=result+height(f,i)*areas[f][mask];
                if (count>1) result=result+difference(f,i,j)*areas[f][mask|(1<<i)];
                if (count==3) result=result+difference(f,k,j)*areas[f][mask|(1<<i)|(1<<j)];
            }
            lower[mask]=std::max(lower[mask],floor_ratio(result.lo,6));
        }
    }
    return lower;
}

} // namespace cube_affine_box

#endif
