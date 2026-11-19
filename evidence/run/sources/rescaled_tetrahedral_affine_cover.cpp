// Direct candidate-threshold test for coordinate-rescaled tetrahedral normals.
// The six parameters are u, z with C=(1,u,uz) and t_i with H_i=uz*t_i.
// Build: c++ -O3 -std=c++20 rescaled_tetrahedral_affine_cover.cpp -o /tmp/affine-cover
// Usage: affine-cover [depth] [budget] [uniform loss in millionths]

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "affine_box_geometry.hpp"

namespace {
using namespace cube_affine_box;
struct Box { std::array<Range,6> variables; int depth; };

class Model {
public:
    explicit Model(const Box& box) : u_(box.variables[0]), z_(box.variables[1]), v_(u_*z_) {
        for (int i=0; i<4; ++i) offsets_[i] = box.variables[i+2];
        total_ = nonnegative(offsets_[0]+offsets_[1]+offsets_[2]+offsets_[3]);
    }

    int exclude(Int target) const {
        if (!caps(target)) return 1;
        const auto areas = make_areas();
        std::array<Range,6> delta{};
        std::array<std::array<Range,4>,6> height{};
        for (unsigned f=0; f<faces.size(); ++f) {
            delta[f]=face_factor(faces[f])*total_;
            for (int i=0; i<4; ++i) height[f][i]=one-intersect(point_height(i,faces[f]),{-scale,scale});
        }
        const auto lower = cell_lower(areas,
            [&](unsigned f,int i) { return height[f][i]; },
            [&](unsigned f,int i,int j) {
                const int sign=((faces[f]>>i)&1)-((faces[f]>>j)&1);
                return sign==0?Range{}:sign>0?delta[f]:-delta[f];
            });
        for (int mask=1; mask<15; ++mask) if (lower[mask]>target) return 2;
        return 0;
    }

private:
    Range u_, z_, v_;
    std::array<Range,4> offsets_;
    Range total_;

    Range c(int kind) const { return kind==1?one:kind==2?u_:v_; }

    // v/(C_a+C_b), with the factor u cancelled for the two small entries.
    Range offset_factor(int a, int b) const {
        if (a!=1 && b!=1) return z_/(one+z_);
        return v_/(c(a)+c(b));
    }

    Range leg(int a, int b, Range offset) const {
        if (a!=1 && b!=1) {
            return fraction((a==2 ? one+z_*offset : z_*(one+offset))/(one+z_));
        }
        return fraction((c(a)+v_*offset)/(c(a)+c(b)));
    }

    bool caps(Int target) const {
        const Range denominator = (one+u_)*(one+v_);
        // Divide first by u+v, using v/(u+v)=z/(1+z).
        const Range linear = v_*(one+v_/(one+z_))/denominator;
        const Range cubic = (v_*v_*(z_/(one+z_)))/(denominator*Range{3*scale,3*scale});
        for (Range offset : offsets_) {
            const auto deficit = linear*offset-cubic*offset*offset*offset;
            if (deficit.lo > 8*target-scale/2 || deficit.hi < scale/2-7*target) return false;
        }
        return true;
    }

    Range pair_sum(int mask) const {
        Range result{};
        for (int i=0; i<4; ++i) if (mask&(1<<i)) result=result+offsets_[i];
        return result;
    }

    Range face_factor(int mask) const {
        const int a=std::countr_zero(static_cast<unsigned>(mask));
        const int b=std::countr_zero(static_cast<unsigned>(mask^(1<<a)));
        const int shift=a^b;
        if (shift==1) return z_/(one+z_);
        return v_/(one+c(shift==2?3:2));
    }

    Range point_height(int omitted, int face) const {
        return (face&(1<<omitted)) ? -face_factor(face)*pair_sum(15^face) : face_factor(face)*pair_sum(face);
    }

    Range overlap(int a, int b, int e, Range t) const {
        // The center difference between the two edge fractions is (A-B)/(A+B).
        Range centered;
        if (a!=1 && b!=1) {
            centered=(one-z_)/(one+z_);
            if (a==3) centered=-centered;
        } else {
            centered=(one-c(a==1?b:a))/(one+c(a==1?b:a));
            if (a!=1) centered=-centered;
        }
        const auto delta=nonnegative(centered+offset_factor(a,b)*t);
        if (delta.hi==0) return {};
        const auto delta_sq=delta*delta;
        // The ratio (A+B)/(A+E) is at most four on positive support:
        // when B=1 and A,E are the two small entries, delta>0 gives
        // A+2v>1, hence A+E>1/2. In other cases the ratio is at most two.
        Range ratio{0,4*scale};
        const auto denominator=c(a)+c(e);
        if (denominator.lo>0) ratio=intersect(ratio,(c(a)+c(b))/denominator);
        if (ratio.lo>ratio.hi) return {}; // No point on the positive branch.
        return fraction(delta_sq*ratio/4);
    }

    Areas make_areas() const {
        return face_areas(
            [&](int face,int i) {
                const int shift=face_shift(face);
                const Range h=face&(1<<i)?offsets_[i]:-offsets_[i];
                Range t=one;
                for (int k=1; k<=3; ++k) if (k!=shift) t=t*leg(shift,k,h);
                return t/2;
            },
            [&](int face,int i,int j) {
                const int shift=face_shift(face);
                return overlap(shift,i^j^shift,i^j,offsets_[i]-offsets_[j]);
            });
    }
};

bool restrict_offsets(Box& box) {
    auto* t=box.variables.data()+2;
    for (int i=1; i<4; ++i) {
        t[0].hi=std::min(t[0].hi,t[i].hi);
        t[i].lo=std::max(t[i].lo,t[0].lo);
    }
    const Int upper=t[0].hi+t[1].hi+t[2].hi+t[3].hi;
    for (int i=0; i<4; ++i) {
        t[i].lo=std::max(t[i].lo,t[i].hi-upper);
        if (t[i].lo>t[i].hi) return false;
    }
    return true;
}

template<class T> T parse(std::string_view arg) {
    T value{};
    const auto [end,error]=std::from_chars(arg.data(),arg.data()+arg.size(),value);
    if (error!=std::errc{} || end!=arg.data()+arg.size()) throw std::invalid_argument("invalid integer");
    return value;
}
} // namespace

int main(int argc, char** argv) {
    const int depth=argc>1?parse<int>(argv[1]):48;
    const auto budget=argc>2?parse<std::uint64_t>(argv[2]):1000000;
    const auto loss=argc>3?parse<int>(argv[3]):0;
    if (depth<0 || depth>120 || loss<0 || loss>3000) return 2;
    const Int target=candidate_upper()+ceil_ratio(Wide{scale}*loss,1000000);
    Box initial{};
    initial.variables={Range{0,scale},Range{0,scale},Range{-scale,scale},Range{-scale,scale},Range{-scale,scale},Range{-scale,scale}};
    std::vector<Box> stack{initial};
    std::array<std::uint64_t,3> excluded{};
    std::uint64_t visited=0, unresolved=0;
    while (!stack.empty() && visited<budget) {
        Box box=stack.back();
        stack.pop_back();
        ++visited;
        const int verdict=restrict_offsets(box)?Model(box).exclude(target):1;
        if (verdict) { ++excluded[verdict]; continue; }
        if (box.depth==depth) {
            if (unresolved==0) {
                std::cerr << "first unresolved box";
                for (auto [a,b]:box.variables) std::cerr << " [" << double(a)/scale << ',' << double(b)/scale << ']';
                std::cerr << '\n';
            }
            ++unresolved;
            continue;
        }
        int axis=box.depth%6;
        for (int tries=0; tries<6 && box.variables[axis].hi-box.variables[axis].lo<2; ++tries) axis=(axis+1)%6;
        const auto [lo,hi]=box.variables[axis];
        if (hi-lo<2) { ++unresolved; continue; }
        const Int midpoint=(lo+hi)/2;
        Box right=box;
        ++box.depth;
        ++right.depth;
        box.variables[axis].hi=midpoint;
        right.variables[axis].lo=midpoint;
        stack.push_back(right);
        stack.push_back(box);
    }
    std::cout << "visited=" << visited << " feasibility=" << excluded[1] << " volume=" << excluded[2]
              << " unresolved=" << unresolved << " pending=" << stack.size() << '\n';
    return unresolved || !stack.empty();
}
