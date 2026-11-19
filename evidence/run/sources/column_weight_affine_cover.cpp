// Candidate-threshold cover for C_ij=c_j and independent offsets.
// Build: c++ -O3 -std=c++20 column_weight_affine_cover.cpp -o /tmp/column-cover
// Usage: column-cover [depth] [budget] [affine-start-depth]
// The local handoff is justified by column_weight_local_certificate.py.

#include <charconv>
#include <iostream>
#include <string_view>
#include <vector>

#include "affine_box_geometry.hpp"
#include "affine_box_linear.hpp"

namespace {
using namespace cube_affine_box;
struct Box { std::array<Range,7> variables; int depth; };

template<class Number>
class Geometry {
public:
    explicit Geometry(const std::array<Number,7>& parameters,Number unit) : one_(unit) {
        const auto s=parameters[0], u=parameters[1], v=parameters[2];
        adjacent_={v,u,s};
        c_={s*u*v,s*u,s,one_};
        sum_=c_[0]+c_[1]+c_[2]+one_;
        for (int i=0; i<4; ++i) offsets_[i]=parameters[3+i];
        signed_total_=c_[1]*(offsets_[0]+offsets_[1])+c_[2]*offsets_[2]+offsets_[3];
    }
    Number signed_total() const { return signed_total_; }

    Number triangle(int face,int cut) const {
        const int corner=cut^face_shift(face);
        const Number t=face&(1<<cut)?offsets_[cut]:-offsets_[cut];
        Number area=one_;
        for (int other=0; other<4; ++other) if (other!=cut && other!=corner) {
            area=area*leg(cut,corner,other,t);
        }
        return area/2;
    }

    Number overlap(int face,int i,int j) const {
        const int a=i^face_shift(face), b=j^face_shift(face);
        const int min_i=min_column(i), min_j=min_column(j);
        const int big=std::max(min_i,min_j);
        const auto difference=ratio(min_i,big)*offsets_[i]-ratio(min_j,big)*offsets_[j];
        return fraction((c_[big]/sum_)*part(big,a,b)*positive_square(difference)/2);
    }

    Number point_height(int omitted,int face) const {
        // Cancel H_omitted before bounding q - (bar H/c_omitted) v_omitted.
        const int sign_j=face&(1<<omitted)?1:-1;
        Number value{};
        for (int i=0; i<4; ++i) if (i!=omitted) {
            const int sign_i=face&(1<<i)?1:-1;
            const Number first=sign_i>0?c_[omitted]:-c_[omitted];
            const Number second=sign_j>0?c_[i]:-c_[i];
            value=value+ratio(min_column(i),omitted)*(first-second)*offsets_[i];
        }
        return value/sum_;
    }

    Number shift(int omitted,int face,Number total) const {
        const Number change=ratio(0,omitted)*total;
        return face&(1<<omitted)?change:-change;
    }

    Number weighted_total() const { return signed_total_/sum_; }

    Number cap_deficit(int i) const {
        std::array<int,3> labels{};
        int used=0;
        for (int j=0; j<4; ++j) if (i!=j) labels[used++]=j;
        const auto [a,b,c]=labels;
        const Number u=ratio(a,b),v=ratio(a,c);
        const Number t=offsets_[i];
        return part(a,b,c)*((one_+u+v)*t-u*v*t*t*t/3)/((one_+u)*(one_+v));
    }

private:
    Number one_;
    std::array<Number,3> adjacent_;
    std::array<Number,4> c_, offsets_;
    Number sum_, signed_total_;

    // c_i/c_j for i<=j, with common vanishing factors removed.
    Number ratio(int i,int j) const {
        Number value=one_;
        for (; i<j; ++i) value=value*adjacent_[i];
        return value;
    }

    Number part(int numerator,int a,int b) const {
        return ratio(numerator,std::max(a,b))/(one_+ratio(std::min(a,b),std::max(a,b)));
    }

    int min_column(int cut) const { return cut==0?1:0; }

    Number leg(int cut,int corner,int other,Number t) const {
        return fraction(part(corner,corner,other)*(one_+ratio(min_column(cut),corner)*t));
    }
};

int exclude(const Box& box,Int target) {
    const Geometry<Range> g(box.variables,one);
    if (g.signed_total().hi<0) return 1;
    for (int i=0; i<4; ++i) {
        const auto e=g.cap_deficit(i);
        if (e.lo>8*target-scale/2 || e.hi<scale/2-7*target) return 1;
    }
    const auto areas=face_areas(
        [&](int face,int cut) { return g.triangle(face,cut); },
        [&](int face,int i,int j) { return g.overlap(face,i,j); });
    std::array<std::array<Range,4>,6> height{}, shift{};
    const auto total=nonnegative(g.weighted_total());
    for (unsigned f=0; f<faces.size(); ++f) {
        for (int i=0; i<4; ++i) {
            height[f][i]=one-intersect(g.point_height(i,faces[f]),{-scale,scale});
            shift[f][i]=g.shift(i,faces[f],total);
        }
    }
    for (int j=0; j<4; ++j) for (int mask=0; mask<16; ++mask) if (!(mask&(1<<j))) {
        Range sum{};
        for (unsigned f=0; f<faces.size(); ++f) {
            sum=sum+height[f][j]*fraction(areas[f][mask]+areas[f][mask|(1<<j)]);
        }
        if (sum.lo>(mask==(15^(1<<j))?6:12)*target) return 2;
    }
    const auto lower=cell_lower(areas,
        [&](unsigned f,int i) { return height[f][i]; },
        [&](unsigned f,int i,int j) { return shift[f][i]-shift[f][j]; });
    for (int mask=1; mask<15; ++mask) if (lower[mask]>target) return 3;
    return 0;
}

std::array<Linear<7>,14> affine_forms(const Box& box) {
    using Affine=Linear<7>;
    std::array<Affine,7> parameters{};
    for (int i=0;i<7;++i) parameters[i]=Affine::variable(box.variables[i],i);
    const auto unit=Affine::constant(scale);
    const Geometry<Affine> g(parameters,unit);
    const auto areas=face_expressions(unit,
        [&](int face,int cut) { return g.triangle(face,cut); },
        [&](int face,int i,int j) { return g.overlap(face,i,j); });
    std::array<std::array<Affine,4>,6> height{}, shift{};
    const auto total=g.weighted_total();
    for (unsigned f=0; f<faces.size(); ++f) for (int i=0; i<4; ++i) {
        height[f][i]=unit-g.point_height(i,faces[f]);
        shift[f][i]=g.shift(i,faces[f],total);
    }
    std::array<Affine,14> forms{};
    for (int mask=1; mask<15; ++mask) {
        std::array<int,3> omit{};
        int count=0;
        for (int i=0;i<4;++i) if (!(mask&(1<<i))) omit[count++]=i;
        const int i=omit[0],j=omit[1],k=omit[2];
        Affine result{};
        for (unsigned f=0; f<faces.size(); ++f) {
            result=result+height[f][i]*areas[f][mask];
            if (count>1) result=result+(shift[f][i]-shift[f][j])*areas[f][mask|(1<<i)];
            if (count==3) result=result+(shift[f][k]-shift[f][j])*areas[f][mask|(1<<i)|(1<<j)];
        }
        forms[mask-1]=result/6;
    }
    return forms;
}

bool candidate_neighborhood(const Box& box) {
    static const Int x_lower=[] {
        Int lower=0,upper=scale;
        while (upper-lower>1) {
            const Int mid=(lower+upper)/2;
            if (Wide{scale+mid}*(scale+mid)<Wide{2}*scale*scale) lower=mid;
            else upper=mid;
        }
        return lower;
    }();
    constexpr Int radius=scale/1600;
    const Range z=box.variables[0], y=z*box.variables[1], x=y*box.variables[2];
    if (z.lo<scale-radius || y.lo<scale-radius ||
        x.lo<x_lower+1-radius || x.hi>x_lower+radius) return false;
    for (int i=0;i<4;++i) {
        const Range offset=(i==0?y:x)*box.variables[3+i];
        if (offset.lo<-radius || offset.hi>radius) return false;
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

int main(int argc,char** argv) {
    const int depth=argc>1?parse<int>(argv[1]):49;
    const auto budget=argc>2?parse<std::uint64_t>(argv[2]):5000000;
    const int affine_depth=argc>3?parse<int>(argv[3]):depth+1;
    if (depth<0 || depth>126) return 2;
    const Int target=candidate_upper();
    Box root{};
    for (int i=0; i<7; ++i) root.variables[i]={i<3?0:-scale,scale};
    std::vector<Box> stack{root};
    std::array<std::uint64_t,6> excluded{};
    std::uint64_t visited=0,unresolved=0;
    std::array<Range,7> open{};
    while (!stack.empty() && visited<budget) {
        Box box=stack.back();
        stack.pop_back();
        ++visited;
        int verdict=exclude(box,target);
        if (!verdict && candidate_neighborhood(box)) verdict=5;
        if (!verdict && box.depth>=affine_depth && exclude_linear<7,14>(affine_forms(box),target)) verdict=4;
        if (verdict) { ++excluded[verdict]; continue; }
        if (box.depth==depth) {
            for (int i=0; i<7; ++i) {
                if (!unresolved) open[i]=box.variables[i];
                else {
                    open[i].lo=std::min(open[i].lo,box.variables[i].lo);
                    open[i].hi=std::max(open[i].hi,box.variables[i].hi);
                }
            }
            ++unresolved;
            continue;
        }
        const int axis=box.depth%7;
        const auto [lo,hi]=box.variables[axis];
        const Int midpoint=(lo+hi)/2;
        Box right=box;
        ++box.depth;
        ++right.depth;
        box.variables[axis].hi=midpoint;
        right.variables[axis].lo=midpoint;
        stack.push_back(right);
        stack.push_back(box);
    }
    std::cout << "visited=" << visited << " caps=" << excluded[1] << " three-cuts=" << excluded[2]
              << " volume=" << excluded[3] << " affine=" << excluded[4] << " local=" << excluded[5]
              << " unresolved=" << unresolved << " pending=" << stack.size() << '\n';
    if (unresolved) {
        for (const auto [lo,hi]:open) std::cout << '[' << lo << ',' << hi << "]/" << scale << ' ';
        std::cout << '\n';
    }
    return unresolved || !stack.empty();
}
