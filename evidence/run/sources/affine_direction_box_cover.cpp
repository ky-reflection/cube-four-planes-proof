// A direct offset cover for independently varying directions in one box.
// Build: c++ -O3 -std=c++20 affine_direction_box_cover.cpp -o /tmp/direction-cover
// Usage: direction-cover [width in thousandths] [depth] [node budget]

#include <charconv>
#include <iostream>
#include <string_view>
#include <vector>

#include "affine_box_linear.hpp"

namespace {
using namespace cube_affine_box;
using Scores = std::array<std::array<Range,4>,4>;
using Affine = Linear<4>;
struct Box { std::array<Range,4> offsets; int depth; };

Affine uncertain_constant(Range value) {
    const Int center=floor_ratio(Wide{value.lo}+value.hi,2);
    return {center,{}, {value.lo-center,value.hi-center},value};
}

std::array<int,3> other_labels(int i) {
    std::array<int,3> labels{};
    int used=0;
    for (int j=0; j<4; ++j) if (j!=i) labels[used++]=j;
    return labels;
}

class DirectionBox {
public:
    explicit DirectionBox(Scores c) : c_(c) {
        for (int i=0; i<4; ++i) {
            const auto [j,k,l]=other_labels(i);
            row_sum_[i]=c_[i][j]+c_[i][k]+c_[i][l];
            cap_bound_[i]=std::min({c_[i][j].hi,c_[i][k].hi,c_[i][l].hi});
            const auto a=c_[i][j], b=c_[i][k], d=c_[i][l];
            const auto den=(a+b)*(a+d)*(b+d);
            cap_linear_[i]=(a*b+a*d+b*d)/den;
            cap_cubic_[i]=one/(den*Range{3*scale,3*scale});
        }
        for (int j=0; j<4; ++j) {
            weights_[j]=tree_weight(j);
            const auto labels=other_labels(j);
            std::array<std::array<Range,4>,4> inverse{};
            for (int i:labels) for (int k:labels) {
                Range numerator{};
                if (i==k) {
                    const int a=next_label(labels,i), b=6-j-i-a;
                    numerator=(c_[a][i]+c_[a][j])*row_sum_[b]+c_[a][b]*(c_[b][i]+c_[b][j]);
                } else {
                    const int l=6-j-i-k;
                    numerator=c_[i][k]*row_sum_[l]+c_[i][l]*c_[l][k];
                }
                inverse[i][k]=numerator/weights_[j];
            }
            for (unsigned f=0; f<faces.size(); ++f) for (int k:labels) {
                Range coefficient{};
                for (int i:labels) {
                    coefficient=coefficient+(faces[f]&(1<<i)?inverse[i][k]:-inverse[i][k]);
                }
                height_coefficients_[f][j][k]=uncertain_constant(coefficient);
            }
        }
        for (int i=0;i<4;++i) for (int j=0;j<4;++j) affine_c_[i][j]=uncertain_constant(c_[i][j]);
    }

    Box root() const {
        Box result{};
        for (int i=0;i<4;++i) result.offsets[i]={-cap_bound_[i],cap_bound_[i]};
        return result;
    }

    int exclude(const Box& box,Int target) const {
        Range total{};
        for (int i=0;i<4;++i) {
            const auto h=box.offsets[i];
            total=total+weights_[i]*h;
            const auto deficit=cap_linear_[i]*h-cap_cubic_[i]*h*h*h;
            if (deficit.lo>8*target-scale/2 || deficit.hi<scale/2-7*target) return 1;
        }
        if (total.hi<0) return 1;
        return exclude_linear<4,14>(forms(box),target)?2:0;
    }

private:
    Scores c_;
    std::array<Range,4> row_sum_{},weights_{},cap_linear_{},cap_cubic_{};
    std::array<Int,4> cap_bound_{};
    std::array<std::array<Affine,4>,4> affine_c_{};
    std::array<std::array<std::array<Affine,4>,4>,6> height_coefficients_{};

    static int next_label(const std::array<int,3>& labels,int except) {
        for (int a:labels) if (a!=except) return a;
        throw std::logic_error("empty minor");
    }

    Range tree_weight(int root) const {
        const auto [a,b,c]=other_labels(root);
        Range result{};
        std::array<int,4> parent{};
        parent[root]=root;
        for (int i:other_labels(a)) for (int j:other_labels(b)) for (int k:other_labels(c)) {
            parent[a]=i;parent[b]=j;parent[c]=k;
            bool tree=true;
            for (int start:{a,b,c}) {
                int at=start;
                for (int step=0;step<3;++step) at=parent[at];
                if (at!=root) tree=false;
            }
            if (tree) result=result+c_[a][i]*c_[b][j]*c_[c][k];
        }
        if (result.lo<=0) throw std::logic_error("nonpositive tree weight");
        return result;
    }

    std::array<Affine,14> forms(const Box& box) const {
        const auto unit=Affine::constant(scale);
        std::array<Affine,4> h{};
        for (int i=0;i<4;++i) h[i]=Affine::variable(box.offsets[i],i);
        const auto& c=affine_c_;
        const auto areas=face_expressions(unit,
            [&](int face,int i) {
                const int corner=i^face_shift(face);
                const auto numerator=c[i][corner]+(face&(1<<i)?h[i]:-h[i]);
                Affine area=unit;
                for (int j=0;j<4;++j) if (j!=i && j!=corner) {
                    area=area*fraction(numerator/(c[i][corner]+c[i][j]));
                }
                return area/2;
            },
            [&](int face,int i,int j) {
                const int a=i^face_shift(face),b=j^face_shift(face);
                const auto denom_i=c[i][a]+c[i][b],denom_j=c[j][a]+c[j][b];
                const auto first=fraction((c[i][a]+h[i])/denom_i);
                const auto second=fraction((c[j][b]-h[j])/denom_j);
                const auto d=(c[i][a]+c[i][j])/denom_i+(c[j][b]+c[j][i])/denom_j;
                return fraction(positive_square(first+second-unit)/d/2);
            });
        std::array<std::array<Affine,4>,6> point{};
        for (unsigned f=0;f<faces.size();++f) for (int j=0;j<4;++j) {
            for (int i=0;i<4;++i) point[f][j]=point[f][j]+height_coefficients_[f][j][i]*h[i];
        }
        std::array<Affine,14> result{};
        for (int mask=1;mask<15;++mask) {
            std::array<int,3> omitted{};
            int count=0;
            for (int i=0;i<4;++i) if (!(mask&(1<<i))) omitted[count++]=i;
            const auto [i,j,k]=omitted;
            Affine volume{};
            for (unsigned f=0;f<faces.size();++f) {
                volume=volume+(unit-point[f][i])*areas[f][mask];
                if (count>1) volume=volume+(point[f][j]-point[f][i])*areas[f][mask|(1<<i)];
                if (count==3) volume=volume+(point[f][j]-point[f][k])*areas[f][mask|(1<<i)|(1<<j)];
            }
            result[mask-1]=volume/6;
        }
        return result;
    }
};

template<class T> T parse(std::string_view arg) {
    T value{};
    const auto [end,error]=std::from_chars(arg.data(),arg.data()+arg.size(),value);
    if (error!=std::errc{} || end!=arg.data()+arg.size()) throw std::invalid_argument("invalid integer");
    return value;
}

} // namespace

int main(int argc,char** argv) {
    const int width=argc>1?parse<int>(argv[1]):5;
    const int depth=argc>2?parse<int>(argv[2]):52;
    const auto budget=argc>3?parse<std::uint64_t>(argv[3]):1000000;
    if (width<0 || width>100 || depth<0 || depth>80) return 2;
    const std::array<std::array<int,4>,4> centers{{{0,1000,620,850},{1000,0,550,930},
                                                   {820,630,0,1000},{580,910,1000,0}}};
    Scores scores{};
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) if (i!=j) {
        const int c=centers[i][j],w=c==1000?0:width;
        scores[i][j]={floor_ratio(Wide{c-w}*scale,1000),ceil_ratio(Wide{c+w}*scale,1000)};
    }
    const DirectionBox direction(scores);
    const Int target=candidate_upper();
    std::vector<Box> stack{direction.root()};
    std::array<std::uint64_t,3> closed{};
    std::uint64_t visited=0,unresolved=0;
    while (!stack.empty() && visited<budget) {
        Box box=stack.back();
        stack.pop_back();
        ++visited;
        const int verdict=direction.exclude(box,target);
        if (verdict) { ++closed[verdict]; continue; }
        if (box.depth==depth) { ++unresolved; continue; }
        const int axis=box.depth%4;
        const auto [lo,hi]=box.offsets[axis];
        const Int midpoint=floor_ratio(Wide{lo}+hi,2);
        Box right=box;
        ++box.depth;
        ++right.depth;
        box.offsets[axis].hi=midpoint;
        right.offsets[axis].lo=midpoint;
        stack.push_back(right);
        stack.push_back(box);
    }
    std::cout << "width=" << width << "/1000 visited=" << visited << " caps=" << closed[1]
              << " volumes=" << closed[2] << " unresolved=" << unresolved << " pending=" << stack.size() << '\n';
    return unresolved || !stack.empty();
}
