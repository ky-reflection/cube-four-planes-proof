#ifndef CUBE_AFFINE_BOX_LINEAR_HPP
#define CUBE_AFFINE_BOX_LINEAR_HPP

#include "affine_box_geometry.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>

namespace cube_affine_box {

template<int N> struct Linear {
    Int center{};
    std::array<Int,N> slope{};
    Range error{}, range{};

    static Linear constant(Int value) { return {value,{}, {}, {value,value}}; }
    static Linear variable(Range values,int axis) {
        const Int center=floor_ratio(Wide{values.lo}+values.hi,2);
        const Int radius=(values.hi-values.lo)/2;
        auto result=constant(center);
        result.slope[axis]=radius;
        result.error={values.lo-center+radius,values.hi-center-radius};
        result.range=values;
        return result;
    }
    Range linear_part() const {
        Int radius=0;
        for (Int a:slope) radius+=std::abs(a);
        return {-radius,radius};
    }
    void tighten() {
        const auto polynomial=Range{center,center}+linear_part();
        range=intersect(range,polynomial+error);
        error=intersect(error,range-polynomial);
        if (range.lo>range.hi || error.lo>error.hi) throw std::logic_error("inconsistent affine range");
    }
};

template<int N> Linear<N> operator+(const Linear<N>& a,const Linear<N>& b) {
    Linear<N> result{a.center+b.center,{},a.error+b.error,a.range+b.range};
    for (int i=0;i<N;++i) result.slope[i]=a.slope[i]+b.slope[i];
    return result;
}
template<int N> Linear<N> operator-(const Linear<N>& a) {
    Linear<N> result{-a.center,{},-a.error,-a.range};
    for (int i=0;i<N;++i) result.slope[i]=-a.slope[i];
    return result;
}
template<int N> Linear<N> operator-(const Linear<N>& a,const Linear<N>& b) { return a+-b; }

template<int N> Linear<N> operator*(const Linear<N>& a,const Linear<N>& b) {
    if ((a.range.lo==0 && a.range.hi==0) || (b.range.lo==0 && b.range.hi==0)) return {};
    Linear<N> result;
    result.center=floor_ratio(Wide{a.center}*b.center,scale);
    for (int i=0;i<N;++i) {
        result.slope[i]=floor_ratio(Wide{a.center}*b.slope[i]+Wide{b.center}*a.slope[i],scale);
    }
    const Range ac{a.center,a.center}, bc{b.center,b.center};
    result.error=Range{-N-1,N+1}+ac*b.error+bc*a.error+(a.range-ac)*(b.range-bc);
    result.range=a.range*b.range;
    result.tighten();
    return result;
}

template<int N> Linear<N> operator/(const Linear<N>& a,const Linear<N>& b) {
    Linear<N> result;
    result.center=floor_ratio(Wide{a.center}*scale,b.center);
    for (int i=0;i<N;++i) {
        result.slope[i]=floor_ratio(Wide{a.slope[i]}*scale-Wide{result.center}*b.slope[i],b.center);
    }
    const Int rounding=(N+1)*ceil_ratio(b.center,scale);
    result.error=(Range{-rounding,rounding}+a.error-Range{result.center,result.center}*b.error
                 -(b.range-Range{b.center,b.center})*result.linear_part())/b.range;
    result.range=a.range/b.range;
    result.tighten();
    return result;
}

template<int N> Linear<N> operator/(const Linear<N>& a,Int divisor) {
    Linear<N> result{floor_ratio(a.center,divisor),{},a.error/divisor+Range{-N-1,N+1},a.range/divisor};
    for (int i=0;i<N;++i) result.slope[i]=floor_ratio(a.slope[i],divisor);
    result.tighten();
    return result;
}

inline Range positive_square(Range a) {
    a=nonnegative(a);
    return a*a;
}

template<int N> Linear<N> positive_square(const Linear<N>& a) {
    if (a.range.hi<=0) return {};
    Linear<N> result;
    const Int c=std::max(Int{0},a.center);
    result.center=floor_ratio(Wide{c}*c,scale);
    for (int i=0;i<N;++i) result.slope[i]=floor_ratio(Wide{2*c}*a.slope[i],scale);
    const Int deviation=std::max(std::abs(a.range.lo-a.center),std::abs(a.range.hi-a.center));
    const Int quadratic=ceil_ratio(Wide{deviation}*deviation,scale);
    result.error=Range{-N-1,N+1}+Range{2*c,2*c}*a.error+Range{0,quadratic};
    result.range=positive_square(a.range);
    result.tighten();
    return result;
}

template<int N> Linear<N> fraction(Linear<N> a) {
    a.range=fraction(a.range);
    a.tighten();
    return a;
}

// Floating-point proposals for nonnegative weights on affine lower forms.
// Any proposal has to pass the exact integer inequality below.
template<int N,int K>
std::array<Int,K> propose_linear_weights(const std::array<Linear<N>,K>& forms,Int target) {
    constexpr int variables=K+N, constraints=1+2*N, columns=variables+constraints+1;
    std::array<std::array<double,columns>,constraints+1> tableau{};
    std::array<int,constraints> basic{};
    for (int row=0;row<constraints;++row) {
        tableau[row][variables+row]=1.;
        basic[row]=variables+row;
    }
    tableau[0].back()=1.;
    for (int i=0;i<K;++i) {
        tableau[0][i]=1.;
        tableau.back()[i]=(target-forms[i].center-forms[i].error.lo)/static_cast<double>(scale);
    }
    for (int j=0;j<N;++j) {
        tableau.back()[K+j]=1.;
        for (int sign : {-1,1}) {
            auto& row=tableau[1+2*j+(sign+1)/2];
            row[K+j]=-1.;
            for (int i=0;i<K;++i) row[i]=sign*forms[i].slope[j]/static_cast<double>(scale);
        }
    }
    for (int step=0;step<120;++step) {
        int enter=0;
        while (enter<columns-1 && tableau.back()[enter]>=-1e-11) ++enter;
        if (enter==columns-1) {
            std::array<Int,K> weights{};
            for (int row=0;row<constraints;++row) if (basic[row]<K) {
                const double value=tableau[row].back();
                if (!std::isfinite(value) || value>2.) return {};
                weights[basic[row]]=std::llround(std::max(0.,value)*65536.);
            }
            return weights;
        }
        int leave=-1;
        double best=std::numeric_limits<double>::infinity();
        for (int row=0;row<constraints;++row) if (tableau[row][enter]>1e-11) {
            const double ratio=std::max(0.,tableau[row].back())/tableau[row][enter];
            if (ratio<best-1e-11 || (leave>=0 && std::abs(ratio-best)<=1e-11 && basic[row]<basic[leave])) {
                leave=row;
                best=ratio;
            }
        }
        if (leave<0) return {};
        auto& pivot=tableau[leave];
        const double reciprocal=1./pivot[enter];
        for (double& entry:pivot) entry*=reciprocal;
        basic[leave]=enter;
        for (int row=0;row<=constraints;++row) if (row!=leave) {
            const double factor=tableau[row][enter];
            if (factor==0.) continue;
            for (int j=0;j<columns;++j) tableau[row][j]-=factor*pivot[j];
            tableau[row][enter]=0.;
        }
    }
    return {};
}

template<int N,int K> bool exclude_linear(const std::array<Linear<N>,K>& forms,Int target) {
    for (const auto& form:forms) if (form.range.lo>target) return true;
    const auto weights=propose_linear_weights<N,K>(forms,target);
    Wide lower=0;
    std::array<Wide,N> slope{};
    for (int i=0;i<K;++i) {
        lower+=Wide{weights[i]}*(forms[i].center+forms[i].error.lo-target);
        for (int j=0;j<N;++j) slope[j]+=Wide{weights[i]}*forms[i].slope[j];
    }
    for (const Wide a:slope) lower-=a<0?-a:a;
    return lower>0;
}

} // namespace cube_affine_box

#endif
