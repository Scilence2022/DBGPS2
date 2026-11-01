// Minimal arithmetic operators used by the upstream RS templates.
// This replaces Boost's operator boilerplate, not the coding algorithms.
#pragma once
namespace boost {
struct empty_base {};
template<class T> struct equality_comparable {
    friend bool operator!=(const T& a,const T& b){return !(a==b);}
};
template<class T,class B=empty_base> struct ring_operators:B {
    friend T operator+(T a,const T& b){return a+=b;}
    friend T operator-(T a,const T& b){return a-=b;}
    friend T operator*(T a,const T& b){return a*=b;}
};
template<class T,class B=empty_base> struct field_operators:ring_operators<T,B> {
    friend T operator/(T a,const T& b){return a/=b;}
};
}
