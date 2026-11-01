// Algebraically identical GF(2^14), using logarithms instead of bitwise products.
#pragma once
#include <array>
#include <bitset>
#include <ostream>
class FastGF14:boost::field_operators<FastGF14,boost::equality_comparable<FastGF14>> {
    struct Tables {
        std::array<unsigned,32766> exp{};
        std::array<unsigned,16384> log{};
        static unsigned multiply(unsigned a,unsigned b){unsigned r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&16384)a^=16553;}return r;}
        Tables(){unsigned x=1;for(unsigned i=0;i<16383;i++){exp[i]=x;log[x]=i;x=multiply(x,66);}for(unsigned i=16383;i<32766;i++)exp[i]=exp[i-16383];}
    };
    static const Tables& tables(){static const Tables t;return t;}
public:
    using uint=unsigned long long;
    static const unsigned m=14;
    static const unsigned prim_poly=16553;
    uint el=0;bool empty=true;
    FastGF14()=default;
    FastGF14(uint x):el(x%2),empty(false){}
    FastGF14(uint x,unsigned):el(x),empty(false){}
    FastGF14& operator+=(const FastGF14& b){el^=b.el;return *this;}
    FastGF14& operator-=(const FastGF14& b){el^=b.el;return *this;}
    FastGF14& operator*=(const FastGF14& b){if(!el||!b.el)el=0;else{const auto&t=tables();el=t.exp[t.log[el]+t.log[b.el]];}return *this;}
    FastGF14& operator/=(const FastGF14& b){return *this*=b.inverse();}
    FastGF14 inverse()const{if(!el)return FastGF14(0);const auto&t=tables();return FastGF14(t.exp[16383-t.log[el]],0);}
    bool operator==(const FastGF14& b)const{return el==b.el;}
    bool operator<(const FastGF14& b)const{return el<b.el;}
    bool iszero()const{return el==0;}
    bool isempty()const{return empty;}
    FastGF14 power(int n)const{if(!el)return *this;const auto&t=tables();long long e=(long long)t.log[el]*n%16383;if(e<0)e+=16383;return FastGF14(t.exp[e],0);}
};
inline FastGF14 pow(const FastGF14&a,int n){return a.power(n);}
inline std::ostream& operator<<(std::ostream&o,const FastGF14&a){return o<<std::bitset<32>(a.el);}
