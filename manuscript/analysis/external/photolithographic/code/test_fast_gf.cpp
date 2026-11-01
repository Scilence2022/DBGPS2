#include "operators_compat.hpp"
#include <cassert>
#include <cmath>
#include <bitset>
#include <random>
#include "GF2M.hpp"
#include "fast_gf14.hpp"
int main(){using R=GF2M<unsigned long long,14,16553>;unsigned checks=0;for(unsigned a=0;a<16384;a++){for(unsigned b:{0U,1U,2U,3U,66U,8192U,16383U}){assert((R(a,0)*R(b,0)).el==(FastGF14(a,0)*FastGF14(b,0)).el);checks++;}assert(R(a,0).inverse().el==FastGF14(a,0).inverse().el);checks++;}std::mt19937 r(20261017);for(int i=0;i<100000;i++){unsigned a=r()%16384,b=r()%16384;assert((R(a,0)*R(b,0)).el==(FastGF14(a,0)*FastGF14(b,0)).el);checks++;}for(int i=0;i<1000;i++){unsigned a=r()%16384;int e=int(r()%257)-128;assert(pow(R(a,0),e).el==pow(FastGF14(a,0),e).el);checks++;}std::cout<<checks<<" GF equivalence checks passed\n";}
