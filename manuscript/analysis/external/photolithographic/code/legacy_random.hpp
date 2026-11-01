// Deterministic reproduction of the libc++ uniform sampler used for synthesis.
// The historical source relied on a library-dependent uniform_int_distribution.
#pragma once
#include <cstdint>
#include <random>
template<class UInt> class LegacyUniform {
    UInt lower, span, mask;
public:
    LegacyUniform(UInt lo,UInt hi):lower(lo),span(hi-lo),mask(1){
        while(mask<span)mask=(mask<<1)|1;
    }
    UInt operator()(std::mt19937& rng){
        if(!span)return lower;
        UInt x;do{x=static_cast<UInt>(rng())&mask;}while(x>span);
        return lower+x;
    }
};
