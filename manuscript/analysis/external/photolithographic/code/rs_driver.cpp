// Portable command-line entry point around the authors' unchanged RS algorithms.
#include "operators_compat.hpp"
#include "legacy_random.hpp"
#ifdef FAST_GF14
#include "fast_gf14.hpp"
#endif
#include <bitset>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include <random>
#include "GF2M.hpp"
#include "DFT.hpp"
#include "ReedSolomon.hpp"
#include "encodedecode.hpp"
int main(int argc,char **argv){
    if(argc<5){std::cerr<<"usage: rs_driver encode|decode input output numblocks\n";return 1;}
    typedef unsigned long long uint;
    typedef GF2M<uint,6,91> GFI;
#ifdef FAST_GF14
    typedef FastGF14 GFO;
#else
    typedef GF2M<uint,14,16553> GFO;
#endif
    DFT_FFT<GFI> di(63,GFI(2,0),7,9);
    DFT_FFT<GFO> d_o(16383,GFO(66,0),129,127);
    typedef RScode<GFI,DFT_FFT<GFI>> IC;
    typedef RScode<GFO,DFT_FFT<GFO>> OC;
    IC ic(20,18,GFI(2,0),di,63);
    OC oc(16383,10977,GFO(66,0),d_o,16383);
    EnDecode<IC,OC> ed(ic,oc,4,6);ed.numblocks=std::stoi(argv[4]);
    std::ifstream f(argv[2],std::ios::binary);
    if(!f)return 2;
    std::vector<std::string> seq;
    std::string content;
    if(std::string(argv[1])=="encode") {
        content.assign(std::istreambuf_iterator<char>(f),{});ed.encode(content,seq);
        std::ofstream out(argv[3]);for(auto &s:seq)out<<s<<'\n';
    } else {
        std::string line;while(std::getline(f,line)){
            auto tab=line.find('\t');if(tab!=std::string::npos)line=line.substr(tab+1);
            if(line.size()==60 && line.find_first_not_of("ACGT")==std::string::npos)seq.push_back(line);
        }
        ed.decode(content,seq,"");
        std::ofstream out(argv[3],std::ios::binary);out.write(content.data(),content.size());
    }
}
