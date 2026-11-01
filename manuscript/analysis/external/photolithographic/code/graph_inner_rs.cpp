// Candidate-level GF(64) decoding. Evaluation references are opened only after
// every output decision has been written. Votes count graph paths, not reads.
#include "graph_inner_rs.h"
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
unsigned gf[64][64], powers[20]; int locations[64][64]; bool ready=false;
unsigned multiply(unsigned a,unsigned b){unsigned r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&64)a^=91;}return r;}
int base(char c){switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;default:return -1;}}
void initialize(){
    if(ready)return;
    for(unsigned a=0;a<64;a++)for(unsigned b=0;b<64;b++)gf[a][b]=multiply(a,b);
    powers[0]=1;for(int j=1;j<20;j++)powers[j]=gf[powers[j-1]][2];
    for(int a=0;a<64;a++)for(int b=0;b<64;b++)locations[a][b]=-1;
    for(unsigned a=1;a<64;a++)for(int j=0;j<20;j++)locations[a][gf[a][powers[j]]]=j;
    ready=true;
}
struct Key {
    uint64_t first=0,second=0;
    bool operator==(const Key &other)const{return first==other.first&&second==other.second;}
};
struct Hash{size_t operator()(const Key &k)const{return std::hash<uint64_t>()(k.first)^(std::hash<uint64_t>()(k.second)+0x9e3779b97f4a7c15ULL+(k.first<<6)+(k.first>>2));}};
struct Counts {uint32_t zero=0,one=0;};
using Map=std::unordered_map<Key,Counts,Hash>;
std::vector<Map> votes; std::string prefix; unsigned domain;
uint64_t paths,zero,one,invalid,outside,reassigned,position_counts[20];
Key pack(const char *s){Key k;for(int j=0;j<60;j++){auto &v=j<32?k.first:k.second;v=(v<<2)|base(s[j]);}return k;}
std::string unpack(Key k){std::string s(60,'A');for(int j=59;j>=0;j--){auto &v=j<32?k.first:k.second;s[j]="ACGT"[v&3];v>>=2;}return s;}
struct Policy{std::string name;unsigned minimum,zero_weight;bool priority,unique;};
std::vector<Policy> policies(){
    std::vector<Policy> p={{"strict",1,1,true,true},{"corrected_unique",1,1,true,true}};
    for(unsigned m:{1U,2U,3U,5U,10U,20U})p.push_back({"corrected_vote"+std::to_string(m),m,1,true,false});
    for(unsigned w:{1U,4U,16U})p.push_back({"pooled_weight"+std::to_string(w),1,w,false,false});
    return p;
}
}

extern "C" int inner_decode(const char *input,char *output,unsigned *index,int *position){
    initialize(); unsigned c[20],s0=0,s1=0;
    for(int j=0;j<20;j++){
        c[j]=0;for(int k=57-3*j;k<60-3*j;k++){int b=base(input[k]);if(b<0)return -1;c[j]=(c[j]<<2)|b;}
        s0^=c[j];s1^=gf[c[j]][powers[j]];
    }
    std::memcpy(output,input,60);output[60]=0;*position=-1;
    int errors=0;
    if(s0||s1){
        if(!s0)return -1;
        int j=locations[s0][s1];if(j<0)return -1;
        c[j]^=s0;*position=j;errors=1;
        for(int x=0;x<3;x++)output[59-3*j-x]="ACGT"[(c[j]>>(2*x))&3];
    }
    *index=0;for(int k=42;k<54;k++)*index=(*index<<2)|base(output[k]);*index^=1071068U;
    return errors;
}
extern "C" void inner_begin(unsigned n,const char *p){
    initialize();domain=n;prefix=p;votes.clear();votes.resize(n);
    paths=zero=one=invalid=outside=reassigned=0;std::memset(position_counts,0,sizeof(position_counts));
}
extern "C" int inner_add(const char *input,unsigned source_index){
    paths++;char corrected[61];unsigned index;int position;
    int e=inner_decode(input,corrected,&index,&position);
    if(e<0){invalid++;return e;}
    if(index>=domain){outside++;return -2;}
    auto &v=votes[index][pack(corrected)];
    if(e==0){v.zero++;zero++;}
    else{v.one++;one++;position_counts[position]++;if(index!=source_index)reassigned++;}
    return e;
}
extern "C" void inner_finish(const char *reference){
    auto ps=policies();std::vector<std::ofstream> outputs;
    for(const auto &p:ps)outputs.emplace_back(prefix+"."+p.name+".strands");
    std::ofstream decisions(prefix+".decisions.tsv");
    decisions<<"id\tzero_codewords\tone_codewords\tbest_one_votes\tbest_one_ties\n";
    std::vector<uint64_t> emitted(ps.size(),0);uint64_t distinct=0;
    for(unsigned i=0;i<domain;i++){
        const auto &m=votes[i];unsigned nz=0,no=0;uint64_t top=0,ties=0;Key zkey;
        for(const auto &p:m){if(p.second.zero){nz++;zkey=p.first;}if(p.second.one){no++;if(p.second.one>top){top=p.second.one;ties=1;}else if(p.second.one==top)ties++;}}
        distinct+=m.size();decisions<<i<<'\t'<<nz<<'\t'<<no<<'\t'<<top<<'\t'<<ties<<'\n';
        for(size_t j=0;j<ps.size();j++){
            const auto &p=ps[j];bool accept=false;Key winner;
            if(p.priority&&nz){if(nz==1){winner=zkey;accept=true;}}
            else if(p.name!="strict"){
                uint64_t maximum=0,ntie=0;unsigned eligible=0;
                for(const auto &v:m){uint64_t score=p.priority?v.second.one:uint64_t(v.second.zero)*p.zero_weight+v.second.one;
                    if(!score)continue;eligible++;
                    if(score>maximum){maximum=score;winner=v.first;ntie=1;}else if(score==maximum)ntie++;
                }
                accept=ntie==1&&maximum>=p.minimum&&(!p.unique||eligible==1);
            }
            if(accept){outputs[j]<<i<<'\t'<<unpack(winner)<<'\n';emitted[j]++;}
        }
    }
    for(auto &f:outputs)f.close();decisions.close();
    std::ofstream stats(prefix+".inner_stats.json");
    stats<<"{\"complete_paths\":"<<paths<<",\"zero_error_paths\":"<<zero<<",\"one_symbol_paths\":"<<one<<",\"uncorrectable_paths\":"<<invalid<<",\"outside_index_domain\":"<<outside<<",\"index_reassigned_paths\":"<<reassigned<<",\"distinct_decoded_codewords\":"<<distinct<<",\"one_symbol_position_counts\":[";
    for(int j=0;j<20;j++){if(j)stats<<',';stats<<position_counts[j];}stats<<"],\"emitted\":{";
    for(size_t j=0;j<ps.size();j++){if(j)stats<<',';stats<<'"'<<ps[j].name<<"\":"<<emitted[j];}stats<<"}}\n";stats.close();
    // The reference is deliberately unavailable to all decision code above.
    if(std::string(reference)!="-"){
        std::ifstream ref(reference);std::string line;std::vector<std::string> target;
        while(std::getline(ref,line))if(!line.empty()&&line[0]!='>')target.push_back(line);
        if(target.size()!=domain)throw std::runtime_error("reference length mismatch");
        std::ofstream diag(prefix+".truth_candidates.tsv");
        diag<<"id\tcorrect_zero_paths\tcorrect_one_paths\tcorrect_candidate_present\n";
        for(unsigned i=0;i<domain;i++){
            auto found=votes[i].find(pack(target[i].c_str()));Counts c;
            if(found!=votes[i].end())c=found->second;
            diag<<i<<'\t'<<c.zero<<'\t'<<c.one<<'\t'<<(c.zero||c.one)<<'\n';
        }
    }
    votes.clear();
}
