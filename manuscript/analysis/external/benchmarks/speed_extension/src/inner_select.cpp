// Shared inner RS(20,18)/GF(64) selection for fair candidate comparisons.
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
static unsigned table[64][64], powers[63];
static unsigned multiply(unsigned a,unsigned b){unsigned r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&64)a^=91;}return r;}
static int base(char c){auto p=std::string("ACGT").find(c);return p==std::string::npos?-1:int(p);}
struct Votes { int minimum=2;std::map<std::string,unsigned> counts; };
static int candidate(std::string &s,unsigned n,unsigned &index){
    std::array<unsigned,20> c;unsigned s0=0,s1=0;
    for(int j=0;j<20;j++){
        c[j]=0;for(int k=57-3*j;k<60-3*j;k++){int b=base(s[k]);if(b<0)return -1;c[j]=(c[j]<<2)|b;}
        s0^=c[j];s1^=table[c[j]][powers[j]];
    }
    int errors=0;
    if(s0 || s1){
        if(!s0)return -1;
        int position=-1;
        for(int j=0;j<20;j++)if(table[s0][powers[j]]==s1){position=j;break;}
        if(position<0)return -1;
        c[position]^=s0;errors=1;
        for(int j=0;j<20;j++)for(int k=0;k<3;k++)s[59-3*j-k]="ACGT"[(c[j]>>(2*k))&3];
    }
    index=0;for(int k=42;k<54;k++)index=(index<<2)|base(s[k]);index^=1071068U;
    return index<n?errors:-1;
}
static unsigned write(const std::vector<Votes>& votes,const std::string &path,int max_errors){
    std::ofstream out(path);unsigned emitted=0;
    for(unsigned i=0;i<votes.size();i++){
        if(votes[i].minimum>max_errors)continue;
        unsigned best=0,ties=0;std::string winner;
        for(const auto &p:votes[i].counts){
            if(p.second>best){best=p.second;ties=1;winner=p.first;}
            else if(p.second==best)ties++;
        }
        if(ties==1){out<<i<<'\t'<<winner<<'\n';emitted++;}
    }
    return emitted;
}
int main(int argc,char **argv){
    if(argc!=4)return 1;
    for(unsigned i=0;i<64;i++)for(unsigned j=0;j<64;j++)table[i][j]=multiply(i,j);
    powers[0]=1;for(int i=1;i<63;i++)powers[i]=table[powers[i-1]][2];
    unsigned n=std::stoul(argv[2]);std::vector<Votes> votes(n);std::ifstream f(argv[1]);std::string line;
    unsigned long long reads=0,accepted=0;
    while(std::getline(f,line)){
        auto tab=line.find('\t');if(tab!=std::string::npos)line=line.substr(tab+1);
        if(line.size()!=60)continue;reads++;
        for(int ori=0;ori<2;ori++){
            std::string s=line;
            if(ori)for(int k=0;k<60;k++){int b=base(line[59-k]);s[k]=b<0?'N':"TGCA"[b];}
            unsigned index;int errors=candidate(s,n,index);if(errors<0)continue;accepted++;
            auto &v=votes[index];if(errors<v.minimum){v.minimum=errors;v.counts.clear();}
            if(errors==v.minimum)v.counts[s]++;
        }
    }
    unsigned strict=write(votes,std::string(argv[3])+".strict.strands",0);
    unsigned corrected=write(votes,std::string(argv[3])+".corrected.strands",1);
    std::cout<<"{\"candidates\":"<<reads<<",\"accepted_orientations\":"<<accepted<<",\"strict_emitted\":"<<strict<<",\"corrected_emitted\":"<<corrected<<"}\n";
}
