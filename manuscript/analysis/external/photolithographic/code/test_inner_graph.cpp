// Verify all 63 nonzero GF(64) symbol errors at all 20 positions against
// independently encoded public reference strands, including index errors.
#include "graph_inner_rs.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
static unsigned gf(unsigned a,unsigned b){unsigned r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&64)a^=91;}return r;}
static std::string damage(std::string s,int j,unsigned error){
    unsigned x=0;for(int k=57-3*j;k<60-3*j;k++)x=(x<<2)|std::string("ACGT").find(s[k]);x^=error;
    for(int k=0;k<3;k++)s[59-3*j-k]="ACGT"[(x>>(2*k))&3];return s;
}
static std::string alternate(std::string s){
    unsigned a=1;for(int j=0;j<6;j++)a=gf(a,2);
    unsigned inverse=1;while(gf(inverse,3)!=1)inverse++;
    unsigned v=gf(a^1,inverse),u=v^1;
    return damage(damage(damage(s,6,1),1,v),0,u);
}
static std::string read(const std::string &p){std::ifstream f(p);return std::string(std::istreambuf_iterator<char>(f),{});}
int main(int argc,char **argv){
    if(argc!=2&&argc!=3)return 1;
    std::ifstream f(argv[1]);std::string s;unsigned checks=0,originals=0;
    std::vector<std::string> samples;
    while(originals<100&&std::getline(f,s)){
        if(s.size()!=60||s[0]=='>')continue;
        char out[61];unsigned index;int position;
        if(inner_decode(s.c_str(),out,&index,&position)!=0||s!=out)return 2;
        samples.push_back(s);
        for(int j=0;j<20;j++)for(unsigned error=1;error<64;error++){
            std::string damaged=s;unsigned symbol=0;
            for(int k=57-3*j;k<60-3*j;k++)symbol=(symbol<<2)|std::string("ACGT").find(s[k]);
            symbol^=error;for(int k=0;k<3;k++)damaged[59-3*j-k]="ACGT"[(symbol>>(2*k))&3];
            unsigned repaired_index;int repaired_position;
            if(inner_decode(damaged.c_str(),out,&repaired_index,&repaired_position)!=1||s!=out||repaired_position!=j||repaired_index!=index)return 3;
            checks++;
        }
        originals++;
    }
    if(originals!=100)return 4;
    std::cout<<checks<<" single-symbol errors corrected exactly across "<<originals<<" public codewords\n";
    if(argc==3){
        std::string root=argv[2];auto a1=alternate(samples[1]),a3=alternate(samples[3]);
        char out[61];unsigned idx;int pos;
        if(inner_decode(a1.c_str(),out,&idx,&pos)!=0||idx!=1||a1==samples[1])return 5;
        for(int blind=0;blind<2;blind++){
            std::string prefix=root+(blind?".blind":".reference");inner_begin(16383,prefix.c_str());
            // Source 1 has an index error; decoding must reassign it to 0.
            inner_add(damage(samples[0],2,1).c_str(),1);
            inner_add(samples[1].c_str(),1);
            for(int j=0;j<3;j++)inner_add(damage(a1,7,1+j).c_str(),1);
            inner_add(damage(samples[2],6,1).c_str(),2);inner_add(damage(samples[2],7,1).c_str(),2);
            inner_add(damage(samples[3],6,1).c_str(),3);inner_add(damage(a3,7,1).c_str(),3);
            inner_finish(blind?"-":argv[1]);
        }
        for(const std::string &policy:{"strict","corrected_unique","corrected_vote1","corrected_vote2","corrected_vote3","corrected_vote5","corrected_vote10","corrected_vote20","pooled_weight1","pooled_weight4","pooled_weight16"})
            if(read(root+".blind."+policy+".strands")!=read(root+".reference."+policy+".strands"))return 6;
        std::string expected="0\t"+samples[0]+"\n1\t"+samples[1]+"\n2\t"+samples[2]+"\n";
        if(read(root+".blind.corrected_vote1.strands")!=expected)return 7;
        if(read(root+".blind.corrected_vote2.strands")!="1\t"+samples[1]+"\n2\t"+samples[2]+"\n")return 8;
        if(read(root+".blind.pooled_weight1.strands")!="0\t"+samples[0]+"\n1\t"+a1+"\n2\t"+samples[2]+"\n")return 9;
        if(read(root+".blind.pooled_weight4.strands")!=expected)return 10;
        std::cout<<"Index reassignment, vote thresholds, ties, weights, and reference-blind output checks passed\n";
    }
}
