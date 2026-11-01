// Stream minimap2 SAM into reference-assisted block FASTQs. Keep full original
// sequences/qualities, restoring reverse-strand SAM records before writing.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>

struct Alignment{int flag=0,block=-1,mapq=0,score=0,matches=0,query=0,reference=0,payload=0,columns=0;};
struct Read{std::string name,sequence,quality;std::vector<Alignment> alignments;bool primary=false;};
static bool broad(const Alignment&a){return a.block>=0&&a.matches>=20&&a.query>=30&&a.payload>=20&&a.score>=25&&a.matches*10>=a.columns*6;}
static bool strict(const Alignment&a){return broad(a)&&a.mapq>=10&&a.matches>=24&&a.query>=35&&a.payload>=25&&a.score>=40&&a.matches*10>=a.columns*7;}
static char complement(char c){switch(c){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';default:return 'N';}}
static void write_read(gzFile out,const Read&r){
    std::string text="@"+r.name+"\n"+r.sequence+"\n+\n"+r.quality+"\n";
    if(gzwrite(out,text.data(),(unsigned)text.size())!=(int)text.size()){std::cerr<<"FASTQ write failed\n";std::exit(3);}
}
int main(int argc,char**argv){
    if(argc!=2)return 1;std::string out=argv[1];
    std::array<gzFile,12> bs{},ss{};std::array<uint64_t,12> nb{},ns{};
    for(int b=0;b<12;b++){
        char name[64];std::snprintf(name,sizeof(name),"/broad/block%02d.fq.gz",b+1);bs[b]=gzopen((out+name).c_str(),"w1");
        std::snprintf(name,sizeof(name),"/strict/block%02d.fq.gz",b+1);ss[b]=gzopen((out+name).c_str(),"w1");
        if(!bs[b]||!ss[b])return 2;
    }
    gzFile ambiguous=gzopen((out+"/ambiguous.fq.gz").c_str(),"w1"),unassigned=gzopen((out+"/unassigned.fq.gz").c_str(),"w1"),inclusive=gzopen((out+"/block12.inclusive.fq.gz").c_str(),"w1");
    if(!ambiguous||!unassigned||!inclusive)return 2;
    uint64_t total=0,na=0,nu=0,ni=0,alignment_records=0,secondary_records=0,reverse_restored=0,strict_rejected=0;
    Read current;
    auto finish=[&](){
        if(current.name.empty())return;
        if(!current.primary||current.sequence.empty()||current.sequence.size()!=current.quality.size()){std::cerr<<"Missing primary sequence: "<<current.name<<"\n";std::exit(4);}
        total++;int best=-1,best_score=-1;unsigned mask=0;
        for(size_t j=0;j<current.alignments.size();j++){const auto&a=current.alignments[j];if(broad(a)&&a.score>best_score){best=(int)j;best_score=a.score;}}
        if(best>=0)for(const auto&a:current.alignments)if(broad(a)&&a.score*5>=best_score*4)mask|=1U<<a.block;
        bool unique=mask&&!(mask&(mask-1));
        if(unique){
            int b=__builtin_ctz(mask);write_read(bs[b],current);nb[b]++;
            bool accept=false;for(const auto&a:current.alignments)if(!(a.flag&2304)&&a.block==b&&strict(a))accept=true;
            if(accept){write_read(ss[b],current);ns[b]++;}else strict_rejected++;
        }else if(mask){write_read(ambiguous,current);na++;}
        else{write_read(unassigned,current);nu++;}
        if(mask&(1U<<11)){write_read(inclusive,current);ni++;}
        if(total%5000000==0)std::cerr<<"READS "<<total<<" broad12 "<<nb[11]<<" strict12 "<<ns[11]<<" ambiguous "<<na<<"\n";
        current=Read();
    };
    char*line=nullptr;size_t cap=0;ssize_t len;
    while((len=getline(&line,&cap,stdin))>=0){
        if(line[0]=='@')continue;
        std::vector<char*>f;char*save=nullptr;
        for(char*p=strtok_r(line,"\t\r\n",&save);p;p=strtok_r(nullptr,"\t\r\n",&save))f.push_back(p);
        if(f.size()<11)return 5;
        if(current.name!=f[0]){finish();current.name=f[0];}
        Alignment a;a.flag=std::atoi(f[1]);a.mapq=std::atoi(f[4]);alignment_records++;
        if(a.flag&256)secondary_records++;
        if(!(a.flag&2304)){
            if(current.primary)return 6;current.primary=true;current.sequence=f[9];current.quality=f[10];
            if(a.flag&16){
                std::reverse(current.sequence.begin(),current.sequence.end());for(char&c:current.sequence)c=complement(c);
                std::reverse(current.quality.begin(),current.quality.end());reverse_restored++;
            }
        }
        if(a.flag&4)continue;
        int id=std::atoi(f[2]);if(id<0||id>=196596)return 7;a.block=id/16383;
        int n=0;for(char*p=f[5];*p;p++){
            if(*p>='0'&&*p<='9'){n=n*10+*p-'0';continue;}
            if(*p=='=')a.matches+=n;
            if(*p=='='||*p=='X'||*p=='M'||*p=='I')a.query+=n;
            if(*p=='='||*p=='X'||*p=='M'||*p=='D')a.reference+=n;
            if(*p=='='||*p=='X'||*p=='M'||*p=='I'||*p=='D')a.columns+=n;
            n=0;
        }
        int start=std::atoi(f[3])-1;a.payload=std::max(0,std::min(42,start+a.reference)-std::max(0,start));
        for(size_t j=11;j<f.size();j++)if(!std::strncmp(f[j],"AS:i:",5))a.score=std::atoi(f[j]+5);
        current.alignments.push_back(a);
    }
    finish();std::free(line);
    for(auto f:bs)if(gzclose(f)!=Z_OK)return 8;
    for(auto f:ss)if(gzclose(f)!=Z_OK)return 8;
    if(gzclose(ambiguous)!=Z_OK||gzclose(unassigned)!=Z_OK||gzclose(inclusive)!=Z_OK)return 8;
    std::ofstream stats(out+"/split_counts.json");
    stats<<"{\"total_reads\":"<<total<<",\"sam_alignment_records\":"<<alignment_records<<",\"secondary_records\":"<<secondary_records
         <<",\"reverse_restored_reads\":"<<reverse_restored<<",\"ambiguous_reads\":"<<na<<",\"unassigned_reads\":"<<nu
         <<",\"block12_inclusive_reads\":"<<ni<<",\"strict_rejected_from_broad\":"<<strict_rejected<<",\"blocks\":[";
    uint64_t assigned=0;
    for(int b=0;b<12;b++){if(b)stats<<',';assigned+=nb[b];stats<<"{\"block\":"<<b+1<<",\"broad_reads\":"<<nb[b]<<",\"strict_reads\":"<<ns[b]<<'}';}
    stats<<"]}\n";stats.close();assert(total==assigned+na+nu);
    std::cerr<<"FINISHED reads "<<total<<" broad_assigned "<<assigned<<" ambiguous "<<na<<" unassigned "<<nu<<"\n";
}
