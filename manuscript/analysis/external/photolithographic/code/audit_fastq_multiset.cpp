// Order-independent SHA-256 record multiset audit (macOS CommonCrypto).
#include <array>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <zlib.h>
#include <CommonCrypto/CommonDigest.h>
#include "../../kseq.h"
KSEQ_INIT(gzFile,gzread)
int main(int argc,char**argv){
    if(argc<2)return 1;uint64_t count=0,bases=0;std::array<uint64_t,4> sum{},xors{};
    for(int j=1;j<argc;j++){
        gzFile f=gzopen(argv[j],"rb");if(!f)return 2;kseq_t*ks=kseq_init(f);int ret;
        while((ret=kseq_read(ks))>=0){
            if(ks->seq.l!=ks->qual.l)return 3;
            CC_SHA256_CTX ctx;unsigned char hash[32];const unsigned char zero=0;
            CC_SHA256_Init(&ctx);CC_SHA256_Update(&ctx,ks->name.s,(CC_LONG)ks->name.l);CC_SHA256_Update(&ctx,&zero,1);
            CC_SHA256_Update(&ctx,ks->seq.s,(CC_LONG)ks->seq.l);CC_SHA256_Update(&ctx,&zero,1);
            CC_SHA256_Update(&ctx,ks->qual.s,(CC_LONG)ks->qual.l);CC_SHA256_Final(hash,&ctx);
            for(int k=0;k<4;k++){
                uint64_t v=0;for(int b=0;b<8;b++)v=(v<<8)|hash[8*k+b];sum[k]+=v;xors[k]^=v;
            }
            count++;bases+=ks->seq.l;
        }
        int err;gzerror(f,&err);if(ret!=-1||(err!=Z_OK&&err!=Z_STREAM_END))return 4;
        kseq_destroy(ks);if(gzclose(f)!=Z_OK)return 5;
    }
    std::cout<<"{\"reads\":"<<count<<",\"bases\":"<<bases<<",\"sha256_record_sum\":\"";
    for(auto v:sum)std::cout<<std::hex<<std::setfill('0')<<std::setw(16)<<v;
    std::cout<<"\",\"sha256_record_xor\":\"";
    for(auto v:xors)std::cout<<std::hex<<std::setfill('0')<<std::setw(16)<<v;
    std::cout<<"\"}\n";
}
