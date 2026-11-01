// Apply the authors' notebook's 55..70 nt filter, then the paper's 60 nt cap.
// Write two deterministic hash subsamples for exploratory validation.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include "../../kseq.h"
KSEQ_INIT(gzFile,gzread)
static uint64_t mix(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
int main(int argc,char **argv){
    if(argc!=3)return 1;gzFile f=gzopen(argv[1],"r");if(!f)return 2;
    kseq_t *ks=kseq_init(f);char fn[1024];
    snprintf(fn,sizeof(fn),"%s.fa.gz",argv[2]);gzFile out=gzopen(fn,"w1");
    snprintf(fn,sizeof(fn),"%s.sample.fa.gz",argv[2]);gzFile sample=gzopen(fn,"w1");
    uint64_t all=0,kept=0,n_sample=0,lens[301]={0};int ret;
    while((ret=kseq_read(ks))>=0){
        ++all;unsigned len=ks->seq.l;if(len<=300)++lens[len];
        if(len<55||len>70)continue;
        if(len>60)len=60;
        gzprintf(out,">%"PRIu64"\n%.*s\n",all,len,ks->seq.s);++kept;
        if(mix(all)%10==0){gzprintf(sample,">%"PRIu64"\n%.*s\n",all,len,ks->seq.s);++n_sample;}
    }
    if(ret!=-1){fprintf(stderr,"invalid input %d\n",ret);return 3;}
    int err;gzerror(f,&err);if(err!=Z_OK && err!=Z_STREAM_END)return 4;
    kseq_destroy(ks);gzclose(f);gzclose(out);gzclose(sample);
    printf("{\"total_reads\":%"PRIu64",\"filtered_reads\":%"PRIu64",\"sample_reads\":%"PRIu64",\"length_histogram\":{",all,kept,n_sample);
    int n=0;for(int i=0;i<=300;++i)if(lens[i])printf("%s\"%d\":%"PRIu64,n++?",":"",i,lens[i]);puts("}}");return 0;
}
