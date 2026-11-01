/* ORACLE DIAGNOSTIC ONLY. Uses references to test whether their exact paths
 * exist in the read-derived graph; never emits or repairs a recovered strand. */
#define main original_dbgps2_main
#include "../../../../DBGPS2.c"
#undef main

static int compare_u64(const void *a, const void *b)
{
    uint64_t x=*(const uint64_t *)a, y=*(const uint64_t *)b;
    return x<y ? -1 : x>y;
}

static unsigned present_path(kc_c4x_t *h, const char *seq, int k, uint64_t mask,
                             const uint64_t *seeds, size_t nseeds)
{
    uint64_t seed31=actgkmer_uint64((unsigned char *)seq+3,31);
    if (!bsearch(&seed31,seeds,nseeds,sizeof(uint64_t),compare_u64)) return 0;
    uint64_t word=actgkmer_uint64((unsigned char *)seq+34-k,k);
    unsigned minimum=KC_MAX;
    for(int pos=34;pos<130;++pos){
        word=((word<<2)|seq_nt4_table[(unsigned char)seq[pos]])&mask;
        unsigned cov=kmer_cov(min_hash_key(word,k),mask,h);
        if(cov<minimum)minimum=cov;
        if(!minimum)return 0;
    }
    unsigned best_child=0;
    for(unsigned base=0;base<4;++base){
        uint64_t next=((word<<2)|base)&mask;
        unsigned cov=kmer_cov(min_hash_key(next,k),mask,h);
        if(cov>best_child)best_child=cov;
    }
    if(best_child<minimum)minimum=best_child;
    unsigned bits=0;
    for(unsigned cutoff=0;cutoff<7;++cutoff)if(minimum>cutoff)bits|=1U<<cutoff;
    return bits;
}

int main(int argc,char **argv)
{
    if(argc!=5 && argc!=6)return 1;
    int k=atoi(argv[3]);uint64_t mask=(1ULL<<(2*k))-1;
    FILE *sf=fopen(argv[2],"r");if(!sf)return 1;
    size_t ns=0,cap=4000000;uint64_t *seeds=malloc(cap*sizeof(*seeds));
    char dir,seed[64];unsigned support;
    while(fscanf(sf," %c %63s %u",&dir,seed,&support)==3){
        if(ns==cap){cap*=2;seeds=realloc(seeds,cap*sizeof(*seeds));}
        seeds[ns++]=actgkmer_uint64((unsigned char *)seed,31);
    }
    fclose(sf);qsort(seeds,ns,sizeof(*seeds),compare_u64);
    kc_c4x_t *h=count_file(argv[1],k,KC_BITS,10000000,4,200);if(!h)return 1;
    gzFile fp=gzopen(argv[4],"r");if(!fp)return 1;kseq_t *ks=kseq_init(fp);
    unsigned long count=0,F[7]={0},R[7]={0},both[7]={0};
    FILE *bitmap=argc==6?fopen(argv[5],"wb"):NULL;
    if(argc==6 && !bitmap)return 1;
    while(kseq_read(ks)>=0){
        if(ks->seq.l!=150)return 1;
        char reverse[151];for(int i=0;i<150;++i)reverse[i]=nt4_seq_table[3-seq_nt4_table[(unsigned char)ks->seq.s[149-i]]];reverse[150]=0;
        unsigned f=present_path(h,ks->seq.s,k,mask,seeds,ns),r=present_path(h,reverse,k,mask,seeds,ns);
        if(bitmap)fputc(f|r,bitmap);
        for(int c=0;c<7;++c){F[c]+=(f>>c)&1;R[c]+=(r>>c)&1;both[c]+=((f|r)>>c)&1;}
        ++count;
    }
    printf("{\"diagnostic\":\"oracle graph-path existence; NOT reconstruction\",\"k\":%d,\"references\":%lu,\"bounds\":[",k,count);
    for(int c=0;c<7;++c)printf("%s{\"cutoff\":%d,\"F\":%lu,\"R\":%lu,\"either_direction\":%lu}",c?",":"",c,F[c],R[c],both[c]);
    puts("]}");
    if(bitmap)fclose(bitmap);
    kseq_destroy(ks);gzclose(fp);free(seeds);c4x_destroy(h);return 0;
}
