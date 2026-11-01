/* Experimental primer-free, internal-index adapter. Reuses DBGPS2's exact
 * counter, canonical graph, and dynamic coverage rule. This is not stock CLI.
 * No design/reference DNA is read until the optional post-assembly diagnostic.
 */
#define main dbgps2_original_main
#include "../../DBGPS2.c"
#undef main
#include <sys/time.h>

#define LIMIT 1050
#define INDEX_XOR 1071068U
static int terminal_layout = 0;
static int prefix_layout = 0;
typedef struct { uint64_t key; unsigned id; } seed_t;
typedef struct { uint64_t left, right; unsigned lc, rc; char dna[61]; } state_t;
static double walltime(void) { struct timeval t; gettimeofday(&t, 0); return t.tv_sec+t.tv_usec/1e6; }
static unsigned mul64(unsigned a, unsigned b) {
    unsigned r=0;
    while (b) { if(b&1) r^=a; b>>=1; a<<=1; if(a&64) a^=91; }
    return r;
}
static int rs_valid(const char *s) {
    unsigned s0=0,s1=0,a=1;
    for(int j=0;j<20;++j) {
        unsigned c=0;
        if(terminal_layout) {
            for(int x=0;x<3;++x)c|=seq_nt4_table[(unsigned char)s[3*j+x]]<<(2*x);
        } else for(int x=57-3*j;x<60-3*j;++x) c=(c<<2)|seq_nt4_table[(unsigned char)s[x]];
        s0^=c; s1^=mul64(c,a); a=mul64(a,2);
    }
    return s0==0 && s1==0;
}
static int compare_seed(const void *a, const void *b) {
    const seed_t *x=a,*y=b;
    return x->id<y->id?-1:x->id>y->id?1:x->key<y->key?-1:x->key>y->key;
}
static int allowed(unsigned prev,unsigned next,int cut,double ratio) {
    int lo=prev/ratio,hi=fmin(KC_MAX,prev*ratio);
    if(lo<cut+1) lo=cut+1;
    return next>=(unsigned)lo && next<=(unsigned)hi;
}
static unsigned cov(kc_c4x_t *h,uint64_t key,int k,uint64_t mask) {
    return kmer_cov(min_hash_key(key,k),mask,h);
}
static void run(kc_c4x_t *h,seed_t *seeds,size_t nseed,int n,int k,int cut,
                double ratio,const char *prefix) {
    char fname[1024];
    snprintf(fname,sizeof(fname),"%s.k%d.c%d.r%.2f.strands",prefix,k,cut,ratio);
    FILE *out=fopen(fname,"w");
    snprintf(fname,sizeof(fname),"%s.k%d.c%d.r%.2f.status.tsv",prefix,k,cut,ratio);
    FILE *status=fopen(fname,"w");
    fprintf(status,"id\tseeds\tcomplete_paths\tvalid_paths\tstatus\n");
    uint64_t mask=(1ULL<<(2*k))-1;
    state_t *a=calloc(LIMIT,sizeof(*a)),*b=calloc(LIMIT,sizeof(*b));
    size_t ptr=0; int recovered=0,overflows=0,ambiguous=0; double start=walltime();
    int anchor_end=prefix_layout?k:terminal_layout?60:54;
    for(int id=0;id<n;++id) {
        while(ptr<nseed && seeds[ptr].id<(unsigned)id) ++ptr;
        size_t lo=ptr;while(ptr<nseed && seeds[ptr].id==(unsigned)id) ++ptr;
        int count=0,overflow=0,valid=0,distinct=0;char winner[61]={0};
        for(size_t z=lo;z<ptr;++z) {
            unsigned c=cov(h,seeds[z].key,k,mask);
            if(c<=(unsigned)cut)continue;
            if(count==LIMIT) {overflow=1;break;}
            state_t *v=&a[count++];memset(v,0,sizeof(*v));
            v->left=v->right=seeds[z].key;v->lc=v->rc=c;
            uint64_t key=v->left;
            for(int j=anchor_end-1;j>=anchor_end-k;--j) {v->dna[j]=nt4_seq_table[key&3];key>>=2;}
        }
        for(int step=0;step<60-k && count && !overflow;++step) {
            int pos=step<anchor_end-k ? anchor_end-k-1-step : k+step;
            int left=step<anchor_end-k,newcount=0;
            for(int j=0;j<count && !overflow;++j)for(int c=0;c<4;++c) {
                uint64_t key=left ? (a[j].left>>2)|((uint64_t)c<<(2*(k-1))) : ((a[j].right<<2)|c)&mask;
                unsigned coverage=cov(h,key,k,mask);
                if(!allowed(left?a[j].lc:a[j].rc,coverage,cut,ratio))continue;
                if(newcount==LIMIT) {overflow=1;break;}
                state_t *v=&b[newcount++];*v=a[j];v->dna[pos]=nt4_seq_table[c];
                if(left){v->left=key;v->lc=coverage;}else{v->right=key;v->rc=coverage;}
            }
            state_t *tmp=a;a=b;b=tmp;count=newcount;
        }
        if(!overflow)for(int j=0;j<count;++j)if(rs_valid(a[j].dna)) {
            ++valid;
            if(!distinct){memcpy(winner,a[j].dna,61);distinct=1;}
            else if(strcmp(winner,a[j].dna))distinct=2;
        }
        const char *why=overflow?"overflow":distinct==1?"unique":distinct>1?"ambiguous":"no_valid_path";
        if(overflow)++overflows;
        else if(distinct==1){fprintf(out,"%d\t%s\n",id,winner);++recovered;}
        else if(distinct>1)++ambiguous;
        fprintf(status,"%d\t%zu\t%d\t%d\t%s\n",id,ptr-lo,count,valid,why);
    }
    fclose(out);fclose(status);free(a);free(b);
    fprintf(stderr,"RESULT k=%d cutoff=%d ratio=%.2f recovered=%d ambiguous=%d overflow=%d wall=%.3f\n",k,cut,ratio,recovered,ambiguous,overflows,walltime()-start);
}
int main(int argc,char **argv) {
    if(argc<6){fprintf(stderr,"usage: adapter reads n k output-prefix truth-fasta-or-- [ratiosCSV [cutoffsCSV]]\n");return 1;}
    int n=atoi(argv[2]),k=atoi(argv[3]);
    const char *style=getenv("ANTK_LAYOUT");
    terminal_layout=style && !strcmp(style,"terminal");
    prefix_layout=style && !strcmp(style,"prefix");
    if(k<12||k>31)return 1;
    double start=walltime();uint64_t mask=(1ULL<<(2*k))-1;
    const char *filter=getenv("ANTK_BLOOM");
    int filter_shift=filter && atoi(filter)>1?atoi(filter):31;
    kc_c4x_t *h=filter && atoi(filter) ? count_files_filtered(1,argv+1,k,KC_BITS,10000000,8,60,filter_shift)
                                     : count_file(argv[1],k,KC_BITS,10000000,8,60);
    if(!h)return 2;
    fprintf(stderr,"COUNT k=%d wall=%.3f\n",k,walltime()-start);
    size_t ns=0,cap=65536;seed_t *seeds=malloc(cap*sizeof(*seeds));
    for(int p=0;p<1<<KC_BITS;++p) {
        kc_c4_t *t=h->h[p];
        for(khint_t j=0;j<kh_end(t);++j)if(kh_exist(t,j) && (kh_key(t,j)&KC_MAX)>1) {
            uint64_t hash=((kh_key(t,j)>>KC_BITS)<<KC_BITS)|p;
            uint64_t key=hash64i(hash,mask);
            for(int ori=0;ori<2;++ori) {
                unsigned id=(key&0xffffffU)^INDEX_XOR;
                if(terminal_layout) {
                    uint64_t v=key;id=0;
                    for(int b=0;b<12;++b){id=(id<<2)|(v&3);v>>=2;}
                }
                if(prefix_layout)id=(key>>(2*(k-12)))&0xffffffU;
                if(id<(unsigned)n) {
                    if(ns==cap){cap*=2;seeds=realloc(seeds,cap*sizeof(*seeds));}
                    seeds[ns++]=(seed_t){key,id};
                }
                key=comp_rev(key,k);
            }
        }
    }
    qsort(seeds,ns,sizeof(*seeds),compare_seed);
    fprintf(stderr,"SEEDS %zu wall_since_start=%.6f\n",ns,walltime()-start);
    if(argc>6) {
        char *ratios=strdup(argv[6]),*rsave=NULL;
        for(char *r=strtok_r(ratios,",",&rsave);r;r=strtok_r(NULL,",",&rsave)) {
            double ratio=atof(r);
            if(!isfinite(ratio)||ratio<1) {fprintf(stderr,"invalid ratio: %s\n",r);return 1;}
            char *cuts=strdup(argc>7?argv[7]:"1"),*csave=NULL;
            for(char *c=strtok_r(cuts,",",&csave);c;c=strtok_r(NULL,",",&csave))
                run(h,seeds,ns,n,k,atoi(c),ratio,argv[4]);
            free(cuts);
        }
        free(ratios);
    }
    else for(int c=1;c<=3;c+=2)for(int r=0;r<2;++r)run(h,seeds,ns,n,k,c,r?10:2.55,argv[4]);
    /* Reference-only diagnostics are evaluated strictly after reconstruction. */
    if(strcmp(argv[5],"-")) {
        char path[1024];snprintf(path,sizeof(path),"%s.k%d.truth_support.tsv",argv[4],k);
        FILE *out=fopen(path,"w");fprintf(out,"id\tmin_kmer_count\tmax_adjacent_ratio\tstrict_rs_valid\n");
        gzFile f=gzopen(argv[5],"r");kseq_t *ks=kseq_init(f);int id=0;
        while(kseq_read(ks)>=0) {
            unsigned prev=0,minc=KC_MAX;double maxr=1;
            for(int pos=0;pos<=60-k;++pos) {
                uint64_t key=actgkmer_uint64((unsigned char*)ks->seq.s+pos,k);
                unsigned c=cov(h,key,k,mask);if(c<minc)minc=c;
                if(pos && c && prev){double r=fmax((double)c/prev,(double)prev/c);if(r>maxr)maxr=r;}
                prev=c;
            }
            fprintf(out,"%d\t%u\t%.5f\t%d\n",id++,minc,maxr,rs_valid(ks->seq.s));
        }
        kseq_destroy(ks);gzclose(f);fclose(out);
    }
    free(seeds);c4x_destroy(h);return 0;
}
