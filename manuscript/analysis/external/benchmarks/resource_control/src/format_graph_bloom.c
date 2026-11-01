/* Reference-free format adapters. Graph counting and bounded traversal are
   included unchanged; metadata layout, discovered anchors, and inner checks differ. */
#include "dbgps_strand.h"
static int mode=0, direction=0;
static int format_check(const dbgps_codec_t*,const dbgps_redundancy_t*,const dbgps_strand_layout_t*,const unsigned char*,unsigned char*,unsigned char*);
#define dbgps_strand_check format_check
#define main original_dbgps_main
#include "fountain_DBGPS2.c"
#undef main
#undef dbgps_strand_check
#include <sys/time.h>
#include <sys/resource.h>
static double cpu_seconds(void){struct rusage r;getrusage(RUSAGE_SELF,&r);return r.ru_utime.tv_sec+r.ru_utime.tv_usec/1e6+r.ru_stime.tv_sec+r.ru_stime.tv_usec/1e6;}
static double wall(void){struct timeval t;gettimeofday(&t,0);return t.tv_sec+t.tv_usec/1e6;}
static void rc(char*s){int n=strlen(s);for(int i=0;i<(n+1)/2;i++){char a=s[i],b=s[n-1-i];s[i]=seq_nt4_table[(unsigned char)b]<4?nt4_seq_table[3-seq_nt4_table[(unsigned char)b]]:'N';s[n-1-i]=seq_nt4_table[(unsigned char)a]<4?nt4_seq_table[3-seq_nt4_table[(unsigned char)a]]:'N';}}
static int cycle(char b){return b=='A'?0:b=='G'?1:b=='C'?2:3;}
static int nbt_valid(const char*s){unsigned v=0;for(int i=1;i<110;i++){unsigned d=(cycle(s[i])-cycle(s[i-1])+3)&3;if(d==3)return 0;if(i==109)return d==0;v=v*3+d;if(i%4==0){if(v>=64)return 0;v=0;}}return 0;}
static int format_check(const dbgps_codec_t*c,const dbgps_redundancy_t*r,const dbgps_strand_layout_t*l,const unsigned char*a,unsigned char*b,unsigned char*s){int v=dbgps_strand_check(c,r,l,a,b,s);if(v<0||mode!=1)return v;unsigned char dna[113];packed_dna_string(b,dna,112);char*body=(char*)dna+2;if(direction)rc(body);return nbt_valid(body)?v:-1;}
KHASHL_MAP_INIT(KH_LOCAL, seedmap_t, seedmap, uint64_t, unsigned, kh_hash_uint64, kh_eq_generic)
typedef struct { uint64_t key;unsigned count;int dir;} seed_t;
static int comp(const void*a,const void*b){const seed_t*x=a,*y=b;return x->dir!=y->dir?x->dir-y->dir:x->key<y->key?-1:x->key>y->key;}
static void add(seedmap_t*m,uint64_t key){int ret;khint_t z=seedmap_put(m,key,&ret);if(ret)kh_val(m,z)=0;kh_val(m,z)++;}
int main(int argc,char**argv){
 const char* wenv=getenv("BENCHMARK_COUNT_WORKERS"); int benchmark_workers=wenv?atoi(wenv):8; if(benchmark_workers!=1&&benchmark_workers!=3&&benchmark_workers!=8)return 2;
 fprintf(stderr,"RESOURCE count_workers=%d pipeline_threads=%d search_workers=1\n",benchmark_workers,benchmark_pipeline_threads());
 if(argc<4)return 1;mode=!strcmp(argv[1],"nbt")?1:!strcmp(argv[1],"fountain")?2:3;
 dbgps_codec_t mapping={0};const dbgps_codec_t*codec=NULL;char error[512];dbgps_redundancy_t red={0};
 if(mode==3){if(dbgps_codec_load(&mapping,argv[3],error,sizeof(error))<0)return 2;codec=&mapping;}
 if(dbgps_redundancy_init(&red,mode==1?"none":mode==2?"rs:2B":"crc16",error,sizeof(error))<0)return 2;
 if(!strcmp(argv[2],"validate")){
  FILE*in=fopen(argv[4],"r"),*out=fopen(argv[5],"w");char seq[1024];int len=mode==1?110:mode==2?152:164;
  dbgps_strand_layout_t l;dbgps_strand_layout(codec,&red,16,mode==1?96:mode==2?128:140,&l);
  while(fscanf(in,"%1023s",seq)==1){if(strlen(seq)!=len||strspn(seq,"ACGT")!=len)continue;
   int valid=0;if(mode==1)valid=nbt_valid(seq);else{unsigned char a[256]={0},b[256],s[256];for(int i=0;i<len;i++)a[i/4]|=seq_nt4_table[(unsigned char)seq[i]]<<(6-2*(i%4));valid=dbgps_strand_check(codec,&red,&l,a,b,s)==0;}
   if(valid)fprintf(out,"%s\n",seq);
  }fclose(in);fclose(out);return 0;
 }
 if(argc!=9){fprintf(stderr,"mode graph reads metadata-map k c r output\n");return 1;}
 const char*reads=argv[2];int k=atoi(argv[4]),cut=atoi(argv[5]);double ratio=atof(argv[6]);const char*outpath=argv[7];int screen=atoi(argv[8]);
 int query_bits=0,filter_shift=0;const char*qb=getenv("SPEED_QUERY_BITS"),*fb=getenv("SPEED_FILTER_SHIFT");
 if(qb&&dbgps_bloom_parse_bits(qb,&query_bits)<0)return 5;
 if(fb&&dbgps_bloom_parse_shift(fb,&filter_shift)<0)return 5;
 double start=wall(),count_cpu_start=cpu_seconds();char*input_file=(char*)reads;
 kc_c4x_t*h=filter_shift?count_files_filtered(1,&input_file,k,KC_BITS,10000000,benchmark_workers,200,filter_shift):count_file(reads,k,KC_BITS,10000000,benchmark_workers,200);
 if(!h)return 3;fprintf(stderr,"COUNT wall=%.6f cpu=%.6f\n",wall()-start,cpu_seconds()-count_cpu_start);
 double buildstart=wall(),bloom_cpu_start=cpu_seconds();if(c4x_set_bloom(h,query_bits)<0)return 5;
 fprintf(stderr,"BLOOM bits=%d bytes=%zu filter_shift=%d wall=%.6f cpu=%.6f\n",query_bits,h->bloom?h->bloom->n_bytes:0,filter_shift,wall()-buildstart,cpu_seconds()-bloom_cpu_start);
 seed_t*seeds=NULL;size_t ns=0;double seedstart=wall(),seed_cpu_start=cpu_seconds();
 if(mode==3){ns=101295684ULL-101010102ULL+1;seeds=calloc(ns,sizeof(*seeds));uint64_t anchor=actgkmer_uint64((unsigned char*)"GCAGAGTAGCATGTC",15);for(size_t i=0;i<ns;i++){uint64_t idx;dbgps_codec_encode_index(codec,101010102+i,16,&idx);seeds[i].key=(anchor<<32)|idx;seeds[i].count=1;}}
 else{
  seedmap_t*maps[2]={seedmap_init(),seedmap_init()};gzFile f=gzopen(reads,"r");kseq_t*q=kseq_init(f);unsigned long nr=0;
  const char*anchors[2]={"GCAACAAGTCAATCCGT","CAAGCAAGCATTCAATT"};
  while(kseq_read(q)>=0){nr++;char seq[2048];if(q->seq.l>=sizeof(seq))return 4;strcpy(seq,q->seq.s);
   if(mode==2){char*p=strstr(seq,"CCTGCAGAGT");if(p&&strlen(p)>=26&&strspn(p+10,"ACGT")>=16){uint64_t key=actgkmer_uint64((unsigned char*)p+26-k,k);add(maps[0],key);}}
   else for(int ori=0;ori<2;ori++){if(ori)rc(seq);for(int d=0;d<2;d++){char*p=strstr(seq,anchors[d]);if(p&&strlen(p)>=31&&strspn(p+17,"ACGT")>=14)add(maps[d],actgkmer_uint64((unsigned char*)p+31-k,k));}}
  }
  kseq_destroy(q);gzclose(f);ns=kh_size(maps[0])+kh_size(maps[1]);seeds=calloc(ns,sizeof(*seeds));size_t z=0;
  for(int d=0;d<2;d++){for(khint_t x=0;x!=kh_end(maps[d]);x++)if(kh_exist(maps[d],x)){seeds[z++]=(seed_t){kh_key(maps[d],x),kh_val(maps[d],x),d};}seedmap_destroy(maps[d]);}qsort(seeds,ns,sizeof(*seeds),comp);
 }
 fprintf(stderr,"SEED n=%zu wall=%.6f cpu=%.6f\n",ns,wall()-seedstart,cpu_seconds()-seed_cpu_start);double searchstart=wall(),search_cpu_start=cpu_seconds();FILE*out=fopen(outpath,"w");uint64_t mask=(1ULL<<(2*k))-1;unsigned long returned=0,queries=0;
 for(size_t i=0;i<ns;i++){if(screen&&(kh_hash_uint64(seeds[i].key)%128))continue;unsigned char*paths[Max_Path_Num];direction=seeds[i].dir;int pos=find_paths_protected(h,paths,seeds[i].key,mask,cut,mode==1?96:mode==2?128:140,k,16,ratio,0,codec,&red);queries++;
  if(pos>=0){char dna[512];packed_dna_string(paths[pos],(unsigned char*)dna,mode==1?112:mode==2?152:164);char*body=dna+(mode==1?2:0);if(mode==1&&direction)rc(body);fprintf(out,"%s\n",body);free(paths[pos]);returned++;}
 }
 fclose(out);fprintf(stderr,"SEARCH queries=%lu returned=%lu wall=%.6f cpu=%.6f\n",queries,returned,wall()-searchstart,cpu_seconds()-search_cpu_start);c4x_destroy(h);free(seeds);dbgps_codec_destroy(&mapping);return 0;
}
