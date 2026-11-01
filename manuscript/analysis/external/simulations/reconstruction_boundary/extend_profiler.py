"""Generate an oracle-diagnostic profiler without changing production search code."""
from pathlib import Path
R=Path(__file__).resolve().parent
s=(R/'src/original_profile.c').read_text()
s=s.replace('#include <sys/time.h>','#include <assert.h>\n#include <sys/time.h>')
s=s.replace('int do_stats = 1;','int do_stats = 1, validate_every = 0;')
s=s.replace('F:NT:", 0)','F:NT:v:", 0)')
s=s.replace("else if (c == 'T') trace_fn = o.arg;","else if (c == 'T') trace_fn = o.arg;\n        else if (c == 'v') validate_every = atoi(o.arg);")
marker='    /* ---- coverage spectrum: reference vs noise k-mers ---- */'
insert=r'''
    /* Oracle diagnostics are never used for seed enumeration or path decisions. */
    const int nstep = g.path_len + 1;
    uint16_t *truth = calloc((size_t)g.n_id * (nstep + 1), sizeof(uint16_t));
    uint16_t *alts = calloc((size_t)g.n_id * nstep * 3, sizeof(uint16_t));
    assert(truth && alts);
    double mean_cov_all = 0, cv_all = 0, mincov_all = 0;
    uint64_t nref_profile = 0;
    int start_pos = pr_len + index_len - k;
    for (uint64_t i = 0; i < g.n_id; ++i) {
        if (!ref_seq[i]) continue;
        uint16_t *v = truth + (size_t)i * (nstep + 1);
        uint64_t parent = kmer_of(ref_seq[i] + start_pos, k);
        v[0] = kmer_cov(min_hash_key(parent,k),g.mask,h);
        double su=0, sq=0; int minv=KC_MAX;
        for (int p=0;p<nstep;++p) {
            uint64_t child=kmer_of(ref_seq[i]+start_pos+p+1,k);
            v[p+1]=kmer_cov(min_hash_key(child,k),g.mask,h);
            if(p<nstep-1){su+=v[p+1]; sq+=(double)v[p+1]*v[p+1];
                if(v[p+1]<minv)minv=v[p+1];}
            int j=0;
            for(int a=0;a<4;++a){
                uint64_t alt=((parent<<2)|a)&g.mask;
                if(alt!=child)alts[((size_t)i*nstep+p)*3+j++]=kmer_cov(min_hash_key(alt,k),g.mask,h);
            }
            parent=child;
        }
        double m=su/(nstep-1), var=sq/(nstep-1)-m*m;
        mean_cov_all+=m; cv_all+=m>0?sqrt(var>0?var:0)/m:0; mincov_all+=minv;
        ++nref_profile;
    }
    /* Compact, lossless per-step evidence, little-endian uint16 on this host. */
    snprintf(fn,sizeof fn,"%s.truth.u16",prefix);fo=fopen(fn,"wb");
    fwrite(truth,sizeof(uint16_t),(size_t)g.n_id*(nstep+1),fo);fclose(fo);
    snprintf(fn,sizeof fn,"%s.alts.u16",prefix);fo=fopen(fn,"wb");
    fwrite(alts,sizeof(uint16_t),(size_t)g.n_id*nstep*3,fo);fclose(fo);
    snprintf(fn,sizeof fn,"%s.oracle.json",prefix);fo=fopen(fn,"w");
    fprintf(fo,"{\"min_id\":%" PRIu64 ",\"n_id\":%" PRIu64 ",\"steps\":%d,\"dtype\":\"little-endian uint16\",\"k\":%d}\n",min_id,g.n_id,nstep,k);fclose(fo);
'''
assert marker in s;s=s.replace(marker,insert+'\n'+marker)
s=s.replace('g.r = rs[ri] <= 0 ? 1e9 : rs[ri];','g.r = rs[ri] <= 0 ? 1e6 : rs[ri];')
s=s.replace('mean_width_ref\\tseconds\\n','mean_width_ref\\tseconds\\tceil_cov\\tceil_end\\tceil_seed\\tceil_adm\\tceil_cont\\tlookahead_losses\\tseed_losses\\tratio_losses\\tsearch_gap\\toverflow_adm\\tambiguous_adm\\twrong_adm\\tno_valid_adm\\tnoise_per_admitted_step\\tmean_true_cov\\tmean_cov_cv\\tmean_mincov\\tnative_checks\\n')
marker='        fprintf(fg, "%d\\t%g\\t%d\\t%" PRIu64'
insert=r'''
        uint64_t uc=0,ue=0,us=0,ua=0,cont=0,oa=0,aa=0,wa=0,va=0,noise=0,nchecked=0;
        for(uint64_t i=0;i<g.n_id;++i){
            if(!g.ref_packed[i])continue;
            uint16_t *v=truth+(size_t)i*(nstep+1);
            int covok=1,endok=0,adm=1,con=1; uint64_t nb=0;
            for(int p=0;p<nstep;++p){
                int a=v[p], b=v[p+1];
                int lo=(int)(a/g.r), hi=a*g.r>KC_MAX?KC_MAX:(int)(a*g.r);
                if(lo<g.cut+1)lo=g.cut+1;
                if(p<nstep-1){
                    if(b<g.cut+1)covok=0;
                    if(b<g.cut+1 || (p>=g.skip && (b<lo || b>hi)))adm=0;
                    if(b<g.cut+1 || (p>=g.skip && (!a || !b || (double)b>a*g.r || (double)b<a/g.r)))con=0;
                    for(int j=0;j<3;++j){int ac=alts[((size_t)i*nstep+p)*3+j];
                        if(p<g.skip ? ac>=g.cut+1 : ac>=lo && ac<=hi)++nb;
                    }
                }else{
                    /* Any terminal successor completes the same reported CRC codeword. */
                    int last_adm=0,last_cont=0;
                    for(int j=-1;j<3;++j){int ac=j<0?b:alts[((size_t)i*nstep+p)*3+j];
                        if(ac>=g.cut+1){
                            endok=1;
                            if(p<g.skip || (ac>=lo && ac<=hi))last_adm=1;
                            if(p<g.skip || (a && (double)ac<=a*g.r && (double)ac>=a/g.r))last_cont=1;
                        }
                    }
                    if(!last_adm)adm=0;if(!last_cont)con=0;
                }
            }
            endok=endok&&covok;int seedok=endok&&(v[0]>0 || g.skip>0);
            uc+=covok;ue+=endok;us+=seedok;ua+=adm;cont+=con;
            if(adm){noise+=nb;oa+=g.outcome[i]==OUT_OVERFLOW;aa+=g.outcome[i]==OUT_AMBIG;wa+=g.outcome[i]==OUT_WRONG;va+=g.outcome[i]==OUT_NOVALID;}
            assert(!adm || seedok);assert(!con || adm);
            assert(g.outcome[i]!=OUT_OK || adm);
            if(validate_every>0 && i%(uint64_t)validate_every==0){
                uint64_t eid; if(codec)dbgps_codec_encode_index(codec,min_id+i,index_len,&eid);else eid=index_encode(min_id+i,index_len);
                unsigned char out[(Max_Path_Len+3)/4];uint32_t nn;uint16_t width;
                int inst=search_one(&g,&g.bufs[0],g.pr|eid,out,&nn,&width);
                unsigned char *paths[Max_Path_Num]={0};
                int native=find_paths_protected(h,paths,g.pr|eid,g.mask,g.cut,payload_len,k,index_len,g.r,skip,codec,&red);
                assert((native>=0)==(inst==OUT_OK));
                if(native>=0){assert(!memcmp(paths[native],out,g.seq_len));free(paths[native]);}
                ++nchecked;
            }
        }
        assert(cnt[OUT_OK]<=ua && cont<=ua && ua<=us && us<=ue && ue<=uc);
'''
assert marker in s;s=s.replace(marker,insert+marker)
s=s.replace('"\\t%.1f\\t%.2f\\t%.2f\\n",','"\\t%.1f\\t%.2f\\t%.2f",')
marker='        fflush(fg);'
insert=r'''
        fprintf(fg,"\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%.9g\t%.9g\t%.9g\t%.9g\t%" PRIu64 "\n",
            uc,ue,us,ua,cont,uc-ue,ue-us,us-ua,ua-cnt[OUT_OK],oa,aa,wa,va,ua?(double)noise/ua/(nstep-1):0,
            mean_cov_all/nref_profile,cv_all/nref_profile,mincov_all/nref_profile,nchecked);
'''
assert marker in s;s=s.replace(marker,insert+marker)
(R/'src/boundary_profile.c').write_text(s)
print('Generated boundary_profile.c')
