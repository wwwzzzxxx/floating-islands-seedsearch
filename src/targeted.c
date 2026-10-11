/* targeted.c — 只对指定的 (小屋基, 神殿偏移) 配对跑高16位变体搜索
 *
 * 用途：pairs.c 已经算出每个配对的"单人覆盖率"，我们只关心覆盖率高的那些，
 *       没必要把 21979 个配对全跑一遍。本工具读入配对清单，只对这些跑完整筛选。
 *
 * 用法: targeted <monpairs_sorted.txt> [取前N个] [-t 线程数] [-out 文件]
 * 输出格式与 stage3 一致，额外附上神殿覆盖率。
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "generator.h"
#include "finders.h"
#include "quadbase.h"
#include "biomes.h"
#define MCV MC_1_18_2

static int P_VIL_D = 256, P_RIVER_SPAN = 3000;
static StructureConfig hutconf, monconf;
static int biome_river;

typedef struct { uint64_t base; int oi, oj, cov; } Pair;

typedef struct {
    Pair *pairs; int npairs;
    int next;
    pthread_mutex_t mtx, rmtx;
    int nhits; FILE *out;
    uint64_t n_var, n_hutbio, n_monbio, n_vil, n_river;
} Ctx;

static inline int d2(int x1,int z1,int x2,int z2){int dx=x1-x2,dz=z1-z2;return dx*dx+dz*dz;}

/* 与 stage3.c 完全相同的河流连通块计算 */
static int river_blob(Generator *g, int span, double *frac, double *cx, double *cz)
{
    int s = span/4; if (s<1) s=1;
    Range r = {4, -s/2, -s/2, s, s};
    int *cache = malloc(sizeof(int)*(size_t)s*s);
    *frac=-1;
    if (!cache) return -1;
    if (genBiomes(g, cache, r)!=0){ free(cache); return -1; }
    int total=s*s, cnt=0, i;
    for (i=0;i<total;i++) if (cache[i]==biome_river) cnt++;
    *frac=(double)cnt/total;
    char *vis=calloc(total,1); int *st=malloc(sizeof(int)*total);
    int best=0; double bx=0,bz=0;
    if (vis&&st){
        for (i=0;i<total;i++){
            if (cache[i]!=biome_river||vis[i]) continue;
            int sp=0,size=0; double sx=0,sz=0;
            st[sp++]=i; vis[i]=1;
            while(sp>0){ int c=st[--sp]; int x=c%s,z=c/s; size++; sx+=x; sz+=z;
                int k,nx,nz;
                for(k=0;k<4;k++){ nx=x+(k==0?1:k==1?-1:0); nz=z+(k==2?1:k==3?-1:0);
                    if(nx<0||nz<0||nx>=s||nz>=s) continue;
                    int ni=nz*s+nx; if(!vis[ni]&&cache[ni]==biome_river){vis[ni]=1;st[sp++]=ni;} } }
            if (size>best){ best=size; bx=sx/size; bz=sz/size; }
        }
    }
    free(vis); free(st); free(cache);
    *cx=(bx-s/2)*4.0; *cz=(bz-s/2)*4.0;
    return best;
}

static void *worker(void *arg)
{
    Ctx *c = (Ctx*)arg;
    Generator g; setupGenerator(&g, MCV, 0);
    static const int rx[4]={0,1,0,1}, rz[4]={0,0,1,1};
    for (;;) {
        pthread_mutex_lock(&c->mtx);
        int idx = c->next++;
        pthread_mutex_unlock(&c->mtx);
        if (idx >= c->npairs) break;
        Pair pr = c->pairs[idx];
        uint64_t base = pr.base & 0xFFFFFFFFFFFFULL;
        int oi = pr.oi, oj = pr.oj;

        /* base 级常量：小屋半径 + 4 个小屋坐标（都只取决于低 48 位） */
        float rh = isQuadBase(hutconf, base, 128);
        if (rh == 0.0f) continue;
        Pos hp[4];
        for (int q=0;q<4;q++) getStructurePos(Swamp_Hut, MCV, base, rx[q], rz[q], &hp[q]);

        for (uint64_t u = 0; u < 65536; u++) {
            uint64_t seed = base | (u << 48);
            c->n_var++;
            applySeed(&g, DIM_OVERWORLD, seed);
            int i, ok = 1;
            for (i=0;i<4;i++){
                if (!isViableStructurePos(Swamp_Hut, &g, hp[i].x, hp[i].z, 0)) { ok=0; break; }
            }
            if (!ok) continue;
            c->n_hutbio++;
            Pos mp[4];
            for (i=0;i<4;i++){
                getStructurePos(Monument, MCV, seed, oi+rx[i], oj+rz[i], &mp[i]);
                if (!isViableStructurePos(Monument, &g, mp[i].x, mp[i].z, 0)) { ok=0; break; }
            }
            if (!ok) continue;
            c->n_monbio++;
            int bd=0x7fffffff, found=0, gx, gz;
            for (gz=-1;gz<=1;gz++) for (gx=-1;gx<=1;gx++){
                Pos vp;
                if(!getStructurePos(Village,MCV,seed,gx,gz,&vp)) continue;
                int dd=d2(vp.x,vp.z,0,0);
                if(dd>=bd) continue;
                if(!isViableStructurePos(Village,&g,vp.x,vp.z,0)) continue;
                bd=dd; found=1;
            }
            if (!found || bd > P_VIL_D*P_VIL_D) continue;
            c->n_vil++;
            double frac, rcx, rcz;
            int blob = river_blob(&g, P_RIVER_SPAN, &frac, &rcx, &rcz);
            if (blob <= 0) continue;
            c->n_river++;
            pthread_mutex_lock(&c->rmtx);
            fprintf(c->out,
                "%lld hut_r=%.1f mon_cov=%d vil_d=%d river%%=%.3f blob=%d river_c=(%.0f,%.0f) mon_off=(%d,%d)\n",
                (long long)(int64_t)seed, rh, pr.cov, (int)sqrt((double)bd),
                frac*100, blob, rcx, rcz, oi, oj);
            c->nhits++;
            pthread_mutex_unlock(&c->rmtx);
        }
    }
    return NULL;
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr,"用法: %s <配对文件> [取前N个] [-t N] [-out f]\n", argv[0]); return 1; }
    const char *pairfile = argv[1];
    int TOPN = 200, NT = 24;
    const char *outpath = "targeted_out.txt";
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i],"-t") && i+1<argc) NT = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-out") && i+1<argc) outpath = argv[++i];
        else if (argv[i][0] != '-') TOPN = atoi(argv[i]);
    }

    Pair *all = malloc(sizeof(Pair)*30000); int nall = 0;
    FILE *f = fopen(pairfile, "r");
    if (!f) { perror(pairfile); return 1; }
    while (nall < 30000) {
        unsigned long long b; int oi, oj, cov;
        if (fscanf(f, "%llu %d %d %d", &b, &oi, &oj, &cov) != 4) break;
        all[nall].base = b; all[nall].oi = oi; all[nall].oj = oj; all[nall].cov = cov;
        nall++;
    }
    fclose(f);
    /* 按覆盖率降序（文件若已排序则不变） */
    for (int i=0;i<nall;i++) for (int j=i+1;j<nall;j++)
        if (all[j].cov > all[i].cov) { Pair t=all[i]; all[i]=all[j]; all[j]=t; }
    if (TOPN > nall) TOPN = nall;
    printf("配对文件 %d 个，取覆盖率前 %d 个（覆盖率 %d ~ %d）\n",
           nall, TOPN, all[0].cov, all[TOPN-1].cov);

    getStructureConfig(Swamp_Hut, MCV, &hutconf);
    getStructureConfig(Monument,  MCV, &monconf);
    { Generator g; setupGenerator(&g, MCV, 0); applySeed(&g, DIM_OVERWORLD, 0);
      biome_river = river;   /* 与 stage3.c 一致：直接用 cubiomes 的 river 群系枚举 */
      printf("river 群系 id = %d\n", biome_river);
    }

    Ctx ctx; memset(&ctx, 0, sizeof ctx);
    ctx.pairs = all; ctx.npairs = TOPN;
    ctx.out = fopen(outpath, "w");
    if (!ctx.out) { perror(outpath); return 1; }
    pthread_mutex_init(&ctx.mtx, NULL); pthread_mutex_init(&ctx.rmtx, NULL);

    pthread_t th[256];
    for (int i = 0; i < NT; i++) pthread_create(&th[i], NULL, worker, &ctx);
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    fclose(ctx.out);

    printf("[统计] 变体检查 %llu | 4小屋swamp %llu | +4神殿ocean %llu | +村庄 %llu | +河流 %llu | 合规 %d\n",
           (unsigned long long)ctx.n_var, (unsigned long long)ctx.n_hutbio,
           (unsigned long long)ctx.n_monbio, (unsigned long long)ctx.n_vil,
           (unsigned long long)ctx.n_river, ctx.nhits);
    return 0;
}
