/* 第三阶段（修正版）：
 *   小屋四联：锚定原点区域 (0,0)-(1,1)（靠近出生点）
 *   神殿四联：允许出现在 (±OFFSET 区域) 范围内——避免与小屋争同一片群系
 *   然后对每个 (base, 神殿偏移) 做 65536 个高16位变体搜索：
 *     4小屋全 swamp + 4神殿全 ocean + 出生点村庄<=256 + 3000x3000 河流连通块
 * 多线程（每线程自己的 Generator）。
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
static int      P_MON_R   = 256;     /* 神殿四联半径 */
static int      P_VIL_D   = 256;     /* 村庄距原点上限 */
static int      P_RIVER_SPAN = 3000;
static int      P_OFFSET  = 6;       /* 神殿四联搜索的区域偏移范围 ±N */
static int      P_THREADS = 32;
static StructureConfig hutconf, monconf;
static int biome_river;

typedef struct {
    uint64_t *bases; int nbases;
    int next;                     /* 下一个待处理 base 索引 */
    pthread_mutex_t mtx;
    /* 结果 */
    pthread_mutex_t rmtx;
    int nhits;
    FILE *out;
    uint64_t n_var, n_hutbio, n_monbio, n_vil, n_river, n_pairs;
} Ctx;

static inline int d2(int x1,int z1,int x2,int z2){int dx=x1-x2,dz=z1-z2;return dx*dx+dz*dz;}

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
    for (;;) {
        pthread_mutex_lock(&c->mtx);
        int idx = c->next++;
        pthread_mutex_unlock(&c->mtx);
        if (idx >= c->nbases) break;
        uint64_t base = c->bases[idx] & 0xFFFFFFFFFFFFULL;

        /* 找该 base 上所有偏移位置的神殿四联 */
        int offs[64][2], no = 0;
        for (int i = -P_OFFSET; i <= P_OFFSET && no < 64; i++)
            for (int j = -P_OFFSET; j <= P_OFFSET && no < 64; j++) {
                uint64_t sp = moveStructure(base, -i, -j);
                if (isQuadBase(monconf, sp, P_MON_R) != 0.0f) { offs[no][0]=i; offs[no][1]=j; no++; }
            }
        if (no == 0) continue;
        pthread_mutex_lock(&c->rmtx); c->n_pairs += no; pthread_mutex_unlock(&c->rmtx);

        /* 对每个偏移做变体搜索 */
        for (int k = 0; k < no; k++) {
            int oi = offs[k][0], oj = offs[k][1];
            for (uint64_t u = 0; u < 65536; u++) {
                uint64_t seed = base | (u << 48);
                c->n_var++;
                float rh = isQuadBase(hutconf, seed, 128);
                if (rh == 0) continue;
                uint64_t msp = moveStructure(seed, -oi, -oj);
                float rm = isQuadBase(monconf, msp, P_MON_R);
                if (rm == 0) continue;
                applySeed(&g, DIM_OVERWORLD, seed);
                Pos hp[4], mp[4];
                static const int rx[4]={0,1,0,1}, rz[4]={0,0,1,1};
                int i, ok = 1;
                for (i=0;i<4;i++){
                    getStructurePos(Swamp_Hut, MCV, seed, rx[i], rz[i], &hp[i]);
                    if (!isViableStructurePos(Swamp_Hut, &g, hp[i].x, hp[i].z, 0)) { ok=0; break; }
                }
                if (!ok) continue;
                c->n_hutbio++;
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
                    "%lld hut_r=%.1f mon_r=%.1f mon_off=(%d,%d) vil_d=%d river%%=%.3f blob=%d river_c=(%.0f,%.0f)\n",
                    (long long)(int64_t)seed, rh, rm, oi, oj, (int)sqrt((double)bd),
                    frac*100, blob, rcx, rcz);
                c->nhits++;
                pthread_mutex_unlock(&c->rmtx);
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr,"用法: %s <hutbases.txt> [-monr R] [-offset N] [-t N] [-out f]\n", argv[0]); return 1; }
    const char *outpath = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i],"-monr") && i+1<argc) P_MON_R = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-vild") && i+1<argc) P_VIL_D = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-riverspan") && i+1<argc) P_RIVER_SPAN = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-offset") && i+1<argc) P_OFFSET = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-t") && i+1<argc) P_THREADS = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-out") && i+1<argc) outpath = argv[++i];
    }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("打开 base 列表失败"); return 1; }
    uint64_t *bases = NULL; int nb = 0, cap = 0; char line[64];
    while (fgets(line, sizeof line, f)) {
        uint64_t v = strtoull(line, 0, 10);
        if (nb >= cap) { cap = cap?cap*2:1024; bases = realloc(bases, cap*sizeof(uint64_t)); }
        bases[nb++] = v;
    }
    fclose(f);
    getStructureConfig(Swamp_Hut, MCV, &hutconf);
    getStructureConfig(Monument,  MCV, &monconf);
    { Generator g; setupGenerator(&g, MCV, 0); applySeed(&g, DIM_OVERWORLD, 0); biome_river = river; }

    Ctx ctx; memset(&ctx, 0, sizeof ctx);
    ctx.bases = bases; ctx.nbases = nb; ctx.next = 0;
    pthread_mutex_init(&ctx.mtx, NULL); pthread_mutex_init(&ctx.rmtx, NULL);
    ctx.out = outpath ? fopen(outpath, "w") : stdout;
    if (!ctx.out) { perror("打开输出失败"); return 1; }

    fprintf(stderr, "base 数 %d，神殿半径 r=%d，偏移 ±%d 区域，线程 %d\n", nb, P_MON_R, P_OFFSET, P_THREADS);
    pthread_t th[512]; int nt = P_THREADS > 512 ? 512 : P_THREADS;
    for (int i = 0; i < nt; i++) pthread_create(&th[i], NULL, worker, &ctx);
    for (int i = 0; i < nt; i++) pthread_join(th[i], NULL);

    printf("\n[统计]\n");
    printf("(base,神殿偏移) 对 : %llu\n", (unsigned long long)ctx.n_pairs);
    printf("变体检查总数       : %llu\n", (unsigned long long)ctx.n_var);
    printf("4小屋全swamp       : %llu\n", (unsigned long long)ctx.n_hutbio);
    printf("+4神殿全ocean      : %llu\n", (unsigned long long)ctx.n_monbio);
    printf("+村庄<=%d格        : %llu\n", P_VIL_D, (unsigned long long)ctx.n_vil);
    printf("+河流连通块        : %llu\n", (unsigned long long)ctx.n_river);
    printf("完全合规种子       : %d\n", ctx.nhits);
    if (outpath) fclose(ctx.out);
    return 0;
}
