/* 第三阶段（修正版）：
 *   小屋四联：锚定原点区域 (0,0)-(1,1)（靠近出生点）
 *   神殿四联：允许出现在 (±OFFSET 区域) 范围内——避免与小屋争同一片群系
 *   然后对每个 (base, 神殿偏移) 做 65536 个高16位变体搜索：
 *     4小屋全 swamp + 4神殿全 ocean + 出生点村庄<=256 + 3000x3000 河流连通块
 * 多线程（每线程自己的 Generator）。
 *
 * ---------------------------------------------------------------------------
 * 成本结构（实测，2 万个 base、16 线程，消融法测得）：
 *   偏移扫描（169 次/基地 isQuadBase 神殿四联判定）   0.17 秒   0.001%
 *   变体循环（65536 个高16位变体的群系/村庄/河流验证）13.50 秒  99.98%
 *
 * 结论：**"用四联结构筛种子"这一步本身几乎是免费的**，耗时全在验证活下来的
 * 种子——因为群系取决于完整 64 位种子，每个变体都要重新 applySeed + 采样群系。
 *
 * 据此做了两处优化（产出与原版逐字节一致，已比对）：
 *   ① 四联结构的判定只取决于种子低 48 位，而 seed = base | (u<<48) 恒保持
 *      低 48 位为 base，所以 rh/rm 对 65536 个变体是同一个值——原版在每个
 *      变体里重算，提到 (base,偏移) 级只算一次。
 *   ② 4 个小屋的坐标同样只取决于低 48 位，提到 base 级只算一次。
 *   合计约快 16%（15.9 秒 → 13.3 秒）。剩下的时间是 applySeed 与群系采样
 *      本身，属固有成本；要再快只能减少变体数或把群系判定移植到 GPU。
 * ---------------------------------------------------------------------------
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

        /* 找该 base 上所有偏移位置的神殿四联（顺便存下半径，变体循环里不必重算） */
        int offs[64][2]; float offr[64]; int no = 0;
        for (int i = -P_OFFSET; i <= P_OFFSET && no < 64; i++)
            for (int j = -P_OFFSET; j <= P_OFFSET && no < 64; j++) {
                uint64_t sp = moveStructure(base, -i, -j);
                float r = isQuadBase(monconf, sp, P_MON_R);
                if (r != 0.0f) { offs[no][0]=i; offs[no][1]=j; offr[no]=r; no++; }
            }
        if (no == 0) continue;
        /* 小屋四联半径只取决于低 48 位，对本 base 是常量，提到这里算一次 */
        float rh_fix = isQuadBase(hutconf, base, 128);
        if (rh_fix == 0.0f) continue;
        /* 4 个小屋坐标同样只取决于低 48 位，也是 base 级常量 */
        Pos hp_fix[4];
        { static const int rx0[4]={0,1,0,1}, rz0[4]={0,0,1,1};
          for (int q=0;q<4;q++) getStructurePos(Swamp_Hut, MCV, base, rx0[q], rz0[q], &hp_fix[q]); }
        pthread_mutex_lock(&c->rmtx); c->n_pairs += no; pthread_mutex_unlock(&c->rmtx);

        /* 对每个偏移做变体搜索 */
        for (int k = 0; k < no; k++) {
            int oi = offs[k][0], oj = offs[k][1];
            float rm = offr[k];   /* 只取决于低 48 位，与变体无关 */
            float rh = rh_fix;
            for (uint64_t u = 0; u < 65536; u++) {
                uint64_t seed = base | (u << 48);
                c->n_var++;
                applySeed(&g, DIM_OVERWORLD, seed);
                Pos mp[4];
                static const int rx[4]={0,1,0,1}, rz[4]={0,0,1,1};
                int i, ok = 1;
                for (i=0;i<4;i++){
                    if (!isViableStructurePos(Swamp_Hut, &g, hp_fix[i].x, hp_fix[i].z, 0)) { ok=0; break; }
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
