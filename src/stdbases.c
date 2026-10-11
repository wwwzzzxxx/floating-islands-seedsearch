/* stdbases.c — 评估"放弃四联小屋、改用 121 个社区标准神殿基"的可行性
 *
 * 对每个标准基 M 和每个区域偏移，反推种子 s = moveStructure(M, oi, oj)，
 * 然后在 65536 个高16位变体里找第一个同时满足：
 *    出生点 256 格内有村庄 + 3000x3000 有河流连通块
 * （不要求四联小屋——这正是"限制条件小"的场景）
 * 找到即早停，所以很快。
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
#include "seedtables.h"
#define MCV MC_1_18_2

static int P_VIL_D = 256, P_RIVER_SPAN = 3000;
static int biome_river;

typedef struct {
    uint64_t mon;      /* 标准基（salt 空间） */
    int qual;          /* viewer 给的覆盖率 */
    int oi, oj;
} Job;

static Job *jobs; static int njobs;
static int NEXT; static pthread_mutex_t MTX = PTHREAD_MUTEX_INITIALIZER;
static FILE *PF; static pthread_mutex_t RTX = PTHREAD_MUTEX_INITIALIZER;
static long nfound, ntested;

static inline int d2(int x1,int z1,int x2,int z2){int dx=x1-x2,dz=z1-z2;return dx*dx+dz*dz;}

static int river_blob(Generator *g, int span, double *frac)
{
    int s = span/4; if (s<1) s=1;
    Range r = {4, -s/2, -s/2, s, s};
    int *cache = malloc(sizeof(int)*(size_t)s*s);
    if (!cache) return -1;
    if (genBiomes(g, cache, r)!=0){ free(cache); return -1; }
    int total=s*s, cnt=0;
    for (int i=0;i<total;i++) if (cache[i]==biome_river) cnt++;
    *frac=(double)cnt/total;
    char *vis=calloc(total,1); int *st=malloc(sizeof(int)*total);
    int best=0;
    if (vis&&st){
        for (int i=0;i<total;i++){
            if (cache[i]!=biome_river||vis[i]) continue;
            int sp=0,size=0; st[sp++]=i; vis[i]=1;
            while(sp>0){ int c=st[--sp]; int x=c%s,z=c/s; size++;
                int k,nx,nz;
                for(k=0;k<4;k++){ nx=x+(k==0?1:k==1?-1:0); nz=z+(k==2?1:k==3?-1:0);
                    if(nx<0||nz<0||nx>=s||nz>=s) continue;
                    int ni=nz*s+nx; if(!vis[ni]&&cache[ni]==biome_river){vis[ni]=1;st[sp++]=ni;} } }
            if (size>best) best=size;
        }
    }
    free(vis); free(st); free(cache);
    return best;
}

static void *worker(void *arg)
{
    (void)arg;
    Generator g; setupGenerator(&g, MCV, 0);
    for (;;) {
        pthread_mutex_lock(&MTX); int idx = NEXT++; pthread_mutex_unlock(&MTX);
        if (idx >= njobs) break;
        Job jb = jobs[idx];
        uint64_t seed48 = moveStructure(jb.mon, jb.oi, jb.oj);

        for (uint64_t u = 0; u < 65536; u++) {
            uint64_t seed = seed48 | (u << 48);
            applySeed(&g, DIM_OVERWORLD, seed);
            /* 村庄 */
            int bd = 0x7fffffff, found = 0;
            for (int gz=-1; gz<=1; gz++) for (int gx=-1; gx<=1; gx++) {
                Pos vp;
                if (!getStructurePos(Village, MCV, seed, gx, gz, &vp)) continue;
                int dd = d2(vp.x, vp.z, 0, 0);
                if (dd >= bd) continue;
                if (!isViableStructurePos(Village, &g, vp.x, vp.z, 0)) continue;
                bd = dd; found = 1;
            }
            if (!found || bd > P_VIL_D*P_VIL_D) continue;
            double frac;
            int blob = river_blob(&g, P_RIVER_SPAN, &frac);
            if (blob <= 0) continue;
            pthread_mutex_lock(&RTX);
            nfound++; ntested += u + 1;
            if (PF) fprintf(PF, "%lld %d %d %d %d %.3f %d\n",
                            (long long)(int64_t)seed, jb.qual, jb.oi, jb.oj,
                            (int)sqrt((double)bd), frac*100, blob);
            pthread_mutex_unlock(&RTX);
            break;   /* 找到一个就够了 */
        }
    }
    return NULL;
}

int main(int argc, char **argv)
{
    int NT = argc > 1 ? atoi(argv[1]) : 24;
    int OFF = argc > 2 ? atoi(argv[2]) : 16;
    { Generator g; setupGenerator(&g, MCV, 0); applySeed(&g, DIM_OVERWORLD, 0);
      biome_river = river; }

    StructureConfig mc; getStructureConfig(Monument, MCV, &mc);
    jobs = malloc(sizeof(Job) * 200000); njobs = 0;
    for (const uint64_t *s = g_qm_90; *s; s++) {
        uint64_t M = (*s - mc.salt) & 0xFFFFFFFFFFFFULL;
        int q = qmonumentQual(*s);
        for (int oi = -OFF; oi <= OFF; oi++)
            for (int oj = -OFF; oj <= OFF; oj++) {
                jobs[njobs].mon = M; jobs[njobs].qual = q;
                jobs[njobs].oi = oi; jobs[njobs].oj = oj;
                njobs++;
            }
    }
    printf("121 个社区标准神殿基 x %d 偏移 = %d 个候选\n", (2*OFF+1)*(2*OFF+1), njobs);
    PF = fopen("stdbases_out.txt", "w");

    pthread_t th[256];
    for (int i = 0; i < NT; i++) pthread_create(&th[i], NULL, worker, NULL);
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    if (PF) fclose(PF);
    printf("找到 %ld 个满足村庄+河流的种子（平均每个候选只试了 %.0f 个变体）\n",
           nfound, njobs ? (double)ntested/njobs : 0);
    return 0;
}
