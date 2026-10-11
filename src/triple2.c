/* triple2.c — triple.c 的查表优化版
 *
 * 关键优化：单个神殿对覆盖率的贡献只取决于"挂机点相对神殿中心的偏移" (u,v)，
 * 与神殿绝对位置无关。所以预先算一张表 T[u][v] = 该偏移下落入 128 格的方块数，
 * 之后每次覆盖率计算退化成几次查表相加。
 *
 * 另外把搜索窗口收缩到"距每个神殿中心 169 格以内"（128+41，41 是 58x58 的半对角线），
 * 超出这个距离的神殿贡献为 0，没必要搜。
 *
 * 输入: monpairs_sorted.txt（base oi oj cov4）
 * 输出: triple_out.txt（base oi oj a4 best3 drop x z）
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
#define MCV MC_1_18_2

#define TR 200                 /* 表半径：偏移范围 [-TR, TR] */
#define TDIM (2*TR+1)
static int16_t *TAB;           /* TDIM x TDIM */
static double RSQ;

static void build_table(void)
{
    TAB = malloc(sizeof(int16_t) * (size_t)TDIM * TDIM);
    if (!TAB) { perror("malloc"); exit(1); }
    for (int u = -TR; u <= TR; u++)
        for (int v = -TR; v <= TR; v++) {
            int cnt = 0;
            for (int px = 0; px < 58; px++)
                for (int pz = 0; pz < 58; pz++) {
                    double ex = px + u, ez = pz + v;
                    if (ex*ex + ez*ez <= RSQ) cnt++;
                }
            TAB[(u+TR)*TDIM + (v+TR)] = (int16_t)cnt;
        }
}

static inline int cov_at(Pos *p, int n, int X, int Z)
{
    int s = 0;
    for (int i = 0; i < n; i++) {
        int u = p[i].x - X, v = p[i].z - Z;
        if (u < -TR || u > TR || v < -TR || v > TR) continue;   /* 贡献为 0 */
        s += TAB[(u+TR)*TDIM + (v+TR)];
    }
    return s;
}

/* 在收缩窗口内穷举 + 候选点，取最大覆盖率 */
static int afk_area(Pos *p, int n, int *outx, int *outz)
{
    int minX = 1<<30, maxX = -(1<<30), minZ = 1<<30, maxZ = -(1<<30);
    for (int i = 0; i < n; i++) {
        if (p[i].x < minX) minX = p[i].x;
        if (p[i].x > maxX) maxX = p[i].x;
        if (p[i].z < minZ) minZ = p[i].z;
        if (p[i].z > maxZ) maxZ = p[i].z;
    }
    /* 完整包围盒：最优挂机点可能在某座神殿 169 格之外（那座贡献 0 但其余更高），
     * 收缩窗口会漏掉这种情形，所以必须搜全盒。查表足够便宜，穷举得起。 */
    /* 对齐 cubiomes：它把包围盒平移 +ax/2(+29) 后再搜；这里放宽到 +-64 覆盖两种约定 */
    int loX = minX - 64, hiX = maxX + 64, loZ = minZ - 64, hiZ = maxZ + 64;
    int best = -1, bx = 0, bz = 0;

    if (loX <= hiX && loZ <= hiZ) {
        for (int X = loX; X <= hiX; X++)
            for (int Z = loZ; Z <= hiZ; Z++) {
                int c = cov_at(p, n, X, Z);
                if (c > best) { best = c; bx = X; bz = Z; }
            }
    }
    /* 候选点兜底（窗口为空或边界情形） */
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            int X = (p[i].x + p[j].x) / 2, Z = (p[i].z + p[j].z) / 2;
            int c = cov_at(p, n, X, Z);
            if (c > best) { best = c; bx = X; bz = Z; }
        }
    { int sx = 0, sz = 0;
      for (int i = 0; i < n; i++) { sx += p[i].x; sz += p[i].z; }
      int X = sx/n, Z = sz/n, c = cov_at(p, n, X, Z);
      if (c > best) { best = c; bx = X; bz = Z; } }

    if (outx) *outx = bx;
    if (outz) *outz = bz;
    return best;
}

typedef struct { uint64_t base; int oi, oj, cov4; } Pair;
static Pair *pairs; static int npairs;
static int NEXT; static pthread_mutex_t MTX = PTHREAD_MUTEX_INITIALIZER;
static FILE *PF; static pthread_mutex_t RTX = PTHREAD_MUTEX_INITIALIZER;
static int hist4[21], hist3[21];

static void *worker(void *arg)
{
    (void)arg;
    static const int rc[4][2] = {{0,0},{0,1},{1,0},{1,1}};
    for (;;) {
        pthread_mutex_lock(&MTX); int idx = NEXT++; pthread_mutex_unlock(&MTX);
        if (idx >= npairs) break;
        Pair pr = pairs[idx];
        uint64_t base = pr.base & 0xFFFFFFFFFFFFULL;
        Pos p[4];
        for (int i = 0; i < 4; i++)
            p[i] = getStructurePos(Monument, MCV, base, pr.oi+rc[i][0], pr.oj+rc[i][1], &p[i]) ? p[i] : p[i];

        int a4x, a4z;
        int a4 = afk_area(p, 4, &a4x, &a4z);
        int best3 = 0, bdrop = -1, b3x = 0, b3z = 0;
        for (int drop = 0; drop < 4; drop++) {
            Pos q[4]; int n = 0;
            for (int i = 0; i < 4; i++) if (i != drop) q[n++] = p[i];
            int tx, tz;
            int a = afk_area(q, 3, &tx, &tz);
            if (a > best3) { best3 = a; bdrop = drop; b3x = tx; b3z = tz; }
        }
        pthread_mutex_lock(&RTX);
        hist4[a4*20/13456 > 20 ? 20 : a4*20/13456]++;
        hist3[best3*20/10092 > 20 ? 20 : best3*20/10092]++;
        if (PF) fprintf(PF, "%llu %d %d %d %d %d %d %d\n",
                        (unsigned long long)base, pr.oi, pr.oj, a4, best3, bdrop, b3x, b3z);
        pthread_mutex_unlock(&RTX);
    }
    return NULL;
}

int main(int argc, char **argv)
{
    int NT = argc > 1 ? atoi(argv[1]) : 24;
    RSQ = 128.0*128.0;
    build_table();

    FILE *f = fopen("monpairs_sorted.txt", "r");
    if (!f) { perror("monpairs_sorted.txt"); return 1; }
    pairs = malloc(sizeof(Pair) * 30000); npairs = 0;
    while (npairs < 30000) {
        unsigned long long b; int oi, oj, cov;
        if (fscanf(f, "%llu %d %d %d", &b, &oi, &oj, &cov) != 4) break;
        pairs[npairs].base = b; pairs[npairs].oi = oi;
        pairs[npairs].oj = oj; pairs[npairs].cov4 = cov;
        npairs++;
    }
    fclose(f);
    PF = fopen("triple_out.txt", "w");
    if (!PF) { perror("triple_out.txt"); return 1; }

    printf("处理 %d 个配对（每配对算 4联 + 4 种 3联，查表法）\n", npairs);
    pthread_t th[256];
    for (int i = 0; i < NT; i++) pthread_create(&th[i], NULL, worker, NULL);
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    fclose(PF);

    printf("\n4 联覆盖率分布（分母 13456）:\n");
    for (int i = 0; i <= 20; i++) if (hist4[i]) printf("  %3d-%3d%% : %6d\n", i*5, i*5+4, hist4[i]);
    printf("\n3 联覆盖率分布（分母 10092）:\n");
    for (int i = 0; i <= 20; i++) if (hist3[i]) printf("  %3d-%3d%% : %6d\n", i*5, i*5+4, hist3[i]);
    return 0;
}
