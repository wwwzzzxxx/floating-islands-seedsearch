/* near.c — 半径 1000 格内找"最佳3座神殿组"，不要求 2x2 四联
 *
 * 需求（用户 2026-10-11 确认）：
 *   - 4 女巫小屋在半径 1000 格内（小屋四联锚定原点区域，天然满足）
 *   - 神殿组（按覆盖率取最佳 3 座）在半径 1000 格内
 *   - 神殿刷怪面积 >= 8500
 *   - 出生点有村庄
 *
 * 与 pairs.c 的区别：pairs.c 用 isQuadBase 要求 4 座全在 256 格内（2x2 四联），
 * 会排除"3 座很近、第 4 座在 256 格外"的排列。本工具对每个 2x2 区域块里的
 * 4 座神殿，直接评估 4 种"取3座"的组合，不做四联判定。
 *
 * 输出: near_out.txt  base oi oj cov x z（cov 为最佳3座覆盖率，分母 10092）
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

#define TR 200
#define TDIM (2*TR+1)
static int16_t *TAB;

static void build_table(void)
{
    TAB = malloc(sizeof(int16_t) * (size_t)TDIM * TDIM);
    for (int u = -TR; u <= TR; u++)
        for (int v = -TR; v <= TR; v++) {
            int c = 0;
            for (int px = 0; px < 58; px++)
                for (int pz = 0; pz < 58; pz++) {
                    double ex = px + u, ez = pz + v;
                    if (ex*ex + ez*ez <= 16384.0) c++;
                }
            TAB[(u+TR)*TDIM + (v+TR)] = (int16_t)c;
        }
}

static inline int cov_at(Pos *p, int n, int X, int Z)
{
    int s = 0;
    for (int i = 0; i < n; i++) {
        int u = p[i].x - X, v = p[i].z - Z;
        if (u < -TR || u > TR || v < -TR || v > TR) continue;
        s += TAB[(u+TR)*TDIM + (v+TR)];
    }
    return s;
}

static int afk_area(Pos *p, int n, int *outx, int *outz)
{
    int minX = 1<<30, maxX = -(1<<30), minZ = 1<<30, maxZ = -(1<<30);
    for (int i = 0; i < n; i++) {
        if (p[i].x < minX) minX = p[i].x;
        if (p[i].x > maxX) maxX = p[i].x;
        if (p[i].z < minZ) minZ = p[i].z;
        if (p[i].z > maxZ) maxZ = p[i].z;
    }
    int best = 0, bx = 0, bz = 0;
    for (int X = minX - 29; X <= maxX + 29; X++)
        for (int Z = minZ - 29; Z <= maxZ + 29; Z++) {
            int c = cov_at(p, n, X, Z);
            if (c > best) { best = c; bx = X; bz = Z; }
        }
    if (outx) *outx = bx;
    if (outz) *outz = bz;
    return best;
}

static uint64_t *bases; static long nbases;
static int OFF = 2;        /* 区域偏移范围 */
static int RADIUS = 1000;  /* 神殿组距原点上限 */
static int MAXD = 360;     /* 三元组两两距离上限（剪枝） */
static int NEXT; static pthread_mutex_t MTX = PTHREAD_MUTEX_INITIALIZER;
static FILE *PF; static pthread_mutex_t RTX = PTHREAD_MUTEX_INITIALIZER;
static int hist[21];
static long nwith;

static void *worker(void *arg)
{
    (void)arg;
    static const int rc[4][2] = {{0,0},{0,1},{1,0},{1,1}};
    for (;;) {
        pthread_mutex_lock(&MTX); long idx = NEXT++; pthread_mutex_unlock(&MTX);
        if (idx >= nbases) break;
        uint64_t base = bases[idx] & 0xFFFFFFFFFFFFULL;
        int best = 0, boi = 0, boj = 0, bax = 0, baz = 0;

        for (int oi = -OFF; oi <= OFF; oi++)
        for (int oj = -OFF; oj <= OFF; oj++) {
            Pos p[4]; int ok = 1;
            for (int i = 0; i < 4; i++) {
                if (!getStructurePos(Monument, MCV, base, oi+rc[i][0], oj+rc[i][1], &p[i])) { ok = 0; break; }
                /* 神殿须在半径 1000 格内 */
                if (p[i].x*p[i].x + p[i].z*p[i].z > (int64_t)RADIUS*RADIUS) { ok = 0; break; }
            }
            if (!ok) continue;
            int64_t md2 = (int64_t)MAXD * MAXD;
            for (int drop = 0; drop < 4; drop++) {
                Pos q[4]; int n = 0;
                for (int i = 0; i < 4; i++) if (i != drop) q[n++] = p[i];
                /* 剪枝：两两距离超限则不可能高覆盖 */
                int64_t d01 = (int64_t)(q[0].x-q[1].x)*(q[0].x-q[1].x) + (int64_t)(q[0].z-q[1].z)*(q[0].z-q[1].z);
                int64_t d02 = (int64_t)(q[0].x-q[2].x)*(q[0].x-q[2].x) + (int64_t)(q[0].z-q[2].z)*(q[0].z-q[2].z);
                int64_t d12 = (int64_t)(q[1].x-q[2].x)*(q[1].x-q[2].x) + (int64_t)(q[1].z-q[2].z)*(q[1].z-q[2].z);
                if (d01 > md2 || d02 > md2 || d12 > md2) continue;
                int tx, tz;
                int a = afk_area(q, 3, &tx, &tz);
                if (a > best) { best = a; boi = oi; boj = oj; bax = tx; baz = tz; }
            }
        }
        pthread_mutex_lock(&RTX);
        if (best > 0) {
            nwith++;
            hist[best*20/10092 > 20 ? 20 : best*20/10092]++;
            if (PF) fprintf(PF, "%llu %d %d %d %d %d\n",
                            (unsigned long long)base, boi, boj, best, bax, baz);
        }
        pthread_mutex_unlock(&RTX);
    }
    return NULL;
}

int main(int argc, char **argv)
{
    int NT = argc > 1 ? atoi(argv[1]) : 16;
    if (argc > 2) OFF = atoi(argv[2]);
    if (argc > 3) RADIUS = atoi(argv[3]);
    if (argc > 4) MAXD = atoi(argv[4]);
    build_table();

    FILE *f = fopen("hutbases_local.txt", "r");
    if (!f) { perror("hutbases_local.txt"); return 1; }
    bases = malloc(sizeof(uint64_t) * 800000); nbases = 0;
    while (nbases < 800000 && fscanf(f, "%llu", (unsigned long long*)&bases[nbases]) == 1) nbases++;
    fclose(f);

    PF = fopen("near_out.txt", "w");
    printf("%ld 个小屋基，神殿区域偏移 ±%d，半径 %d 格内，三元组两两<=%d\n",
           nbases, OFF, RADIUS, MAXD);
    pthread_t th[256];
    for (int i = 0; i < NT; i++) pthread_create(&th[i], NULL, worker, NULL);
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    if (PF) fclose(PF);

    printf("\n%ld 个基找到符合条件的 3 座组\n", nwith);
    printf("覆盖率分布（分母 10092）:\n");
    for (int i = 0; i <= 20; i++) if (hist[i]) printf("  %3d-%3d%% : %6ld\n", i*5, i*5+4, hist[i]);
    return 0;
}
