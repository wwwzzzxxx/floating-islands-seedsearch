/* 第四阶段：对候选种子做下界要塞「四联路口」分析
 * 路口定义（对齐 CrossroadFinder 的形状）：
 *   QUAD_SQUARE : 4 个路口成 2x2 方阵，相邻间距 19 格
 *   QUAD_LINE   : 4 个路口共线，相邻间距 19 格
 *   TRIPLE_*    : 3 个路口（次级目标）
 * 搜索原点附近 ±RANGE 格内的要塞。
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "generator.h"
#include "finders.h"
#include "quadbase.h"
#include "biomes.h"

#define MCV MC_1_18_2

static inline int d2(int x1,int z1,int x2,int z2){int dx=x1-x2,dz=z1-z2;return dx*dx+dz*dz;}

/* 判断是否为“路口”构件 */
static int is_crossing(int type)
{
    /* cubiomes 要塞构件类型（finders.h 枚举顺序） */
    enum {
        F_BRIDGE_STRAIGHT = 0, F_BRIDGE_CROSSING, F_BRIDGE_FORTIFIED_CROSSING,
        F_BRIDGE_STAIRS, F_BRIDGE_SPAWNER, F_BRIDGE_CORRIDOR_ENTRANCE,
        F_CORRIDOR_STRAIGHT, F_CORRIDOR_CROSSING, F_CORRIDOR_TURN_RIGHT,
        F_CORRIDOR_TURN_LEFT, F_CORRIDOR_STAIRS, F_CORRIDOR_T_CROSSING,
        F_CORRIDOR_NETHER_WART, F_END
    };
    return type == F_BRIDGE_CROSSING || type == F_BRIDGE_FORTIFIED_CROSSING ||
           type == F_CORRIDOR_CROSSING || type == F_CORRIDOR_T_CROSSING;
}

/* 在一组路口里找图案：
 *  *oline = 最长等距(19格)共线联数
 *  *osq   = 是否存在 2x2 方阵（19格间距）—— 这才是经典「四联路口」
 * 返回 max(oline, osq?4:0)
 */
static int find_pattern(Piece *list, int n, double *ox, double *oz, int *span, int *osq, int *oline)
{
    /* 用包围盒中心作为旋转不变量位置 */
    Pos3 cr[256]; int nc = 0;
    for (int i = 0; i < n && nc < 256; i++)
        if (is_crossing(list[i].type)) {
            cr[nc].x = (list[i].bb0.x + list[i].bb1.x) / 2;
            cr[nc].z = (list[i].bb0.z + list[i].bb1.z) / 2;
            cr[nc].y = 0;
            nc++;
        }

    /* 2x2 方阵 */
    int sq = 0; double sqx=0, sqz=0;
    for (int i = 0; i < nc; i++) {
        int x0 = cr[i].x, z0 = cr[i].z;
        int has[4] = {0,0,0,0};
        for (int m = 0; m < nc; m++) {
            if (cr[m].x==x0    && cr[m].z==z0)    has[0]=1;
            if (cr[m].x==x0+19 && cr[m].z==z0)    has[1]=1;
            if (cr[m].x==x0    && cr[m].z==z0+19) has[2]=1;
            if (cr[m].x==x0+19 && cr[m].z==z0+19) has[3]=1;
        }
        if (has[0]&&has[1]&&has[2]&&has[3]) { sq = 4; sqx = x0+9.5; sqz = z0+9.5; }
    }

    /* 最长等距共线 */
    int bestline = 0; double bx=0, bz=0; int bspan=0;
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            if (i == j) continue;
            int dx = cr[j].x - cr[i].x, dz = cr[j].z - cr[i].z;
            if (dx != 0 && dz != 0) continue;
            if (dx == 0 && dz == 0) continue;
            int stepx = (dx > 0) - (dx < 0), stepz = (dz > 0) - (dz < 0);
            int len = (dx != 0) ? (dx > 0 ? dx : -dx) : (dz > 0 ? dz : -dz);
            int cnt = 0, k;
            for (k = 0; k <= len/19 + 1; k++) {
                int px = cr[i].x + stepx*19*k, pz = cr[i].z + stepz*19*k;
                int found = 0;
                for (int m = 0; m < nc; m++)
                    if (cr[m].x == px && cr[m].z == pz) { found = 1; break; }
                if (found) cnt++; else break;
            }
            if (cnt > bestline) {
                bestline = cnt;
                bx = cr[i].x + stepx*19.0*(cnt-1)/2.0;
                bz = cr[i].z + stepz*19.0*(cnt-1)/2.0;
                bspan = 19*(cnt-1);
            }
        }
    }
    *osq = sq; *oline = bestline;
    if (sq) { *ox = sqx; *oz = sqz; *span = 19; return 4; }
    *ox = bx; *oz = bz; *span = bspan;
    return bestline;
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr,"用法: %s <seeds.txt> [-range N] [-out f]\n", argv[0]); return 1; }
    int range = 4000;
    const char *outpath = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i],"-range") && i+1<argc) range = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-out") && i+1<argc) outpath = argv[++i];
    }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("打开种子列表失败"); return 1; }
    uint64_t *seeds = NULL; int ns = 0, cap = 0; char line[128];
    while (fgets(line, sizeof line, f)) {
        char *p = line; while (*p==' ') p++;
        if (*p < '0' || *p > '9') { if (*p != '-') continue; }
        uint64_t v = strtoull(p, 0, 10);
        if (ns >= cap) { cap = cap?cap*2:1024; seeds = realloc(seeds, cap*sizeof(uint64_t)); }
        seeds[ns++] = v;
    }
    fclose(f);
    fprintf(stderr, "读入 %d 个种子，搜索要塞范围 ±%d 格\n", ns, range);

    StructureConfig fconf;
    getStructureConfig(Fortress, MCV, &fconf);
    Generator gn; setupGenerator(&gn, MCV, 0); applySeed(&gn, DIM_NETHER, 0);
    int regsize = fconf.regionSize << 4;   /* 432 */
    int rmax = range/regsize + 2;

    FILE *out = outpath ? fopen(outpath,"w") : stdout;
    Piece *list = malloc(sizeof(Piece)*512);
    int nquad = 0, ntriple = 0, nany = 0;
    for (int s = 0; s < ns; s++) {
        uint64_t seed = seeds[s];
        applySeed(&gn, DIM_NETHER, seed);
        int best = 0, bestline = 0, bestsq = 0; double bx=0,bz=0; int bspan=0, bdist=0; long bdd=0x7fffffff;
        for (int rz = -rmax; rz <= rmax; rz++)
        for (int rx = -rmax; rx <= rmax; rx++) {
            Pos fp;
            if (!getStructurePos(Fortress, MCV, seed, rx, rz, &fp)) continue;
            long dd = d2(fp.x, fp.z, 0, 0);
            if (dd > range*range) continue;
            /* 1.18+ 堡垒与要塞同区域同位置二选一，必须先确认要塞真的生成 */
            if (!isViableStructurePos(Fortress, &gn, fp.x, fp.z, 0)) continue;
            int n = getFortressPieces(list, 512, MCV, seed, fp.x>>4, fp.z>>4);
            if (n <= 0) continue;
            double ox, oz; int span, sq, line;
            int cnt = find_pattern(list, n, &ox, &oz, &span, &sq, &line);
            if (line > bestline) bestline = line;
            if (sq && (!bestsq || dd < bdd)) {
                best = 4; bestsq = 1; bx = ox; bz = oz; bdd = dd;
                bdist = (int)sqrt((double)dd);
            }
            else if (!bestsq && line >= 4 && line >= bestline && dd < bdd) {
                bx = ox; bz = oz; bdd = dd; bdist = (int)sqrt((double)dd);
            }
        }
        if (best >= 4) {
            nquad++;
            fprintf(out, "%lld SQUARE pos=(%.0f,%.0f) dist=%d maxline=%d\n",
                    (long long)(int64_t)seed, bx, bz, bdist, bestline);
        } else if (bestline >= 4) {
            nquad++;
            fprintf(out, "%lld LINE%d pos=(%.0f,%.0f) dist=%d\n",
                    (long long)(int64_t)seed, bestline, bx, bz, bdist);
        } else if (bestline == 3) {
            ntriple++;
        } else nany++;
        if ((s+1) % 200 == 0) { fprintf(stderr, "  已处理 %d/%d（四联 %d，三联 %d）\n", s+1, ns, nquad, ntriple); }
    }
    printf("\n[要塞路口统计] 共 %d 个种子\n", ns);
    printf("  有四联路口 : %d\n", nquad);
    printf("  仅三联     : %d\n", ntriple);
    printf("  无多联路口 : %d\n", nany);
    if (outpath) fclose(out);
    return 0;
}
