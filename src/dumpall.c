/* dumpall.c — 导出每个种子的完整坐标（机器可读，一行一种子）
 *
 * 输入: 种子清单，每行 "seed ... 覆盖率 ... (oi,oj)"（兼容 seeds.txt）
 * 输出: seeds_full.txt，每行：
 *   seed hut_r cov vil_d river_blob
 *   hut:(x,z)x4 | mon:(x,z)x4 | afk3:(x,z) | vil:(x,z) | river_c:(x,z) | fort:(x,z) dist
 *
 * afk3 = 4 座神殿中取最佳 3 座时的最优挂机点（口径对齐 SunnySlopes pairLevel=3）。
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

#define TR 200
#define TDIM (2*TR+1)
static int16_t *TAB;
static int biome_river;

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

/* 最佳 3 座的覆盖率与挂机点 */
static int best3(Pos *p, int *ax, int *az)
{
    int best = 0, bx = 0, bz = 0;
    for (int drop = 0; drop < 4; drop++) {
        Pos q[4]; int n = 0;
        for (int i = 0; i < 4; i++) if (i != drop) q[n++] = p[i];
        int minX = 1<<30, maxX = -(1<<30), minZ = 1<<30, maxZ = -(1<<30);
        for (int i = 0; i < 3; i++) {
            if (q[i].x < minX) minX = q[i].x;
            if (q[i].x > maxX) maxX = q[i].x;
            if (q[i].z < minZ) minZ = q[i].z;
            if (q[i].z > maxZ) maxZ = q[i].z;
        }
        for (int X = minX - 29; X <= maxX + 29; X++)
            for (int Z = minZ - 29; Z <= maxZ + 29; Z++) {
                int c = cov_at(q, 3, X, Z);
                if (c > best) { best = c; bx = X; bz = Z; }
            }
    }
    if (ax) *ax = bx;
    if (az) *az = bz;
    return best;
}

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

static int is_crossing(int type)
{
    enum { F_BRIDGE_CROSSING=1, F_BRIDGE_FORTIFIED_CROSSING=2,
           F_CORRIDOR_CROSSING=7, F_CORRIDOR_T_CROSSING=11 };
    return type==F_BRIDGE_CROSSING || type==F_BRIDGE_FORTIFIED_CROSSING ||
           type==F_CORRIDOR_CROSSING || type==F_CORRIDOR_T_CROSSING;
}

static int find_fortress(uint64_t seed, int range, int *fx, int *fz, int *fdist)
{
    Generator gn; setupGenerator(&gn, MCV, 0);
    applySeed(&gn, DIM_NETHER, seed);
    StructureConfig fc; getStructureConfig(Fortress, MCV, &fc);
    int regsize = fc.regionSize << 4;
    int rmax = range/regsize + 2;
    long best = 0x7fffffff;
    for (int rz = -rmax; rz <= rmax; rz++)
    for (int rx = -rmax; rx <= rmax; rx++) {
        Pos fp;
        if (!getStructurePos(Fortress, MCV, seed, rx, rz, &fp)) continue;
        long dd = d2(fp.x, fp.z, 0, 0);
        if (dd > (long)range*range) continue;
        if (!isViableStructurePos(Fortress, &gn, fp.x, fp.z, 0)) continue;
        Piece *L = malloc(sizeof(Piece)*512);
        int n = getFortressPieces(L, 512, MCV, seed, fp.x>>4, fp.z>>4);
        if (n > 0) {
            int c[256][2], nc = 0;
            for (int i = 0; i < n && nc < 256; i++)
                if (is_crossing(L[i].type)) {
                    c[nc][0] = (L[i].bb0.x + L[i].bb1.x) / 2;
                    c[nc][1] = (L[i].bb0.z + L[i].bb1.z) / 2;
                    nc++;
                }
            for (int i = 0; i < nc; i++) {
                int h[4] = {0,0,0,0};
                for (int m = 0; m < nc; m++) {
                    if (c[m][0]==c[i][0]     && c[m][1]==c[i][1])     h[0]=1;
                    if (c[m][0]==c[i][0]+19  && c[m][1]==c[i][1])     h[1]=1;
                    if (c[m][0]==c[i][0]     && c[m][1]==c[i][1]+19)  h[2]=1;
                    if (c[m][0]==c[i][0]+19  && c[m][1]==c[i][1]+19)  h[3]=1;
                }
                if (h[0]&&h[1]&&h[2]&&h[3] && dd < best) {
                    best = dd;
                    *fx = c[i][0] + 9; *fz = c[i][1] + 9;
                    *fdist = (int)sqrt((double)dd);
                }
            }
        }
        free(L);
    }
    return best < 0x7fffffff;
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "用法: %s <种子清单> [输出文件] [要塞范围]\n", argv[0]); return 1; }
    const char *outpath = argc > 2 ? argv[2] : "seeds_full.txt";
    int range = argc > 3 ? atoi(argv[3]) : 4000;
    build_table();
    biome_river = river;   /* 必须显式赋值：静态初值 0 是 ocean，会把它当成河流来数 */

    FILE *in = fopen(argv[1], "r");
    if (!in) { perror(argv[1]); return 1; }
    FILE *out = fopen(outpath, "w");
    if (!out) { perror(outpath); return 1; }
    fprintf(out, "# seed hut_r cov vil_d river_blob | hut:(x,z)x4 | mon:(x,z)x4 | afk3 | vil | river_c | fort\n");

    char line[512];
    int n = 0;
    while (fgets(line, sizeof line, in)) {
        if (line[0] == '#') continue;
        /* seed hr cov ... (oi,oj) —— 第4列是覆盖率 */
        long long seedsigned;
        double hr = 0, cov = 0;
        int oi = 0, oj = 0;
        int k = sscanf(line, "%lld %lf %lf", &seedsigned, &hr, &cov);
        if (k < 1) continue;
        char *p1 = strchr(line, '(');
        char *p2 = p1 ? strchr(p1, ')') : NULL;
        if (p1 && p2 && p2 > p1) {
            char buf[64]; int len = (int)(p2 - p1 - 1);
            if (len > 0 && len < 63) {
                memcpy(buf, p1+1, len); buf[len] = 0;
                char *comma = strchr(buf, ',');
                if (comma) { *comma = 0; oi = atoi(buf); oj = atoi(comma+1); }
            }
        }
        uint64_t seed = (uint64_t)seedsigned;

        Generator g; setupGenerator(&g, MCV, 0);
        applySeed(&g, DIM_OVERWORLD, seed);   /* 必须紧跟 setupGenerator，与 coords.c 一致 */
        static const int rc[4][2] = {{0,0},{0,1},{1,0},{1,1}};
        Pos hp[4], mp[4];
        for (int i = 0; i < 4; i++) {
            getStructurePos(Swamp_Hut, MCV, seed, rc[i][0], rc[i][1], &hp[i]);
            getStructurePos(Monument, MCV, seed, oi+rc[i][0], oj+rc[i][1], &mp[i]);
        }
        int afx, afz, cov3 = best3(mp, &afx, &afz);

        /* 注意：必须先做结构可行性检查再算河流。
         * isViableStructurePos 会初始化生成器的群系采样路径，
         * 跳过它直接 genBiomes 会得到不同的结果（实测连通块 2507 -> 5610）。 */
        for (int i = 0; i < 4; i++) {
            isViableStructurePos(Swamp_Hut, &g, hp[i].x, hp[i].z, 0);
            isViableStructurePos(Monument, &g, mp[i].x, mp[i].z, 0);
        }
        int bd = 0x7fffffff, vx = 0, vz = 0, found = 0;
        for (int gz=-1; gz<=1; gz++) for (int gx=-1; gx<=1; gx++) {
            Pos vp;
            if (!getStructurePos(Village, MCV, seed, gx, gz, &vp)) continue;
            int dd = d2(vp.x, vp.z, 0, 0);
            if (dd >= bd) continue;
            if (!isViableStructurePos(Village, &g, vp.x, vp.z, 0)) continue;
            bd = dd; vx = vp.x; vz = vp.z; found = 1;
        }
        int vild = found ? (int)sqrt((double)bd) : -1;
        double frac, rcx, rcz;
        int blob = river_blob(&g, 3000, &frac, &rcx, &rcz);

        int fx = 0, fz = 0, fd = 0;
        int hasf = find_fortress(seed, range, &fx, &fz, &fd);

        fprintf(out, "%lld %.1f %.0f %d %d |", seedsigned, hr, cov, vild, blob);
        for (int i = 0; i < 4; i++) fprintf(out, " (%d,%d)", hp[i].x, hp[i].z);
        fprintf(out, " |");
        for (int i = 0; i < 4; i++) fprintf(out, " (%d,%d)", mp[i].x, mp[i].z);
        fprintf(out, " | (%d,%d) | (%d,%d) | (%.0f,%.0f) |", afx, afz, vx, vz, rcx, rcz);
        if (hasf) fprintf(out, " (%d,%d) %d\n", fx, fz, fd);
        else      fprintf(out, " none\n");

        /* 人类可读的 Markdown 行 */
        printf("| `%lld` | `%d 200 %d` / `%d 200 %d` / `%d 200 %d` / `%d 200 %d` | "
               "`%d 200 %d` / `%d 200 %d` / `%d 200 %d` / `%d 200 %d` | "
               "**`%d 200 %d`** | `%d 200 %d`（%d 格） | `%.0f 200 %.0f`（%d） | "
               "`%d 200 %d`（%d 格） |\n",
               seedsigned,
               hp[0].x, hp[0].z, hp[1].x, hp[1].z, hp[2].x, hp[2].z, hp[3].x, hp[3].z,
               mp[0].x, mp[0].z, mp[1].x, mp[1].z, mp[2].x, mp[2].z, mp[3].x, mp[3].z,
               afx, afz, vx, vz, vild, rcx, rcz, blob,
               hasf ? fx : 0, hasf ? fz : 0, fd);
        n++;
    }
    fclose(in); fclose(out);
    fprintf(stderr, "已导出 %d 个种子 -> %s\n", n, outpath);
    return 0;
}
