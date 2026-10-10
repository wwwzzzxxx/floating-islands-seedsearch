/* 输出推荐种子的全部关键坐标，方便 tp 过去看 */
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
static int biome_river;
static inline int d2(int x1,int z1,int x2,int z2){int dx=x1-x2,dz=z1-z2;return dx*dx+dz*dz;}

static void river_blob(Generator *g, int span, double *frac, int *blob, double *cx, double *cz)
{
    int s = span/4; if (s<1) s=1;
    Range r = {4, -s/2, -s/2, s, s};
    int *cache = malloc(sizeof(int)*(size_t)s*s);
    *blob=-1; *frac=-1;
    if (!cache) return;
    if (genBiomes(g, cache, r)!=0){ free(cache); return; }
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
    *blob = best;
    *cx=(bx-s/2)*4.0; *cz=(bz-s/2)*4.0;
}

int main(int argc, char **argv)
{
    /* 用法: coords <seed> <神殿偏移oi> <神殿偏移oj> */
    if (argc < 4) { fprintf(stderr,"用法: %s <seed> <oi> <oj>\n", argv[0]); return 1; }
    uint64_t seed = strtoull(argv[1],0,10);
    int oi = atoi(argv[2]), oj = atoi(argv[3]);
    Generator g; setupGenerator(&g, MCV, 0); applySeed(&g, DIM_OVERWORLD, seed);
    biome_river = river;
    StructureConfig fc; getStructureConfig(Fortress, MCV, &fc);

    printf("======== 种子 %lld ========\n", (long long)(int64_t)seed);
    printf("（tp 用 /tp <x> <y> <z>；Y 建议先 200 再往下看，浮岛世界虚空在下）\n\n");

    printf("[主世界] 出生点附近 (0,0)\n");
    printf("  原点 /tp 0 200 0\n\n");

    printf("[主世界] 四联女巫小屋（ swamp 群系，4 个都在 128 格内）\n");
    static const int rx[4]={0,1,0,1}, rz[4]={0,0,1,1};
    for (int i=0;i<4;i++){
        Pos p; getStructurePos(Swamp_Hut, MCV, seed, rx[i], rz[i], &p);
        int v = isViableStructurePos(Swamp_Hut, &g, p.x, p.z, 0);
        printf("  小屋%d  /tp %d 200 %d   (%s)\n", i+1, p.x, p.z, v?"在 swamp ✓":"不在 swamp ✗");
    }

    printf("\n[主世界] 四联海底神殿（区域偏移 (%d,%d)，4 个都在 256 格内）\n", oi, oj);
    for (int i=0;i<4;i++){
        Pos p; getStructurePos(Monument, MCV, seed, oi+rx[i], oj+rz[i], &p);
        int v = isViableStructurePos(Monument, &g, p.x, p.z, 0);
        printf("  神殿%d  /tp %d 200 %d   (%s)\n", i+1, p.x, p.z, v?"在 ocean ✓":"不在 ocean ✗");
    }

    printf("\n[主世界] 出生点村庄\n");
    int bd=0x7fffffff; Pos bv; int found=0;
    for (int gz=-1;gz<=1;gz++) for (int gx=-1;gx<=1;gx++){
        Pos vp;
        if(!getStructurePos(Village,MCV,seed,gx,gz,&vp)) continue;
        int dd=d2(vp.x,vp.z,0,0);
        if(dd>=bd) continue;
        if(!isViableStructurePos(Village,&g,vp.x,vp.z,0)) continue;
        bd=dd; bv=vp; found=1;
    }
    if (found) printf("  村庄  /tp %d 200 %d   距原点 %.0f 格\n", bv.x, bv.z, sqrt((double)bd));
    else printf("  未找到\n");

    double frac, rcx, rcz; int blob;
    river_blob(&g, 3000, &frac, &blob, &rcx, &rcz);
    printf("\n[主世界] 超大河流群（3000x3000 扫描）\n");
    printf("  覆盖率 %.2f%%   最大连通块 %d 个群系单元(约 %.0fx%.0f 格)\n", frac*100, blob, sqrt((double)blob)*4, sqrt((double)blob)*4);
    printf("  河流重心  /tp %.0f 200 %.0f\n", rcx, rcz);

    printf("\n[下界] 要塞四联路口（±4000 格内最优）\n");
    Piece *L = malloc(sizeof(Piece)*512);
    Generator gn; setupGenerator(&gn, MCV, 0); applySeed(&gn, DIM_NETHER, seed);
    int regsize = fc.regionSize<<4, rmax = 4000/regsize+2;
    int bestline=0, bestsq=0; int bdist=0;
    double sqx=0, sqz=0; long bsq_dd=0x7fffffff;      /* 方阵：最近的位置与距离 */
    double lnx=0, lnz=0; long bln_dd=0x7fffffff;      /* 共线：最长/最近 */
    int nfort=0, nbast=0;
    for (int rz2=-rmax; rz2<=rmax; rz2++) for (int rx2=-rmax; rx2<=rmax; rx2++){
        Pos fp;
        if (!getStructurePos(Fortress, MCV, seed, rx2, rz2, &fp)) continue;
        long dd = d2(fp.x, fp.z, 0, 0);          /* 平方距离 */
        if (dd > 4000L*4000L) continue;
        /* 1.18+ 堡垒与要塞同区域同位置，二选一：必须确认要塞真的可行 */
        if (!isViableStructurePos(Fortress, &gn, fp.x, fp.z, 0)) { nbast++; continue; }
        nfort++;
        int n = getFortressPieces(L, 512, MCV, seed, fp.x>>4, fp.z>>4);
        if (n<=0) continue;
        /* 收集路口中心 */
        int c[256][2], nc=0;
        for (int i=0;i<n&&nc<256;i++){
            int t=L[i].type;
            if (t==1||t==2||t==7||t==11){ c[nc][0]=(L[i].bb0.x+L[i].bb1.x)/2; c[nc][1]=(L[i].bb0.z+L[i].bb1.z)/2; nc++; }
        }
        for (int i=0;i<nc;i++){
            /* 共线（19 格等距）*/
            for (int s=0;s<2;s++){
                int cnt=1;
                for (int k=1;k<20;k++){
                    int px=c[i][0]+(s?19*k:0), pz=c[i][1]+(s?0:19*k);
                    int f=0; for(int m=0;m<nc;m++) if(c[m][0]==px&&c[m][1]==pz){f=1;break;}
                    if(f) cnt++; else break;
                }
                if (cnt>bestline || (cnt==bestline && cnt>=4 && dd<bln_dd)){
                    bestline=cnt;
                    lnx=c[i][0]+(s?19.0*(cnt-1)/2:0); lnz=c[i][1]+(s?0:19.0*(cnt-1)/2);
                    bln_dd=dd;
                }
            }
            /* 2x2 方阵：保留最近的那个 */
            int h[4]={0,0,0,0};
            for(int m=0;m<nc;m++){
                if(c[m][0]==c[i][0]&&c[m][1]==c[i][1])h[0]=1;
                if(c[m][0]==c[i][0]+19&&c[m][1]==c[i][1])h[1]=1;
                if(c[m][0]==c[i][0]&&c[m][1]==c[i][1]+19)h[2]=1;
                if(c[m][0]==c[i][0]+19&&c[m][1]==c[i][1]+19)h[3]=1;
            }
            if(h[0]&&h[1]&&h[2]&&h[3] && dd<bsq_dd){
                bestsq=1; sqx=c[i][0]+9.5; sqz=c[i][1]+9.5; bsq_dd=dd; bdist=(int)sqrt((double)dd);
            }
        }
    }
    printf("  (4000格内: 真要塞 %d 个, 被堡垒挤掉 %d 个)\n", nfort, nbast);
    if (bestsq) printf("  有 2x2 四联方阵！/tp %.0f 200 %.0f   （要塞距原点 %d 格，同时最长共线 %d 联）\n", sqx, sqz, bdist, bestline);
    else if (bestline>=4) printf("  最长共线 %d 联  /tp %.0f 200 %.0f\n", bestline, lnx, lnz);
    else printf("  只找到 %d 联\n", bestline);
    return 0;
}
