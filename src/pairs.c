/* pairs.c — 导出所有"与小屋四联兼容的神殿四联"及其单人覆盖率
 * 用法: pairs [线程数] [神殿半径] [偏移范围]
 * 输出: monpairs.txt，每行 "base oi oj coverage"
 *       coverage = getOptimalAfk 算出的、单人挂机点能覆盖的神殿刷怪方块数（上限 13456）
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include "generator.h"
#include "finders.h"
#include "quadbase.h"
#define MCV MC_1_18_2

static uint64_t *bases; static long nbases;
static int P_OFF = 16, P_MON_R = 256;
static StructureConfig mc;
static int NEXT;
static pthread_mutex_t MTX = PTHREAD_MUTEX_INITIALIZER;
static FILE *PF;
static pthread_mutex_t RTX = PTHREAD_MUTEX_INITIALIZER;
static long nquad;

static void *worker(void *arg)
{
    (void)arg;
    static const int rc[4][2] = {{0,0},{0,1},{1,0},{1,1}};
    for (;;) {
        pthread_mutex_lock(&MTX); long idx = NEXT++; pthread_mutex_unlock(&MTX);
        if (idx >= nbases) break;
        uint64_t base = bases[idx] & 0xFFFFFFFFFFFFULL;
        for (int oi = -P_OFF; oi <= P_OFF; oi++)
        for (int oj = -P_OFF; oj <= P_OFF; oj++) {
            uint64_t sp = moveStructure(base, -oi, -oj);
            if (isQuadBase(mc, sp, P_MON_R) == 0.0f) continue;
            Pos p[4];
            for (int i = 0; i < 4; i++)
                p[i] = getStructurePos(Monument, MCV, base, oi+rc[i][0], oj+rc[i][1], &p[i]) ? p[i] : p[i];
            int sp2;
            getOptimalAfk(p, 58, 0, 58, &sp2);
            pthread_mutex_lock(&RTX);
            if (PF) fprintf(PF, "%llu %d %d %d\n",
                            (unsigned long long)base, oi, oj, sp2);
            nquad++;
            pthread_mutex_unlock(&RTX);
        }
    }
    return NULL;
}

int main(int argc, char **argv)
{
    int NT = argc > 1 ? atoi(argv[1]) : 24;
    int R  = argc > 2 ? atoi(argv[2]) : 256;
    if (argc > 3) P_OFF = atoi(argv[3]);

    FILE *f = fopen("hutbases_local.txt", "r");
    if (!f) { perror("hutbases_local.txt"); return 1; }
    bases = malloc(sizeof(uint64_t) * 800000); nbases = 0;
    while (nbases < 800000 && fscanf(f, "%llu", (unsigned long long*)&bases[nbases]) == 1) nbases++;
    fclose(f);
    getStructureConfig(Monument, MCV, &mc);

    PF = fopen("monpairs.txt", "w");
    if (!PF) { perror("monpairs.txt"); return 1; }

    printf("小屋基 %ld 个 x %d 偏移，神殿半径<=%d\n",
           nbases, (2*P_OFF+1)*(2*P_OFF+1), R);
    pthread_t th[256];
    for (int i = 0; i < NT; i++) pthread_create(&th[i], NULL, worker, NULL);
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    fclose(PF);
    printf("共 %ld 个兼容的神殿四联，已写入 monpairs.txt\n", nquad);
    return 0;
}
