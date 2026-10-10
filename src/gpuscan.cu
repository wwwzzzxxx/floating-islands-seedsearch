/* GPU 四联结构 base 扫描器（cubiomes quadbase 的 CUDA 移植）
 * 用法: ./gpuscan -start 0 -count N -hutr 128 [-monr 256] [-mode both|hut|count]
 *   -hutr 128  女巫小屋四联半径（128 走 Feature24 精确路径，其它走 Feature）
 *   -monr R    海底神殿四联半径（<=0 表示不检查神殿，纯枚举小屋 base）
 * 输出: 计数 + 命中种子列表
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <cuda_runtime.h>

#define CK(x) do { cudaError_t e=(x); if (e != cudaSuccess) { \
    fprintf(stderr, "CUDA 错误 %s @%d: %s\n", #x, __LINE__, cudaGetErrorString(e)); exit(1);} } while(0)

/* ---------- 移植自 cubiomes rng.h ---------- */
__device__ __forceinline__ int nextInt24(uint64_t &S)
{
    const uint64_t M = (1ULL << 48) - 1;
    const uint64_t K = 0x5deece66dULL;
    uint64_t a = M;
    uint64_t c = K * S;
    c += 11; a &= c;
    S = a;
    a = (uint64_t)((int64_t)a >> 17);
    c = 0xaaaaaaabULL * a;
    c = (uint64_t)((int64_t)c >> 36);
    return (int)a - (int)(c << 3) * 3;      /* a % 24 */
}

/* ---------- 移植自 cubiomes quadbase.h: getEnclosingRadius ---------- */
__device__ __forceinline__ float getEnclosingRadius(
    int x0, int z0, int x1, int z1, int x2, int z2, int x3, int z3,
    int ax, int ay, int az, int reg, int gap)
{
    x0 = (x0 << 4);
    z0 = (z0 << 4);
    x1 = ((reg + x1) << 4) + ax;
    z1 = ((reg + z1) << 4) + az;
    x2 = ((reg + x2) << 4) + ax;
    z2 = (z2 << 4);
    x3 = (x3 << 4);
    z3 = ((reg + z3) << 4) + az;

    int sqrad = 0x7fffffff;
    int cbx0 = (x1 > x2 ? x1 : x2) - gap;
    int cbz0 = (z1 > z3 ? z1 : z3) - gap;
    int cbx1 = (x0 < x3 ? x0 : x3) + gap;
    int cbz1 = (z0 < z2 ? z0 : z2) + gap;
    for (int z = cbz0; z <= cbz1; z++)
        for (int x = cbx0; x <= cbx1; x++) {
            int sq = 0, s;
            s = (x-x0)*(x-x0) + (z-z0)*(z-z0); if (s > sq) sq = s;
            s = (x-x1)*(x-x1) + (z-z1)*(z-z1); if (s > sq) sq = s;
            s = (x-x2)*(x-x2) + (z-z2)*(z-z2); if (s > sq) sq = s;
            s = (x-x3)*(x-x3) + (z-z3)*(z-z3); if (s > sq) sq = s;
            if (sq < sqrad) sqrad = sq;
        }
    return sqrad < 0x7fffffff ? sqrtf(sqrad + ay*ay/4.0f) : 0xffff;
}

/* ---------- 移植自 cubiomes quadbase.h: isQuadBaseFeature24（小屋, r=128 精确路径） ---------- */
__device__ __forceinline__ float quadHut24(uint64_t seed)
{
    const uint64_t SALT = 14357620ULL;                 /* Swamp_Hut, MC 1.18.2 */
    const uint64_t K = 0x5deece66dULL;
    seed += SALT;
    uint64_t s00 = seed;
    uint64_t s11 = 341873128712ULL + 132897987541ULL + seed;
    int x0, z0, x1, z1, x2, z2, x3, z3, x, z;

    s00 ^= K;
    x0 = nextInt24(s00); if (x0 < 20) return 0;
    z0 = nextInt24(s00); if (z0 < 20) return 0;
    s11 ^= K;
    x1 = nextInt24(s11); if (x1 > x0 - 20) return 0;
    z1 = nextInt24(s11); if (z1 > z0 - 20) return 0;
    x = x1 + 32 - x0; z = z1 + 32 - z0;
    if (x*x + z*z > 255) return 0;

    uint64_t s01 = 341873128712ULL + seed;
    uint64_t s10 = 132897987541ULL + seed;
    s01 ^= K;
    x2 = nextInt24(s01); if (x2 >= 4) return 0;
    z2 = nextInt24(s01); if (z2 < 20) return 0;
    s10 ^= K;
    x3 = nextInt24(s10); if (x3 < 20) return 0;
    z3 = nextInt24(s10); if (z3 >= 4) return 0;
    x = x2 + 32 - x3; z = z3 - z2 + 32;
    if (x*x + z*z > 255) return 0;

    float sqrad = getEnclosingRadius(x0,z0,x1,z1,x2,z2,x3,z3, 8,8,10, 32, 128);
    return sqrad < 128 ? sqrad : 0;
}

/* ---------- 移植自 cubiomes quadbase.h: isQuadBaseFeature（小屋, 任意半径） ---------- */
__device__ __forceinline__ float quadHutR(uint64_t seed, int radius)
{
    const uint64_t SALT = 14357620ULL;
    const uint64_t M = (1ULL << 48) - 1, K = 0x5deece66dULL, b = 0xb;
    const int R = 32, C = 24;
    seed += SALT;
    uint64_t s00 = seed;
    uint64_t s11 = 341873128712ULL + 132897987541ULL + seed;
    int x0,z0,x1,z1,x2,z2,x3,z3,x,z;
    int cd = radius / 8;
    int rm = R - (int)sqrtf((float)(cd*cd - (R-C+1)*(R-C+1)));
    uint64_t s;

    s = s00 ^ K;
    s = (s * K + b) & M; x0 = (int)(s >> 17) % C; if (x0 <= rm) return 0;
    s = (s * K + b) & M; z0 = (int)(s >> 17) % C; if (z0 <= rm) return 0;
    s = s11 ^ K;
    s = (s * K + b) & M; x1 = (int)(s >> 17) % C; if (x1 >= x0 - rm) return 0;
    s = (s * K + b) & M; z1 = (int)(s >> 17) % C; if (z1 >= z0 - rm) return 0;
    x = x1 + R - x0; z = z1 + R - z0;
    if (x*x + z*z > cd*cd) return 0;

    uint64_t s01 = 341873128712ULL + seed;
    uint64_t s10 = 132897987541ULL + seed;
    s = s01 ^ K;
    s = (s * K + b) & M; x2 = (int)(s >> 17) % C; if (x2 >= C - rm) return 0;
    s = (s * K + b) & M; z2 = (int)(s >> 17) % C; if (z2 <= rm) return 0;
    s = s10 ^ K;
    s = (s * K + b) & M; x3 = (int)(s >> 17) % C; if (x3 <= rm) return 0;
    s = (s * K + b) & M; z3 = (int)(s >> 17) % C; if (z3 >= C - rm) return 0;
    x = x2 + R - x3; z = z3 - z2 + R;
    if (x*x + z*z > cd*cd) return 0;

    float sqrad = getEnclosingRadius(x0,z0,x1,z1,x2,z2,x3,z3, 8,8,10, R, radius);
    return sqrad < radius ? sqrad : 0;
}

/* ---------- 移植自 cubiomes quadbase.h: isQuadBaseLarge（神殿, 任意半径） ---------- */
__device__ __forceinline__ float quadMon(uint64_t seed, int radius)
{
    const uint64_t SALT = 10387313ULL;                 /* Monument, MC 1.18.2 */
    const uint64_t M = (1ULL << 48) - 1, K = 0x5deece66dULL, b = 0xb;
    seed += SALT;
    uint64_t s00 = seed;
    uint64_t s01 = 341873128712ULL + seed;
    uint64_t s10 = 132897987541ULL + seed;
    uint64_t s11 = 341873128712ULL + 132897987541ULL + seed;
    const int R = 32, C = 27, ax = 58, ay = 23, az = 58;
    int rm = (int)(2 * R + ((ax < az ? ax : az) - 2*radius + 7) / 8);
    uint64_t s; int p;
    int x0,z0,x1,z1,x2,z2,x3,z3;

    s = s00 ^ K;
    s = (s * K + b) & M; p =  (int)(s >> 17) % C;
    s = (s * K + b) & M; p += (int)(s >> 17) % C; if (p <= rm) return 0;
    x0 = p;
    s = (s * K + b) & M; p =  (int)(s >> 17) % C;
    s = (s * K + b) & M; p += (int)(s >> 17) % C; if (p <= rm) return 0;
    z0 = p;

    s = s11 ^ K;
    s = (s * K + b) & M; p =  (int)(s >> 17) % C;
    s = (s * K + b) & M; p += (int)(s >> 17) % C; if (p > x0 - rm) return 0;
    x1 = p;
    s = (s * K + b) & M; p =  (int)(s >> 17) % C;
    s = (s * K + b) & M; p += (int)(s >> 17) % C; if (p > z0 - rm) return 0;
    z1 = p;

    s = ((x1-x0)>>1)*((x1-x0)>>1) + ((z1-z0)>>1)*((z1-z0)>>1);
    if (s > (uint64_t)4*radius*radius) return 0;

    s = s01 ^ K;
    s = (s * K + b) & M; p =  (int)(s >> 17) % C;
    s = (s * K + b) & M; p += (int)(s >> 17) % C; if (p > x0 - rm) return 0;
    x2 = p;
    s = (s * K + b) & M; p =  (int)(s >> 17) % C;
    s = (s * K + b) & M; p += (int)(s >> 17) % C; if (p <= rm) return 0;
    z2 = p;

    s = s10 ^ K;
    s = (s * K + b) & M; p =  (int)(s >> 17) % C;
    s = (s * K + b) & M; p += (int)(s >> 17) % C; if (p <= rm) return 0;
    x3 = p;
    s = (s * K + b) & M; p =  (int)(s >> 17) % C;
    s = (s * K + b) & M; p += (int)(s >> 17) % C; if (p > z0 - rm) return 0;
    z3 = p;

    float sqrad = getEnclosingRadius(x0>>1,z0>>1, x1>>1,z1>>1, x2>>1,z2>>1, x3>>1,z3>>1,
                                     ax, ay, az, R, radius);
    return sqrad < radius ? sqrad : 0;
}

__device__ __forceinline__ float quadHut(uint64_t seed, int radius)
{
    return (radius == 128) ? quadHut24(seed) : quadHutR(seed, radius);
}

/* ---------------- kernel ---------------- */
__global__ void scanKernel(uint64_t start, uint64_t count,
                           int hutR, int monR, int storeHut,
                           unsigned long long *hutHits, unsigned long long *bothHits,
                           unsigned int *nHut, unsigned int *nMon, unsigned int *nBoth,
                           unsigned int cap)
{
    uint64_t i = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t stride = (uint64_t)gridDim.x * blockDim.x;
    for (; i < count; i += stride) {
        uint64_t seed = start + i;
        float h = quadHut(seed, hutR);
        float m = (monR > 0) ? quadMon(seed, monR) : 0.0f;
        if (h != 0.0f) {
            if (storeHut) {
                unsigned int k = atomicAdd(nHut, 1u);
                if (k < cap) hutHits[k] = (unsigned long long)seed;
            } else atomicAdd(nHut, 1u);
        }
        if (m != 0.0f) atomicAdd(nMon, 1u);
        if (h != 0.0f && m != 0.0f) {
            unsigned int k = atomicAdd(nBoth, 1u);
            if (k < cap) bothHits[k] = (unsigned long long)seed;
        }
    }
}

int main(int argc, char **argv)
{
    uint64_t start = 0, count = 1000000;
    int hutR = 128, monR = 0, storeHut = 0;
    unsigned int cap = 8000000;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i],"-start") && i+1<argc) start = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i],"-count") && i+1<argc) count = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i],"-hutr")  && i+1<argc) hutR = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-monr")  && i+1<argc) monR = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-storehut")) storeHut = 1;
        else if (!strcmp(argv[i],"-cap") && i+1<argc) cap = atoi(argv[++i]);
    }
    int dev = 0; CK(cudaSetDevice(dev));
    cudaDeviceProp prop; CK(cudaGetDeviceProperties(&prop, dev));
    printf("GPU: %s  SM数=%d  扫描范围=[%llu,+%llu) 小屋r=%d 神殿r=%d\n",
           prop.name, prop.multiProcessorCount,
           (unsigned long long)start, (unsigned long long)count, hutR, monR);

    unsigned long long *dHut=NULL, *dBoth=NULL;
    unsigned int *dNHut=NULL, *dNMon=NULL, *dNBoth=NULL;
    CK(cudaMalloc(&dNHut, 4)); CK(cudaMalloc(&dNMon, 4)); CK(cudaMalloc(&dNBoth, 4));
    CK(cudaMemset(dNHut, 0, 4)); CK(cudaMemset(dNMon, 0, 4)); CK(cudaMemset(dNBoth, 0, 4));
    if (storeHut) { CK(cudaMalloc(&dHut, (size_t)cap*8)); CK(cudaMemset(dHut, 0, (size_t)cap*8)); }
    CK(cudaMalloc(&dBoth, (size_t)cap*8)); CK(cudaMemset(dBoth, 0, (size_t)cap*8));

    int blocks = prop.multiProcessorCount * 32;
    int threads = 256;
    cudaEvent_t t0, t1; CK(cudaEventCreate(&t0)); CK(cudaEventCreate(&t1));
    CK(cudaEventRecord(t0));
    scanKernel<<<blocks, threads>>>(start, count, hutR, monR, storeHut, dHut, dBoth, dNHut, dNMon, dNBoth, cap);
    CK(cudaDeviceSynchronize());
    CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
    float ms; CK(cudaEventElapsedTime(&ms, t0, t1));

    unsigned int nHut=0, nMon=0, nBoth=0;
    CK(cudaMemcpy(&nHut, dNHut, 4, cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(&nMon, dNMon, 4, cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(&nBoth, dNBoth, 4, cudaMemcpyDeviceToHost));
    printf("RESULT count=%llu hutR=%d monR=%d nHut=%u nMon=%u nBoth=%u\n",
           (unsigned long long)count, hutR, monR, nHut, nMon, nBoth);
    printf("耗时 %.3f 秒   吞吐 %.4g 种子/秒\n", ms/1000.0, count/(ms/1000.0));
    printf("小屋四联命中 : %u   1/%.4g\n", nHut, nHut?(double)count/nHut:0.0);
    if (monR > 0) {
        printf("神殿四联命中 : %u   1/%.4g\n", nMon, nMon?(double)count/nMon:0.0);
        printf("两者都四联   : %u   1/%.4g\n", nBoth, nBoth?(double)count/nBoth:0.0);
    }

    if (nBoth && nBoth <= 200000) {
        unsigned long long *h = (unsigned long long*)malloc((size_t)nBoth*8);
        CK(cudaMemcpy(h, dBoth, (size_t)nBoth*8, cudaMemcpyDeviceToHost));
        printf("--- 两者都四联的种子 ---\n");
        for (unsigned int i = 0; i < nBoth; i++) printf("%lld\n", (long long)(int64_t)h[i]);
        free(h);
    } else if (nBoth) printf("(命中 %u 个，过多不打印)\n", nBoth);

    if (storeHut && nHut && nHut <= 2000000) {
        unsigned long long *h = (unsigned long long*)malloc((size_t)nHut*8);
        CK(cudaMemcpy(h, dHut, (size_t)nHut*8, cudaMemcpyDeviceToHost));
        FILE *f = fopen("hutbases.txt", "w");
        for (unsigned int i = 0; i < nHut; i++) fprintf(f, "%lld\n", (long long)(int64_t)h[i]);
        fclose(f); free(h);
        printf("小屋 base 列表已写入 hutbases.txt（%u 个）\n", nHut);
    } else if (storeHut) printf("(小屋 base %u 个，超过打印/落盘上限)\n", nHut);
    return 0;
}
