#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>
static double tnow(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

#ifdef DOUBLE
typedef double R;
#else
typedef float R;
#endif

__constant__ int DI, DJ, NT, NB, MAXW;
__constant__ int N, LT, MT, NTT;
__constant__ R CMU, C2, CMU1, CMU2, EE, CK, HINUM, SMNUM, ALVIS;

typedef struct {
  int i0, i1, j0, j1, k0, k1;
} Box;

enum {
  U_,
  V_,
  W_,
  P_,
  PP_,
  DK_,
  DE_,
  TM_,
  VISE_,
  DEN_,
  DU_,
  DV_,
  DW_,
  GEN_,
  MCF_,
  AP_,
  SU_,
  SP_,
  SUK_,
  SPK_,
  AE_,
  AW_,
  AN_,
  AS_,
  AT_,
  AB_,
  X_,
  Y_,
  Z_,
  TJO_,
  CX_,
  CY_,
  CZ_,
  EX_,
  EY_,
  EZ_,
  SX_,
  SY_,
  SZ_,
  TXXE_,
  TXXW_,
  TYYN_,
  TYYS_,
  TZZT_,
  TZZB_,
  TYXE_,
  TYXW_,
  TYZT_,
  TYZB_,
  TXYN_,
  TXYS_,
  TXZT_,
  TXZB_,
  TZXE_,
  TZXW_,
  TZYN_,
  TZYS_,
  F1_,
  FOLD_,
  NA_
};

typedef struct {
  R *a[NA_];
  int *ibc, *jbc, *kbc, *ity, *nw;
  R *yn, *yn1, *sx, *sy, *sz, *anv, *anw, *ypl, *tau;
  R *visc, *flowin, *cfac, *err, *alf, *c1;
  int *act, *flg;
} S;
__constant__ S G;

#define IDX(i, j, k) ((i) - 1 + DI * (((j) - 1) + DJ * ((k) - 1)))
#define A(n, i, j, k) G.a[n##_][(size_t)c * NT + IDX(i, j, k)]
#define CELL(B)                                                                                    \
  int tid = blockIdx.x * blockDim.x + threadIdx.x;                                                 \
  int ni = B.i1 - B.i0 + 1, nj = B.j1 - B.j0 + 1, nk = B.k1 - B.k0 + 1, nn = ni * nj * nk;         \
  if (ni <= 0 || nj <= 0 || nk <= 0 || tid >= nn * NB)                                             \
    return;                                                                                        \
  int c = tid / nn, r = tid % nn;                                                                  \
  int i = B.i0 + r % ni, j = B.j0 + (r / ni) % nj, k = B.k0 + r / (ni * nj)
#define DECODE(B)                                                                                  \
  CELL(B);                                                                                         \
  if (!G.act[c])                                                                                   \
  return

__device__ __forceinline__ R *field(int which, int c) {
  return G.a[which] + (size_t)c * NT;
}
#define F(i, j, k) f[IDX(i, j, k)]

#ifdef DOUBLE
__device__ void atomicMaxD(double *a, double v) {
  unsigned long long *p = (unsigned long long *)a, old = *p, as;
  do {
    as = old;
    if (__longlong_as_double(as) >= v)
      return;
    old = atomicCAS(p, as, __double_as_longlong(v));
  } while (old != as);
}
#else
__device__ void atomicMaxD(float *a, float v) {
  atomicMax((int *)a, __float_as_int(v));
}
#endif

__global__ void k_diffusivity(Box B, R sigf) {
  DECODE(B);
  A(F1, i, j, k) = A(VISE, i, j, k) / sigf;
}

__global__ void __launch_bounds__(256) k_coef(Box B, int fv) {
  DECODE(B);
  R *f = field(fv, c);
  R GAE = 0.5 * (A(F1, i + 1, j, k) + A(F1, i, j, k)),
    GAW = 0.5 * (A(F1, i - 1, j, k) + A(F1, i, j, k));
  R GAN = 0.5 * (A(F1, i, j + 1, k) + A(F1, i, j, k)),
    GAS = 0.5 * (A(F1, i, j - 1, k) + A(F1, i, j, k));
  R GAT = 0.5 * (A(F1, i, j, k + 1) + A(F1, i, j, k)),
    GAB = 0.5 * (A(F1, i, j, k - 1) + A(F1, i, j, k));
  R DENE = 0.5 * (A(DEN, i + 1, j, k) + A(DEN, i, j, k)),
    DENW = 0.5 * (A(DEN, i - 1, j, k) + A(DEN, i, j, k));
  R DENN = 0.5 * (A(DEN, i, j + 1, k) + A(DEN, i, j, k)),
    DENS = 0.5 * (A(DEN, i, j - 1, k) + A(DEN, i, j, k));
  R DENT = 0.5 * (A(DEN, i, j, k + 1) + A(DEN, i, j, k)),
    DENB = 0.5 * (A(DEN, i, j, k - 1) + A(DEN, i, j, k));
#define AV(v, di, dj, dk) (0.5 * (A(v, i + di, j + dj, k + dk) + A(v, i, j, k)))
  R UE = AV(U, 1, 0, 0), UW = AV(U, -1, 0, 0), UN = AV(U, 0, 1, 0), US = AV(U, 0, -1, 0),
    UT = AV(U, 0, 0, 1), UB = AV(U, 0, 0, -1);
  R VE = AV(V, 1, 0, 0), VW = AV(V, -1, 0, 0), VN = AV(V, 0, 1, 0), VS = AV(V, 0, -1, 0),
    VT = AV(V, 0, 0, 1), VB = AV(V, 0, 0, -1);
  R WE = AV(W, 1, 0, 0), WW = AV(W, -1, 0, 0), WN = AV(W, 0, 1, 0), WS = AV(W, 0, -1, 0),
    WT = AV(W, 0, 0, 1), WB = AV(W, 0, 0, -1);
  R cx = A(CX, i, j, k), cy = A(CY, i, j, k), cz = A(CZ, i, j, k), ex = A(EX, i, j, k),
    ey = A(EY, i, j, k), ez = A(EZ, i, j, k);
  R sx = A(SX, i, j, k), sy = A(SY, i, j, k), sz = A(SZ, i, j, k), tj = A(TJO, i, j, k);
  R CE = DENE * (UE * cx + VE * cy + WE * cz), CW = DENW * (UW * cx + VW * cy + WW * cz);
  R CN = DENN * (UN * ex + VN * ey + WN * ez), CS = DENS * (US * ex + VS * ey + WS * ez);
  R CT = DENT * (UT * sx + VT * sy + WT * sz), CB = DENB * (UB * sx + VB * sy + WB * sz);
  R TXXE = A(TXXE, i, j, k), TXXW = A(TXXW, i, j, k), TYYN = A(TYYN, i, j, k),
    TYYS = A(TYYS, i, j, k), TZZT = A(TZZT, i, j, k), TZZB = A(TZZB, i, j, k);
  R TYXE = A(TYXE, i, j, k), TYXW = A(TYXW, i, j, k), TYZT = A(TYZT, i, j, k),
    TYZB = A(TYZB, i, j, k), TXYN = A(TXYN, i, j, k), TXYS = A(TXYS, i, j, k);
  R TXZT = A(TXZT, i, j, k), TXZB = A(TXZB, i, j, k), TZXE = A(TZXE, i, j, k),
    TZXW = A(TZXW, i, j, k), TZYN = A(TZYN, i, j, k), TZYS = A(TZYS, i, j, k);
  R DDE = GAE * TXXE + GAN * TXYN - GAS * TXYS + GAT * TXZT - GAB * TXZB;
  R DDW = GAW * TXXW + GAS * TXYS - GAN * TXYN - GAT * TXZT + GAB * TXZB;
  R DDN = GAN * TYYN + GAE * TYXE - GAW * TYXW + GAT * TYZT - GAB * TYZB;
  R DDS = GAS * TYYS + GAW * TYXW - GAE * TYXE - GAT * TYZT + GAB * TYZB;
  R DDT = GAT * TZZT + GAE * TZXE - GAW * TZXW + GAN * TZYN - GAS * TZYS;
  R DDB = GAB * TZZB + GAW * TZXW - GAE * TZXE - GAN * TZYN + GAS * TZYS;
  R CPO = fabs(CE - CW + CN - CS + CT - CB), cf = G.cfac[c];
  A(AE, i, j, k) = (fmax(fabs(cf * CE), DDE) - 0.5 * CE) * tj;
  A(AW, i, j, k) = (fmax(fabs(cf * CW), DDW) + 0.5 * CW) * tj;
  A(AN, i, j, k) = (fmax(fabs(cf * CN), DDN) - 0.5 * CN) * tj;
  A(AS, i, j, k) = (fmax(fabs(cf * CS), DDS) + 0.5 * CS) * tj;
  A(AT, i, j, k) = (fmax(fabs(cf * CT), DDT) - 0.5 * CT) * tj;
  A(AB, i, j, k) = (fmax(fabs(cf * CB), DDB) + 0.5 * CB) * tj;
  R DDNE = GAE * TYXE + GAN * TXYN, DDSE = -GAE * TYXE - GAS * TXYS,
    DDNW = -GAW * TYXW - GAN * TXYN, DDSW = GAW * TYXW + GAS * TXYS;
  R DDET = GAE * TZXE + GAT * TXZT, DDEB = -GAE * TZXE - GAB * TXZB,
    DDWT = -GAW * TZXW - GAT * TXZT, DDWB = GAW * TZXW + GAB * TXZB;
  R DDNT = GAN * TZYN + GAT * TYZT, DDNB = -GAN * TZYN - GAB * TYZB,
    DDST = -GAS * TZYS - GAT * TYZT, DDSB = GAS * TZYS + GAB * TYZB;
#define CROSS(Q)                                                                                   \
  (CPO * Q(i, j, k) + DDNE * Q(i + 1, j + 1, k) + DDSE * Q(i + 1, j - 1, k) +                      \
   DDNW * Q(i - 1, j + 1, k) + DDSW * Q(i - 1, j - 1, k) + DDET * Q(i + 1, j, k + 1) +             \
   DDEB * Q(i + 1, j, k - 1) + DDWT * Q(i - 1, j, k + 1) + DDWB * Q(i - 1, j, k - 1) +             \
   DDNT * Q(i, j + 1, k + 1) + DDNB * Q(i, j + 1, k - 1) + DDST * Q(i, j - 1, k + 1) +             \
   DDSB * Q(i, j - 1, k - 1))
#define QF(a, b, d) F(a, b, d)
#define QV(a, b, d) A(V, a, b, d)
#define QW(a, b, d) A(W, a, b, d)
#define QK(a, b, d) A(DK, a, b, d)
  A(SU, i, j, k) = CROSS(QF);
  A(DV, i, j, k) = CROSS(QV);
  A(DW, i, j, k) = CROSS(QW);
  A(PP, i, j, k) = CROSS(QK);
  A(SUK, i, j, k) = CPO;
  A(SPK, i, j, k) = -CPO;
  R UX = (UE - UW) * cx + (UN - US) * ex + (UT - UB) * sx,
    UY = (UE - UW) * cy + (UN - US) * ey + (UT - UB) * sy,
    UZ = (UE - UW) * cz + (UN - US) * ez + (UT - UB) * sz;
  R VX = (VE - VW) * cx + (VN - VS) * ex + (VT - VB) * sx,
    VY = (VE - VW) * cy + (VN - VS) * ey + (VT - VB) * sy,
    VZ = (VE - VW) * cz + (VN - VS) * ez + (VT - VB) * sz;
  R WX = (WE - WW) * cx + (WN - WS) * ex + (WT - WB) * sx,
    WY = (WE - WW) * cy + (WN - WS) * ey + (WT - WB) * sy,
    WZ = (WE - WW) * cz + (WN - WS) * ez + (WT - WB) * sz;
  A(GEN, i, j, k) =
      A(VISE, i, j, k) * ((UY + VX) * (UY + VX) + (VZ + WY) * (VZ + WY) + (WX + UZ) * (WX + UZ) +
                          ((G.flg[c] & 1) ? (R)0 : 2 * (UX * UX + VY * VY + WZ * WZ)));
}

__global__ void __launch_bounds__(256) k_src_u(Box B) {
  DECODE(B);
#define F1A(a, b, d) A(F1, a, b, d)
  R GAE = 0.5 * (F1A(i + 1, j, k) + F1A(i, j, k)), GAW = 0.5 * (F1A(i - 1, j, k) + F1A(i, j, k));
  R GAN = 0.5 * (F1A(i, j + 1, k) + F1A(i, j, k)), GAS = 0.5 * (F1A(i, j - 1, k) + F1A(i, j, k));
  R GAT = 0.5 * (F1A(i, j, k + 1) + F1A(i, j, k)), GAB = 0.5 * (F1A(i, j, k - 1) + F1A(i, j, k));
#define PA(a, b, d) A(P, a, b, d)
  R PCXI = (PA(i + 1, j, k) + PA(i + 1, j + 1, k) + PA(i + 1, j, k + 1) + PA(i + 1, j + 1, k + 1) -
            PA(i, j, k) - PA(i, j + 1, k) - PA(i, j, k + 1) - PA(i, j + 1, k + 1)) *
           0.25;
  R PEDA = (PA(i, j + 1, k) + PA(i + 1, j + 1, k) + PA(i, j + 1, k + 1) + PA(i + 1, j + 1, k + 1) -
            PA(i, j, k) - PA(i + 1, j, k) - PA(i, j, k + 1) - PA(i + 1, j, k + 1)) *
           0.25;
  R PSCI = (PA(i, j, k + 1) + PA(i + 1, j, k + 1) + PA(i, j + 1, k + 1) + PA(i + 1, j + 1, k + 1) -
            PA(i, j, k) - PA(i + 1, j, k) - PA(i, j + 1, k) - PA(i + 1, j + 1, k)) *
           0.25;
#define Q(q, a, b, d) q##p[(a) + (b) * DI + (d) * DI * DJ]
#define GRADS(q)                                                                                   \
  const R *q##p = &A(q, i, j, k);                                                                  \
  R q##CE = Q(q, 1, 0, 0) - Q(q, 0, 0, 0), q##CW = Q(q, 0, 0, 0) - Q(q, -1, 0, 0);                 \
  R q##CN = (Q(q, 1, 1, 0) + Q(q, 1, 0, 0) - Q(q, -1, 1, 0) - Q(q, -1, 0, 0)) * 0.25;              \
  R q##CS = (Q(q, 1, 0, 0) + Q(q, 1, -1, 0) - Q(q, -1, 0, 0) - Q(q, -1, -1, 0)) * 0.25;            \
  R q##CT = (Q(q, 1, 0, 1) + Q(q, 1, 0, 0) - Q(q, -1, 0, 1) - Q(q, -1, 0, 0)) * 0.25;              \
  R q##CB = (Q(q, 1, 0, 0) + Q(q, 1, 0, -1) - Q(q, -1, 0, 0) - Q(q, -1, 0, -1)) * 0.25;            \
  R q##EE = (Q(q, 1, 1, 0) + Q(q, 0, 1, 0) - Q(q, 1, -1, 0) - Q(q, 0, -1, 0)) * 0.25;              \
  R q##EW = (Q(q, 0, 1, 0) + Q(q, -1, 1, 0) - Q(q, 0, -1, 0) - Q(q, -1, -1, 0)) * 0.25;            \
  R q##EN = Q(q, 0, 1, 0) - Q(q, 0, 0, 0), q##ES = Q(q, 0, 0, 0) - Q(q, 0, -1, 0);                 \
  R q##ET = (Q(q, 0, 1, 1) + Q(q, 0, 1, 0) - Q(q, 0, -1, 1) - Q(q, 0, -1, 0)) * 0.25;              \
  R q##EB = (Q(q, 0, 1, 0) + Q(q, 0, 1, -1) - Q(q, 0, -1, 0) - Q(q, 0, -1, -1)) * 0.25;            \
  R q##SE = (Q(q, 1, 0, 1) + Q(q, 0, 0, 1) - Q(q, 1, 0, -1) - Q(q, 0, 0, -1)) * 0.25;              \
  R q##SW = (Q(q, 0, 0, 1) + Q(q, -1, 0, 1) - Q(q, 0, 0, -1) - Q(q, -1, 0, -1)) * 0.25;            \
  R q##SN = (Q(q, 0, 1, 1) + Q(q, 0, 0, 1) - Q(q, 0, 1, -1) - Q(q, 0, 0, -1)) * 0.25;              \
  R q##SS = (Q(q, 0, 0, 1) + Q(q, 0, -1, 1) - Q(q, 0, 0, -1) - Q(q, 0, -1, -1)) * 0.25;            \
  R q##ST = Q(q, 0, 0, 1) - Q(q, 0, 0, 0), q##SB = Q(q, 0, 0, 0) - Q(q, 0, 0, -1);
  GRADS(U);
  GRADS(V);
  GRADS(W);
  R cx = A(CX, i, j, k), cy = A(CY, i, j, k), cz = A(CZ, i, j, k), ex = A(EX, i, j, k),
    ey = A(EY, i, j, k), ez = A(EZ, i, j, k);
  R sx = A(SX, i, j, k), sy = A(SY, i, j, k), sz = A(SZ, i, j, k);
  R su = A(SU, i, j, k) - PCXI * cx - PEDA * ex - PSCI * sx;
  R dv = A(DV, i, j, k) - PCXI * cy - PEDA * ey - PSCI * sy;
  R dw = A(DW, i, j, k) - PCXI * cz - PEDA * ez - PSCI * sz;
#define FACE(c0)                                                                                   \
  R c0##E = (A(c0, i + 1, j, k) + A(c0, i, j, k)) * 0.5,                                           \
    c0##W = (A(c0, i - 1, j, k) + A(c0, i, j, k)) * 0.5,                                           \
    c0##N = (A(c0, i, j + 1, k) + A(c0, i, j, k)) * 0.5,                                           \
    c0##S = (A(c0, i, j, k) + A(c0, i, j - 1, k)) * 0.5,                                           \
    c0##T = (A(c0, i, j, k + 1) + A(c0, i, j, k)) * 0.5,                                           \
    c0##B = (A(c0, i, j, k) + A(c0, i, j, k - 1)) * 0.5;
#define SOC(q, Ca, Ea, Sa, x, e, s)                                                                \
  (x * (GAE * (q##CE * Ca##E + q##EE * Ea##E + q##SE * Sa##E) -                                    \
        GAW * (q##CW * Ca##W + q##EW * Ea##W + q##SW * Sa##W)) +                                   \
   e * (GAN * (q##CN * Ca##N + q##EN * Ea##N + q##SN * Sa##N) -                                    \
        GAS * (q##CS * Ca##S + q##ES * Ea##S + q##SS * Sa##S)) +                                   \
   s * (GAT * (q##CT * Ca##T + q##ET * Ea##T + q##ST * Sa##T) -                                    \
        GAB * (q##CB * Ca##B + q##EB * Ea##B + q##SB * Sa##B)))
  FACE(CX);
  FACE(EX);
  FACE(SX);
  FACE(CY);
  FACE(EY);
  FACE(SY);
  FACE(CZ);
  FACE(EZ);
  FACE(SZ);
  su += SOC(U, CX, EX, SX, cx, ex, sx) + SOC(V, CX, EX, SX, cy, ey, sy) +
        SOC(W, CX, EX, SX, cz, ez, sz);
  dv += SOC(U, CY, EY, SY, cx, ex, sx) + SOC(V, CY, EY, SY, cy, ey, sy) +
        SOC(W, CY, EY, SY, cz, ez, sz);
  dw += SOC(U, CZ, EZ, SZ, cx, ex, sx) + SOC(V, CZ, EZ, SZ, cy, ey, sy) +
        SOC(W, CZ, EZ, SZ, cz, ez, sz);
  R tj = A(TJO, i, j, k);
  A(SU, i, j, k) = su * tj;
  A(DV, i, j, k) = dv;
  A(DW, i, j, k) = dw;
  A(SP, i, j, k) = A(SPK, i, j, k) * tj;
}
__global__ void k_src_vw(Box B, int w) {
  DECODE(B);
  A(SU, i, j, k) = (w ? A(DW, i, j, k) : A(DV, i, j, k)) * A(TJO, i, j, k);
}
__global__ void k_src_k(Box B) {
  DECODE(B);
  R tj = A(TJO, i, j, k), d = A(DEN, i, j, k);
  A(SU, i, j, k) = (A(GEN, i, j, k) + A(PP, i, j, k)) * tj;
  A(SP, i, j, k) = (A(SPK, i, j, k) - CMU * d * d * A(DK, i, j, k) / A(VISE, i, j, k)) * tj;
}
__global__ void k_src_e(Box B) {
  DECODE(B);
  R tj = A(TJO, i, j, k), d = A(DEN, i, j, k);
  A(SU, i, j, k) = (A(SU, i, j, k) +
                    G.c1[c] * CMU * A(GEN, i, j, k) * d * d * A(DK, i, j, k) / A(VISE, i, j, k)) *
                   tj;
  A(SP, i, j, k) = (A(SPK, i, j, k) - C2 * d * A(DE, i, j, k) / (A(DK, i, j, k) + SMNUM)) * tj;
}

__device__ void wall_function(int ie, int c, int w, int i, int j, int k, int ii, int jj, int kk,
                              int i2, int j2, int k2, R *f, R *an) {
  size_t o = (size_t)c * MAXW + w;
  R YP = G.yn[o], YP1 = G.yn1[o], DK1 = A(DK, i, j, k), DEN1 = A(DEN, i, j, k),
    TJ = A(TJO, i, j, k), visc = G.visc[c];
  R *SU1 = &A(SU, i, j, k), *SP1 = &A(SP, i, j, k);
  R F1v = f[IDX(ii, jj, kk)];
  if (ie == 1) {
    R SQ = sqrt(DK1), YPL = DEN1 * SQ * CMU1 * YP / visc, TM;
    G.ypl[o] = YPL;
    TM = YPL <= 11.63 ? visc / YP : DEN1 * CMU1 * SQ * CK / log(EE * YPL);
    G.tau[o] = -TM;
    R PTA = TJ * G.sx[o] * TM / YP1;
    *SP1 -= PTA;
    *SU1 += PTA * F1v;
    R P1 = sqrt(1 - G.sx[o] * G.sx[o]);
    G.anv[o] = *an * sqrt(1 - G.sy[o] * G.sy[o]);
    G.anw[o] = *an * sqrt(1 - G.sz[o] * G.sz[o]);
    *an = P1 * *an;
  } else if (ie == 2 || ie == 3) {
    R TM = -G.tau[o], s = ie == 2 ? G.sy[o] : G.sz[o], PTA = TJ * s * TM / YP1;
    *SP1 -= PTA;
    *SU1 += PTA * F1v;
    *an = ie == 2 ? G.anv[o] : G.anw[o];
  } else if (ie == 5) {
    R SQ = sqrt(DK1), YPL = G.ypl[o];
    R DITM =
        YPL <= 11.63 ? DEN1 * CMU2 * SQ * YPL / YP : DEN1 * CMU2 * SQ * log(EE * YPL) / (CK * YP);
    R du = A(U, i, j, k) - A(U, ii, jj, kk), dv = A(V, i, j, k) - A(V, ii, jj, kk),
      dw = A(W, i, j, k) - A(W, ii, jj, kk),
      P1 = (G.flg[c] & 2) ? (G.sx[o] * du) * (G.sx[o] * du) + (G.sy[o] * dv) * (G.sy[o] * dv) +
                                (G.sz[o] * dw) * (G.sz[o] * dw)
                          : du * du + dv * dv + dw * dw;
    du = A(U, i2, j2, k2) - A(U, ii, jj, kk);
    dv = A(V, i2, j2, k2) - A(V, ii, jj, kk);
    dw = A(W, i2, j2, k2) - A(W, ii, jj, kk);
    R P2 = du * du + dv * dv + dw * dw;
    R YP2 = 2 * YP1 - YP, G1 = A(GEN, i, j, k), VI = A(VISE, i, j, k),
      GENR = G1 - VI * P2 / YP2 / YP2;
    if (GENR <= 0)
      GENR = 0;
    R SPK1 = A(SPK, i, j, k), tau = G.tau[o];
    *SU1 = *SU1 + TJ * (P1 * tau * tau / VI + GENR - G1);
    *SP1 = TJ * (SPK1 - DITM);
    *an = 0;
  } else if (ie == 6) {
    R TERM = CMU2 / (CK * YP);
    *SU1 = HINUM * TERM * pow(DK1, 1.5);
    *SP1 = -HINUM;
  }
}
__global__ void k_wall(int ie, int fv) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= NB * MAXW)
    return;
  int c = tid / MAXW, w0 = tid % MAXW;
  if (!G.act[c] || w0 >= G.nw[c])
    return;
  size_t base = (size_t)c * MAXW;
  if (w0 > 0 && G.ibc[base + w0 - 1] == G.ibc[base + w0] &&
      G.jbc[base + w0 - 1] == G.jbc[base + w0] && G.kbc[base + w0 - 1] == G.kbc[base + w0])
    return;
  int lk[7] = {0, AN_, AS_, AE_, AW_, AT_, AB_};
  R *f = field(fv, c);
  for (int w = w0; w < G.nw[c] && G.ibc[base + w] == G.ibc[base + w0] &&
                   G.jbc[base + w] == G.jbc[base + w0] && G.kbc[base + w] == G.kbc[base + w0];
       w++) {
    size_t o = base + w;
    int i = G.ibc[o], j = G.jbc[o], k = G.kbc[o], t = G.ity[o];
    int di = (t == 3) - (t == 4), dj = (t == 1) - (t == 2), dk = (t == 5) - (t == 6);
    int q = lk[t];
    wall_function(ie, c, w, i, j, k, i + di, j + dj, k + dk, i - di, j - dj, k - dk, f,
                  &G.a[q][(size_t)c * NT + IDX(i, j, k)]);
  }
}
__global__ void k_outlet_links(Box B) {
  DECODE(B);
  if (G.flg[c] & 16)
    return;
  A(AE, i, j, k) = 0;
}
__global__ void k_assemble(Box B, int fv, int ie) {
  DECODE(B);
  R *f = field(fv, c);
  R alf = G.alf[c * 8 + ie];
  A(F1, i, j, k) = F(i, j, k);
  R ANAB = A(AE, i, j, k) + A(AW, i, j, k) + A(AN, i, j, k) + A(AS, i, j, k) + A(AT, i, j, k) +
           A(AB, i, j, k);
  R ap = ANAB - A(SP, i, j, k), pduv = 1;
  if (A(MCF, i, j, k) >= 1) {
    ap = alf;
    A(AN, i, j, k) = A(AS, i, j, k) = A(AE, i, j, k) = A(AW, i, j, k) = A(AT, i, j, k) =
        A(AB, i, j, k) = 0;
    A(SU, i, j, k) = F(i, j, k);
    pduv = 0;
  }
  R P1 = 1.2 * ap;
  ap = ap / alf;
  A(AP, i, j, k) = ap;
  A(SU, i, j, k) += pduv * (1 - alf) * ap * F(i, j, k);
  if (ie == 1)
    A(DU, i, j, k) = A(TJO, i, j, k) * pduv / (P1 - ANAB);
}
__device__ void wmax(int c, int slot, R v) {
  unsigned m = __activemask();
  int c0 = __shfl_sync(m, c, __ffs(m) - 1);
  int same = __all_sync(m, c == c0);
  if (same) {
    for (int o = 16; o > 0; o >>= 1) {
      R w = __shfl_xor_sync(m, v, o);
      v = fmax(v, w);
    }
    if ((threadIdx.x & 31) == __ffs(m) - 1)
      atomicMaxD(&G.err[c * 8 + slot], v);
  } else
    atomicMaxD(&G.err[c * 8 + slot], v);
}
__global__ void k_line(Box B, int fv, int color) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x, ni = B.i1 - B.i0 + 1, nj = B.j1 - B.j0 + 1,
      nn = ni * nj;
  if (tid >= nn * NB)
    return;
  int c = tid / nn, r = tid % nn, i = B.i0 + r % ni, j = B.j0 + r / ni;
  if (!G.act[c] || ((i + j) & 1) != color)
    return;
  R *f = field(fv, c);
  R a[16], cc[16];
  int ks = B.k0, nt = B.k1;
  a[ks - 1] = 0;
  cc[ks - 1] = F(i, j, ks - 1);
  for (int k = ks; k <= nt; k++) {
    R At = A(AT, i, j, k), Bb = A(AB, i, j, k),
      Cc = A(AE, i, j, k) * F(i + 1, j, k) + A(AN, i, j, k) * F(i, j + 1, k) + A(SU, i, j, k) +
           A(AW, i, j, k) * F(i - 1, j, k) + A(AS, i, j, k) * F(i, j - 1, k);
    R term = 1.0 / (A(AP, i, j, k) - Bb * a[k - 1] + 1e-30);
    a[k] = At * term;
    cc[k] = (Cc + Bb * cc[k - 1]) * term;
  }
  for (int k = nt; k >= ks; k--)
    F(i, j, k) = a[k] * F(i, j, k + 1) + cc[k];
}
__global__ void __launch_bounds__(1024) k_lines(Box B, int fv, int nsw) {
  extern __shared__ R sh[];
  int c = blockIdx.x;
  if (c >= NB || !G.act[c])
    return;
  int ni = B.i1 - B.i0 + 1, nj = B.j1 - B.j0 + 1, nk = B.k1 - B.k0 + 1, np = ni * nj, fi = ni + 2,
      fj = nj + 2, fp = fi * fj;
  R *cf = sh, *fs = sh + 8 * nk * np;
  R *f = field(fv, c);
#define CIX(q, ii, jj, kk) cf[((q) * nk + (kk)) * np + (jj) * ni + (ii)]
#define FIX(ii, jj, kk) fs[(kk) * fp + (jj) * fi + (ii)]
  int ids[8] = {AT_, AB_, AE_, AN_, SU_, AW_, AS_, AP_};
  enum { cT, cB, cE, cN, cSU, cW, cS, cP };
  for (int t = threadIdx.x; t < 8 * nk * np; t += blockDim.x) {
    int q = t / (nk * np), r = t % (nk * np), kk = r / np, jj = (r % np) / ni, ii = r % ni;
    cf[t] = G.a[ids[q]][(size_t)c * NT + IDX(B.i0 + ii, B.j0 + jj, B.k0 + kk)];
  }
  for (int t = threadIdx.x; t < (nk + 2) * fp; t += blockDim.x) {
    int kk = t / fp, r = t % fp, jj = r / fi, ii = r % fi;
    fs[t] = f[IDX(B.i0 - 1 + ii, B.j0 - 1 + jj, B.k0 - 1 + kk)];
  }
  __syncthreads();
  for (int s = 0; s < nsw; s++)
    for (int color = 0; color < 2; color++) {
      int nh = (ni + 1) / 2;
      for (int t = threadIdx.x; t < nh * nj; t += blockDim.x) {
        int jj = t / nh, ii = 2 * (t % nh) + ((color + B.i0 + B.j0 + jj) & 1);
        if (ii >= ni)
          continue;
        R a[16], cc[16];
        a[0] = 0;
        cc[0] = FIX(ii + 1, jj + 1, 0);
#pragma unroll
        for (int kk = 0; kk < nk; kk++) {
          R At = CIX(cT, ii, jj, kk), Bb = CIX(cB, ii, jj, kk);
          R Cc = CIX(cE, ii, jj, kk) * FIX(ii + 2, jj + 1, kk + 1) +
                 CIX(cN, ii, jj, kk) * FIX(ii + 1, jj + 2, kk + 1) + CIX(cSU, ii, jj, kk) +
                 CIX(cW, ii, jj, kk) * FIX(ii, jj + 1, kk + 1) +
                 CIX(cS, ii, jj, kk) * FIX(ii + 1, jj, kk + 1);
          R term = (R)1 / (CIX(cP, ii, jj, kk) - Bb * a[kk] + (R)1e-30);
          a[kk + 1] = At * term;
          cc[kk + 1] = (Cc + Bb * cc[kk]) * term;
        }
#pragma unroll
        for (int kk = nk - 1; kk >= 0; kk--)
          FIX(ii + 1, jj + 1, kk + 1) = a[kk + 1] * FIX(ii + 1, jj + 1, kk + 2) + cc[kk + 1];
      }
      __syncthreads();
    }
  for (int t = threadIdx.x; t < nk * np; t += blockDim.x) {
    int kk = t / np, r = t % np, jj = r / ni, ii = r % ni;
    f[IDX(B.i0 + ii, B.j0 + jj, B.k0 + kk)] = FIX(ii + 1, jj + 1, kk + 1);
  }
}
__global__ void k_err(Box B, int fv, int slot) {
  CELL(B);
  R v = 0;
  if (G.act[c]) {
    R *f = field(fv, c);
    v = fabs(F(i, j, k) - A(F1, i, j, k));
    if (!isfinite(v))
      v = INFINITY;
  }
  wmax(c, slot, v);
}

__global__ void k_pclear(Box B) {
  DECODE(B);
  A(PP, i, j, k) = 0;
  A(SUK, i, j, k) = 0;
  A(SPK, i, j, k) = 0;
  A(SP, i, j, k) = 0;
}
__global__ void k_pflux(Box B, int dir) {
  DECODE(B);
  int di = dir == 0, dj = dir == 1, dk = dir == 2;
  int i2 = i - di, j2 = j - dj, k2 = k - dk;
#define Q4(v, a1, b1, c1)                                                                          \
  (0.25 * (A(v, i, j, k) + A(v, i2, j2, k2) + A(v, i + (a1), j + (b1), k + (c1)) +                 \
           A(v, i2 + (a1), j2 + (b1), k2 + (c1))))
#define FV(v, a1, b1, c1)                                                                          \
  ((a1) == di && (b1) == dj && (c1) == dk      ? A(v, i, j, k)                                     \
   : (a1) == -di && (b1) == -dj && (c1) == -dk ? A(v, i2, j2, k2)                                  \
                                               : Q4(v, a1, b1, c1))
#define FACES(v, x)                                                                                \
  R x##e = FV(v, 1, 0, 0), x##w = FV(v, -1, 0, 0), x##n = FV(v, 0, 1, 0), x##s = FV(v, 0, -1, 0),  \
    x##t = FV(v, 0, 0, 1), x##b = FV(v, 0, 0, -1)
  FACES(DEN, D);
  FACES(U, U);
  FACES(V, V);
  FACES(W, W);
#define H2(v) (0.5 * (A(v, i, j, k) + A(v, i2, j2, k2)))
  R CXQ = H2(CX), EXQ = H2(EX), SXQ = H2(SX), CYQ = H2(CY), EYQ = H2(EY), SYQ = H2(SY),
    CZQ = H2(CZ), EZQ = H2(EZ), SZQ = H2(SZ);
  R CE = De * (Ue * CXQ + Ve * CYQ + We * CZQ), CW = Dw * (Uw * CXQ + Vw * CYQ + Ww * CZQ);
  R CN = Dn * (Un * EXQ + Vn * EYQ + Wn * EZQ), CS = Ds * (Us * EXQ + Vs * EYQ + Ws * EZQ);
  R CT = Dt * (Ut * SXQ + Vt * SYQ + Wt * SZQ), CB = Db * (Ub * SXQ + Vb * SYQ + Wb * SZQ);
  R val = -(CE - CW + CN - CS + CT - CB) * (A(TJO, i2, j2, k2) + A(TJO, i, j, k)) * 0.5;
  if (dir == 0)
    A(SUK, i, j, k) = val;
  else if (dir == 1)
    A(SPK, i, j, k) = val;
  else
    A(SP, i, j, k) = val;
}
__global__ void k_pcoef(Box B) {
  DECODE(B);
  int II = i - 1, JJ = j - 1, KK = k - 1;
#define DUA(a, b, d) A(DU, a, b, d)
  R GAE = (DUA(i, j, k) + DUA(i, JJ, k) + DUA(i, j, KK) + DUA(i, JJ, KK)) * 0.25;
  R GAW = (DUA(II, j, k) + DUA(II, JJ, k) + DUA(II, j, KK) + DUA(II, JJ, KK)) * 0.25;
  R GAN = (DUA(i, j, k) + DUA(II, j, k) + DUA(i, j, KK) + DUA(II, j, KK)) * 0.25;
  R GAS = (DUA(i, JJ, k) + DUA(II, JJ, k) + DUA(i, JJ, KK) + DUA(II, JJ, KK)) * 0.25;
  R GAT = (DUA(i, j, k) + DUA(II, j, k) + DUA(i, JJ, k) + DUA(II, JJ, k)) * 0.25;
  R GAB = (DUA(i, j, KK) + DUA(II, j, KK) + DUA(i, JJ, KK) + DUA(II, JJ, KK)) * 0.25;
#define T2(v) ((A(v, i, j, k) + A(v, II, JJ, KK)) * 0.5)
  R TQ = T2(TJO);
  A(AE, i, j, k) = TQ * GAE * T2(TXXE);
  A(AW, i, j, k) = TQ * GAW * T2(TXXW);
  A(AN, i, j, k) = TQ * GAN * T2(TYYN);
  A(AS, i, j, k) = TQ * GAS * T2(TYYS);
  A(AT, i, j, k) = TQ * GAT * T2(TZZT);
  A(AB, i, j, k) = TQ * GAB * T2(TZZB);
  R P1 = GAN != 0, P2 = GAS != 0, P3 = GAE != 0, P4 = GAW != 0, P5 = GAT != 0, P6 = GAB != 0;
  R PQ = 0.5 * ((P1 + P2) * (P3 + P4 + P5 + P6) + (P3 + P4) * (P1 + P2 + P5 + P6) +
                (P5 + P6) * (P1 + P2 + P3 + P4));
  if (PQ == 0)
    PQ = 1;
  A(SU, i, j, k) = (P1 * P5 * A(SUK, i, j, k) + P5 * P2 * A(SUK, i, j - 1, k) +
                    P2 * P6 * A(SUK, i, j - 1, k - 1) + P6 * P1 * A(SUK, i, j, k - 1) +
                    P3 * P5 * A(SPK, i, j, k) + P5 * P4 * A(SPK, i - 1, j, k) +
                    P4 * P6 * A(SPK, i - 1, j, k - 1) + P6 * P3 * A(SPK, i, j, k - 1) +
                    P1 * P3 * A(SP, i, j, k) + P3 * P2 * A(SP, i, j - 1, k) +
                    P2 * P4 * A(SP, i - 1, j - 1, k) + P4 * P1 * A(SP, i - 1, j, k)) /
                   PQ;
}
__global__ void k_pdiag(Box B) {
  DECODE(B);
  A(AP, i, j, k) = A(AE, i, j, k) + A(AW, i, j, k) + A(AN, i, j, k) + A(AS, i, j, k) +
                   A(AT, i, j, k) + A(AB, i, j, k);
  if (A(MCF, i, j, k) >= 2) {
    A(AP, i, j, k) = 1;
    A(AE, i, j, k) = A(AW, i, j, k) = A(AN, i, j, k) = A(AS, i, j, k) = A(AT, i, j, k) =
        A(AB, i, j, k) = 0;
    A(SU, i, j, k) = 0;
  }
}
__global__ void k_pblock(Box B, R *wk) {
  extern __shared__ R sh[];
  int c = blockIdx.x;
  if (c >= NB || !G.act[c])
    return;
  int lane = threadIdx.x & 63, part = threadIdx.x >> 6, np_ = blockDim.x >> 6;
  R *a = sh, *b = a + DI, *cc = b + DI, *d = cc + DI, *pa = d + DI, *pb = pa + np_ * DI,
    *pc = pb + np_ * DI;
  for (int i = B.i0 + lane; i <= B.i1; i += 64) {
    R sa = 0, sb = 0, sc = 0;
    for (int j = B.j0 + part; j <= B.j1; j += np_)
      for (int k = B.k0; k <= B.k1; k++) {
        sa -= A(AE, i, j, k);
        sb -= A(AW, i, j, k);
        sc += A(SU, i, j, k);
      }
    pa[part * DI + i] = sa;
    pb[part * DI + i] = sb;
    pc[part * DI + i] = sc;
  }
  __syncthreads();
  for (int i = B.i0 + threadIdx.x; i <= B.i1; i += blockDim.x) {
    R sa = 0, sb = 0, sc = 0;
    for (int q = 0; q < np_; q++) {
      sa += pa[q * DI + i];
      sb += pb[q * DI + i];
      sc += pc[q * DI + i];
    }
    a[i] = sa;
    b[i] = sb;
    cc[i] = sc;
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    int is = B.i0, lt = B.i1;
    for (int i = is; i <= lt; i++)
      d[i] = (-a[i] - b[i]) * 1.5;
    for (int i = is + 1; i <= lt; i++) {
      R t = b[i] / d[i - 1];
      d[i] -= t * a[i - 1];
      cc[i] -= t * cc[i - 1];
    }
    d[lt] = cc[lt] / d[lt];
    for (int ii = is + 1; ii <= lt; ii++) {
      int i = is + lt - ii;
      d[i] = (cc[i] - a[i] * d[i + 1]) / d[i];
    }
  }
  __syncthreads();
  R *w = wk + (size_t)c * DI;
  for (int i = B.i0 + threadIdx.x; i <= B.i1; i += blockDim.x)
    w[i] = d[i];
}
__global__ void k_pblock_set(Box B, R *wk) {
  DECODE(B);
  R *f = field(PP_, c);
  R *d = wk + (size_t)c * DI;
  if ((j == B.j0 && k == B.k0) || (j > B.j0 && k > B.k0))
    F(i, j, k) = d[i];
}
__global__ void k_pupdate(Box B) {
  CELL(B);
  R v = 0;
  if (G.act[c] && A(MCF, i, j, k) < 2) {
    R *f = field(PP_, c);
    R p1 = F(i, j, k) - f[IDX(LT - ((LT + i + j + MT / 2) & 1), MT / 2, 2)];
    A(P, i, j, k) += G.alf[c * 8 + 7] * p1;
    v = fabs(p1);
    if (!isfinite(v))
      v = INFINITY;
  }
  wmax(c, 0, v);
}
__global__ void k_correct(Box B) {
  DECODE(B);
  if (A(MCF, i, j, k) >= 1)
    return;
  R *f = field(PP_, c);
  int II = i + 1, JJ = j + 1, KK = k + 1;
  R PE = (F(II, j, k) + F(II, JJ, k) + F(II, j, KK) + F(II, JJ, KK)) * 0.25,
    PW = (F(i, j, k) + F(i, JJ, k) + F(i, j, KK) + F(i, JJ, KK)) * 0.25;
  R PN = (F(i, JJ, k) + F(II, JJ, k) + F(i, JJ, KK) + F(II, JJ, KK)) * 0.25,
    P5 = (F(i, j, k) + F(II, j, k) + F(i, j, KK) + F(II, j, KK)) * 0.25;
  R PT = (F(i, j, KK) + F(II, j, KK) + F(i, JJ, KK) + F(II, JJ, KK)) * 0.25,
    P3 = (F(i, j, k) + F(II, j, k) + F(i, JJ, k) + F(II, JJ, k)) * 0.25;
  R a = PE - PW, e = PN - P5, s = PT - P3, du = A(DU, i, j, k);
  A(U, i, j, k) -= du * (a * A(CX, i, j, k) + e * A(EX, i, j, k) + s * A(SX, i, j, k));
  A(V, i, j, k) -= du * (a * A(CY, i, j, k) + e * A(EY, i, j, k) + s * A(SY, i, j, k));
  A(W, i, j, k) -= du * (a * A(CZ, i, j, k) + e * A(EZ, i, j, k) + s * A(SZ, i, j, k));
}
__global__ void k_positive(Box B) {
  DECODE(B);
  if (A(DK, i, j, k) <= SMNUM)
    A(DK, i, j, k) = SMNUM;
}
__global__ void k_viscosity(Box B) {
  DECODE(B);
  R *dk = &A(DK, i, j, k), *de = &A(DE, i, j, k);
  if (*dk <= SMNUM)
    *dk = SMNUM;
  if (*de <= SMNUM)
    *de = SMNUM;
  R tv = *de <= SMNUM ? G.visc[c] : A(DEN, i, j, k) * CMU * *dk * *dk / *de + G.visc[c];
  A(VISE, i, j, k) += ALVIS * (tv - A(VISE, i, j, k));
}
__global__ void k_symmetry(Box B) {
  DECODE(B);
#define CP(v)                                                                                      \
  A(v, i, j, 1) = A(v, i, j, 2);                                                                   \
  A(v, i, j, N) = A(v, i, j, N - 1);
  CP(U);
  CP(V);
  CP(TM);
  CP(DK);
  CP(DE);
  A(W, i, j, 1) = 0;
  A(W, i, j, N) = 0;
}
__global__ void k_outflow(R *flow) {
  int c = blockIdx.x;
  if (c >= NB || !G.act[c])
    return;
  int i = LT, nj = MT, nk = NTT;
  double fl = 0, ar = 0;
  for (int t = threadIdx.x; t < nj * nk; t += blockDim.x) {
    int j = 2 + t % nj, k = 2 + t / nj;
#define C4(v) (A(v, i, j, k) + A(v, i, j - 1, k) + A(v, i, j, k - 1) + A(v, i, j - 1, k - 1))
    R uc = C4(U) * 0.25f, dc = C4(DEN) * 0.25f;
    R p1 = (A(X, i, j, k) + A(X, i, j, k - 1) - A(X, i, j - 1, k) - A(X, i, j - 1, k - 1)) * 0.5f,
      p2 = (A(Y, i, j, k) + A(Y, i, j, k - 1) - A(Y, i, j - 1, k) - A(Y, i, j - 1, k - 1)) * 0.5f,
      p3 = (A(Z, i, j, k) + A(Z, i, j, k - 1) - A(Z, i, j - 1, k) - A(Z, i, j - 1, k - 1)) * 0.5f;
    R q1 = (A(X, i, j, k) + A(X, i, j - 1, k) - A(X, i, j, k - 1) - A(X, i, j - 1, k - 1)) * 0.5f,
      q2 = (A(Y, i, j, k) + A(Y, i, j - 1, k) - A(Y, i, j, k - 1) - A(Y, i, j - 1, k - 1)) * 0.5f,
      q3 = (A(Z, i, j, k) + A(Z, i, j - 1, k) - A(Z, i, j, k - 1) - A(Z, i, j - 1, k - 1)) * 0.5f;
    R area = sqrt(p1 * p1 + p2 * p2 + p3 * p3) * sqrt(q1 * q1 + q2 * q2 + q3 * q3);
    fl += dc * area * uc;
    ar += dc * area;
  }
  __shared__ double s1[256], s2[256];
  s1[threadIdx.x] = fl;
  s2[threadIdx.x] = ar;
  __syncthreads();
  for (int o = blockDim.x / 2; o > 0; o >>= 1) {
    if (threadIdx.x < o) {
      s1[threadIdx.x] += s1[threadIdx.x + o];
      s2[threadIdx.x] += s2[threadIdx.x + o];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0)
    flow[c] = (s1[0] - G.flowin[c]) / s2[0];
}
__global__ void k_outlet(Box B, R *uinc) {
  DECODE(B);
  if (G.flg[c] & 16)
    return;
  int I = LT;
  A(U, I + 1, j, k) = A(U, I, j, k) - uinc[c];
  A(V, I + 1, j, k) = A(V, I, j, k);
  A(W, I + 1, j, k) = A(W, I, j, k);
  A(DK, I + 1, j, k) = A(DK, I, j, k);
  A(DE, I + 1, j, k) = A(DE, I, j, k);
}

#define MAXDK 16

static const char *argv0 = "cns2d";

static void warn(const char *fmt, ...) {
  va_list a;
  fprintf(stderr, "%s: ", argv0);
  va_start(a, fmt);
  vfprintf(stderr, fmt, a);
  va_end(a);
  fputc('\n', stderr);
}

static void die(const char *fmt, ...) {
  va_list a;
  fprintf(stderr, "%s: ", argv0);
  va_start(a, fmt);
  vfprintf(stderr, fmt, a);
  va_end(a);
  fputc('\n', stderr);
  exit(1);
}

static void *emalloc(size_t n) {
  void *p = calloc(n ? n : 1, 1);
  if (!p)
    die("out of memory: %zu bytes", n);
  return p;
}

static void *erealloc(void *p, size_t n) {
  p = realloc(p, n);
  if (!p)
    die("out of memory: %zu bytes", n);
  return p;
}

static FILE *eopen(const char *path, const char *mode) {
  FILE *f = fopen(path, mode);
  if (!f)
    die("%s: %s", path, strerror(errno));
  return f;
}

static void egets(char *s, int n, FILE *f, const char *path) {
  if (!fgets(s, n, f))
    die("%s: %s", path, ferror(f) ? strerror(errno) : "unexpected end of file");
}

static void cucheck(cudaError_t e, const char *call, int line) {
  if (e != cudaSuccess)
    die("line %d: %s: %s", line, call, cudaGetErrorString(e));
}
#define CK(call) cucheck((call), #call, __LINE__)

typedef struct {
  int i, j, k, t;
  float yn, yn1, sx, sy, sz;
} Wall;

typedef struct {
  char dir[1024];
  R hyb, c1;
  double relax, relaxp;
  int flags, block[6];
  int l, m, n, ig, nlimt, nswp;
  float visc, flowin, alr[8];
  Wall *walls;
  int nwall;
} Case;

static int hDI, hDJ, hDK, hNT, hNB, hMAXW;
static int hL, hM, hN, hLT, hMT, hNTT, hIG;
static Case *cases;
static S hs;
static R *d_wk, *d_flow;
#define PBLOCK_SHARED (sizeof(R) * 28 * hDI)
static cudaStream_t stream;
static int nsweep[8];
static R sigma[8];

static int nblk(long n) {
  return (int)((n + 255) / 256);
}
static long cnt(Box b) {
  long n = (long)(b.i1 - b.i0 + 1) * (b.j1 - b.j0 + 1) * (b.k1 - b.k0 + 1);
  return n > 0 ? n * hNB : 0;
}
#define RUN(kern, b, ...)                                                                          \
  do {                                                                                             \
    long n_ = cnt(b);                                                                              \
    if (n_)                                                                                        \
      kern<<<nblk(n_), 256, 0, stream>>>(b, ##__VA_ARGS__);                                        \
  } while (0)

static void lines(Box b, int fv, int nsw) {
  static int shared_set;
  int ni = b.i1 - b.i0 + 1, nj = b.j1 - b.j0 + 1, nk = b.k1 - b.k0 + 1;
  size_t sz = sizeof(R) * ((size_t)8 * nk * ni * nj + (size_t)(nk + 2) * (ni + 2) * (nj + 2));
  if (sz <= 227 * 1024) {
    if (!shared_set) {
      cudaFuncSetAttribute(k_lines, cudaFuncAttributeMaxDynamicSharedMemorySize, 227 * 1024);
      shared_set = 1;
    }
    k_lines<<<hNB, 1024, sz, stream>>>(b, fv, nsw);
  } else
    for (int s = 0; s < nsw; s++)
      for (int color = 0; color < 2; color++)
        k_line<<<nblk((long)ni * nj * hNB), 256, 0, stream>>>(b, fv, color);
}

static void solveq(int ie, int fv) {
  Box in = {2, hLT, 2, hMT, 2, hNTT};
  Box all = {1, hL, 1, hM, 1, hN};
  Box outlet = {hLT, hLT, 2, hMT, 2, hNTT};
  if (ie == 1 || ie == 6) {
    RUN(k_diffusivity, all, sigma[ie]);
    RUN(k_coef, in, fv);
  }
  if (ie == 1)
    RUN(k_src_u, in);
  else if (ie == 2)
    RUN(k_src_vw, in, 0);
  else if (ie == 3)
    RUN(k_src_vw, in, 1);
  else if (ie == 5)
    RUN(k_src_k, in);
  else if (ie == 6)
    RUN(k_src_e, in);
  if (hIG == 2)
    k_wall<<<nblk((long)hNB * hMAXW), 256, 0, stream>>>(ie, fv);
  RUN(k_outlet_links, outlet);
  RUN(k_assemble, in, fv, ie);
  lines(in, fv, nsweep[ie]);
  RUN(k_err, in, fv, ie);
}

static void solvep(void) {
  Box in = {2, hL, 2, hM, 2, hN};
  Box cells = {2, hLT, 2, hMT, 2, hNTT};
  RUN(k_pclear, in);
  RUN(k_pflux, in, 0);
  RUN(k_pflux, in, 1);
  RUN(k_pflux, in, 2);
  RUN(k_pcoef, in);
  RUN(k_pdiag, in);
  k_pblock<<<hNB, 512, PBLOCK_SHARED, stream>>>(in, d_wk);
  RUN(k_pblock_set, in, d_wk);
  lines(in, PP_, nsweep[0]);
  RUN(k_pupdate, in);
  RUN(k_correct, cells);
}

static void iterate(void) {
  Box sym = {1, hL, 2, hMT, 1, 1};
  Box plane = {1, 1, 2, hMT, 2, hNTT};
  Box all = {1, hL, 1, hM, 1, hN};
  RUN(k_symmetry, sym);
  k_outflow<<<hNB, 256, 0, stream>>>(d_flow);
  RUN(k_outlet, plane, d_flow);
  solveq(1, U_);
  solveq(2, V_);
  solveq(3, W_);
  if (hIG == 2) {
    solveq(5, DK_);
    RUN(k_positive, all);
    solveq(6, DE_);
  }
  solvep();
  if (hIG == 2)
    RUN(k_viscosity, all);
}

#define ID(i, j, k) ((i) - 1 + hDI * ((j) - 1 + hDJ * ((k) - 1)))

static float fx(const char *s, int a, int w) {
  char b[32];
  int n = (int)strlen(s);
  while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r'))
    n--;
  if (a >= n)
    return 0;
  if (a + w > n)
    w = n - a;
  memcpy(b, s + a, w);
  b[w] = 0;
  return (float)atof(b);
}
static int ix(const char *s, int a, int w) {
  return (int)fx(s, a, w);
}

static const int stored[] = {X_,    Y_,    Z_,    TJO_,  CX_,   CY_,   CZ_,   EX_,   EY_,
                             EZ_,   SX_,   SY_,   SZ_,   TXXE_, TXXW_, TYYN_, TYYS_, TZZT_,
                             TZZB_, TYXE_, TYXW_, TYZT_, TYZB_, TXYN_, TXYS_, TXZT_, TXZB_,
                             TZXE_, TZXW_, TZYN_, TZYS_, MCF_,  U_,    V_,    W_,    P_,
                             DK_,   DE_,   VISE_, DEN_,  TM_};
#define NSTORED ((int)(sizeof stored / sizeof stored[0]))

static R *host[NA_];
static int *ibc, *jbc, *kbc, *ity, *nw;
static R *hyn, *hyn1, *sxw, *syw, *szw;
static Wall *wl;
static int wcap;

static float *g[NA_], *vm[10];
static float *X, *Y, *Z;
static int *MC, L1, L2, M1, M2, N1, N2, nwall;

static void walval(int is, int it, int js, int jt, int ks, int kt, float *F) {
  int i, j, k, w;
  for (j = js; j <= jt; j++)
    for (k = ks; k <= kt; k++) {
      F[ID(is - 1, j, k)] = F[ID(is, j, k)];
      F[ID(it + 1, j, k)] = F[ID(it, j, k)];
    }
  for (i = is - 1; i <= it + 1; i++)
    for (k = ks; k <= kt; k++) {
      F[ID(i, js - 1, k)] = F[ID(i, js, k)];
      F[ID(i, jt + 1, k)] = F[ID(i, jt, k)];
    }
  for (i = is - 1; i <= it + 1; i++)
    for (j = js - 1; j <= jt + 1; j++) {
      F[ID(i, j, ks - 1)] = F[ID(i, j, ks)];
      F[ID(i, j, kt + 1)] = F[ID(i, j, kt)];
    }
  for (w = 0; w < nwall; w++) {
    int i = wl[w].i, j = wl[w].j, k = wl[w].k, o = ID(i, j, k);
    switch (wl[w].t) {
    case 1:
      F[ID(i, j + 1, k)] = F[o];
      break;
    case 2:
      F[ID(i, j - 1, k)] = F[o];
      break;
    case 3:
      F[ID(i + 1, j, k)] = F[o];
      break;
    case 4:
      F[ID(i - 1, j, k)] = F[o];
      break;
    case 5:
      F[ID(i, j, k + 1)] = F[o];
      break;
    case 6:
      F[ID(i, j, k - 1)] = F[o];
      break;
    }
  }
}

static void wall(int i, int j, int k, int t, float p1, float p2, float p3, int a1, int b1, int c1,
                 int a2, int b2, int c2) {
  float pq = sqrtf(p1 * p1 + p2 * p2 + p3 * p3), q1, q2, q3, aa, bb, cc = 1.0f, coth, yn;
  int o = ID(i, j, k), o1 = ID(a1, b1, c1), o2 = ID(a2, b2, c2), n;
  p1 = p1 / pq;
  p2 = p2 / pq;
  p3 = p3 / pq;
  if (nwall == wcap) {
    wcap = wcap ? 2 * wcap : 256;
    wl = (Wall *)erealloc(wl, wcap * sizeof(Wall));
  }
  n = nwall++;
  wl[n].i = i;
  wl[n].j = j;
  wl[n].k = k;
  wl[n].t = t;
  wl[n].sx = sqrtf(1.f - p1 * p1);
  wl[n].sy = sqrtf(1.f - p2 * p2);
  wl[n].sz = sqrtf(1.f - p3 * p3);
  q1 = X[o] - X[o1];
  q2 = Y[o] - Y[o1];
  q3 = Z[o] - Z[o1];
  aa = sqrtf((q1 - p1) * (q1 - p1) + (q2 - p2) * (q2 - p2) + (q3 - p3) * (q3 - p3));
  bb = sqrtf(q1 * q1 + q2 * q2 + q3 * q3);
  coth = (bb * bb + cc * cc - aa * aa) / (2 * bb * cc);
  yn = bb * fabsf(coth);
  wl[n].yn = yn;
  q1 = X[o] - X[o2];
  q2 = Y[o] - Y[o2];
  q3 = Z[o] - Z[o2];
  bb = sqrtf(q1 * q1 + q2 * q2 + q3 * q3);
  aa = sqrtf((q1 - p1) * (q1 - p1) + (q2 - p2) * (q2 - p2) + (q3 - p3) * (q3 - p3));
  coth = (bb * bb + cc * cc - aa * aa) / (2 * bb * cc);
  wl[n].yn1 = (bb * fabsf(coth) + yn) * 0.5f;
}

#define WCROSS(A1, A2, B1, B2)                                                                     \
  p1 = (Y[A1] - Y[A2]) * (Z[B1] - Z[B2]) - (Z[A1] - Z[A2]) * (Y[B1] - Y[B2]);                      \
  p2 = (Z[A1] - Z[A2]) * (X[B1] - X[B2]) - (X[A1] - X[A2]) * (Z[B1] - Z[B2]);                      \
  p3 = (X[A1] - X[A2]) * (Y[B1] - Y[B2]) - (Y[A1] - Y[A2]) * (X[B1] - X[B2])

static void wallj(int i, int j, int k, int t, int j1, int j2) {
  int i1 = i + 1, i2 = i - 1, k1 = k + 1, k2 = k - 1;
  float p1, p2, p3;
  if (t == 1) {
    if (i == L1)
      i2 = i;
    if (i == L2)
      i1 = i;
    if (k == N1)
      k2 = k;
    if (k == N2)
      k1 = k;
  }
  WCROSS(ID(i1, j, k1), ID(i1, j, k2), ID(i1, j, k2), ID(i2, j, k));
  wall(i, j, k, t, p1, p2, p3, i, j1, k, i, j2, k);
}
static void walli(int i, int j, int k, int t, int i1, int i2) {
  int j1 = j + 1, j2 = j - 1, k1 = k + 1, k2 = k - 1;
  float p1, p2, p3;
  if (j == M1)
    j2 = j;
  if (j == M2)
    j1 = j;
  if (k == N1)
    k2 = k;
  if (k == N2)
    k1 = k;
  WCROSS(ID(i, j1, k1), ID(i, j1, k2), ID(i, j1, k2), ID(i, j2, k));
  wall(i, j, k, t, p1, p2, p3, i1, j, k, i2, j, k);
}
static void wallk(int i, int j, int k, int t, int k1, int k2) {
  int i1 = i + 1, i2 = i - 1, j1 = j + 1, j2 = j - 1;
  float p1, p2, p3;
  if (i == L1)
    i2 = i;
  if (i == L2)
    i1 = i;
  if (j == M1)
    j2 = j;
  if (j == M2)
    j1 = j;
  WCROSS(ID(i1, j1, k), ID(i1, j2, k), ID(i1, j2, k), ID(i2, j, k));
  wall(i, j, k, t, p1, p2, p3, i, j, k1, i, j, k2);
}

static float edge(float *A, int i, int j, int k, int a) {
  int o = ID(i, j, k), oj = ID(i, j - 1, k), ok = ID(i, j, k - 1), ojk = ID(i, j - 1, k - 1),
      oi = ID(i - 1, j, k), oij = ID(i - 1, j - 1, k), oik = ID(i - 1, j, k - 1),
      oijk = ID(i - 1, j - 1, k - 1);
  if (a == 0)
    return (A[o] + A[oj] + A[ok] + A[ojk] - A[oi] - A[oij] - A[oik] - A[oijk]) * 0.25f;
  if (a == 1)
    return (A[o] + A[oi] + A[ok] + A[oik] - A[oj] - A[oij] - A[ojk] - A[oijk]) * 0.25f;
  return (A[o] + A[oi] + A[oj] + A[oij] - A[ok] - A[oik] - A[ojk] - A[oijk]) * 0.25f;
}

#define WH2(A, nb) (((A)[nb] + (A)[o]) * 0.5f)

static void read_control(Case *k) {
  char s[512], p[1100];
  FILE *f;
  int q;
  snprintf(p, sizeof p, "%s/control", k->dir);
  f = eopen(p, "r");
  egets(s, sizeof s, f, p);
  k->nlimt = ix(s, 0, 5);
  k->ig = ix(s, 5, 5);
  k->nswp = ix(s, 10, 5);
  egets(s, sizeof s, f, p);
  for (q = 0; q < 8; q++)
    k->alr[q] = fx(s, 7 * q, 7);
  egets(s, sizeof s, f, p);
  k->visc = 1.f / fx(s, 0, 7);
  fclose(f);
  if (k->ig != 1 && k->ig != 2)
    die("%s: model %d is neither 1 (laminar) nor 2 (turbulent)", p, k->ig);
  if (k->nlimt < 1)
    die("%s: iteration limit %d", p, k->nlimt);
  if (k->nswp < 1)
    die("%s: %d pressure sweeps", p, k->nswp);
  if (!isfinite(k->visc) || k->visc <= 0)
    die("%s: bad Reynolds number", p);
}

static void setup(int c) {
  Case *K = &cases[c];
  char s[512], p[1100];
  FILE *ff;
  int i, j, k, q, n, NN = hNT, L, M, N, LO, MO, NO, LT, MT, NT, wdist;
  float visc, alvis, cmu = 0.09f, smnum = 1.E-30f, flowin = 0.0f;
  float *U = g[U_], *V = g[V_], *W = g[W_], *DK = g[DK_], *DE = g[DE_], *VISE = g[VISE_];
  float *DEN = g[DEN_], *TJO = g[TJO_], *CX = g[CX_], *CY = g[CY_], *CZ = g[CZ_], *EX = g[EX_];
  float *EY = g[EY_], *EZ = g[EZ_], *SX = g[SX_], *SY = g[SY_], *SZ = g[SZ_];
  float *cell[10] = {CX, CY, CZ, EX, EY, EZ, SX, SY, SZ, TJO};
  read_control(K);
  visc = K->visc;
  alvis = K->alr[6];
  snprintf(p, sizeof p, "%s/field", K->dir);
  ff = eopen(p, "r");
  egets(s, sizeof s, ff, p);
  L = ix(s, 0, 5);
  M = ix(s, 5, 5);
  N = ix(s, 10, 5);
  L1 = ix(s, 15, 5);
  L2 = ix(s, 20, 5);
  M1 = ix(s, 25, 5);
  M2 = ix(s, 30, 5);
  N1 = ix(s, 35, 5);
  N2 = ix(s, 40, 5);
  K->l = L;
  K->block[0] = L1;
  K->block[1] = L2;
  K->block[2] = M1;
  K->block[3] = M2;
  K->block[4] = N1;
  K->block[5] = N2;
  K->m = M;
  K->n = N;
  LO = L + 1;
  MO = M + 1;
  NO = N + 1;
  LT = L - 1;
  MT = M - 1;
  NT = N - 1;
  if (L < 3 || M < 3 || N < 3)
    die("%s: grid %dx%dx%d is too small", p, L, M, N);
  if (LO > hDI || MO > hDJ || NO > hDK)
    die("%s: grid %dx%dx%d does not fit %dx%dx%d", p, L, M, N, hDI, hDJ, hDK);
  for (q = 0; q < NSTORED; q++)
    memset(g[stored[q]], 0, sizeof(float) * NN);
  for (q = 0; q < 10; q++)
    memset(vm[q], 0, sizeof(float) * NN);
  memset(MC, 0, sizeof(int) * NN);
  wl = NULL;
  wcap = 0;
  X = g[X_];
  Y = g[Y_];
  Z = g[Z_];
  for (i = 1; i <= LO; i++)
    for (j = 1; j <= MO; j++)
      for (k = 1; k <= NO; k++) {
        DEN[ID(i, j, k)] = 1.0f;
        VISE[ID(i, j, k)] = visc;
      }
  for (k = 1; k <= N; k++)
    for (i = 1; i <= L; i++)
      for (j = 1; j <= M; j++) {
        int o = ID(i, j, k);
        do {
          egets(s, sizeof s, ff, p);
        } while (strlen(s) < 100);
        X[o] = fx(s, 0, 8);
        Y[o] = fx(s, 8, 8);
        Z[o] = fx(s, 16, 8);
        U[o] = fx(s, 24, 12);
        V[o] = fx(s, 36, 12);
        W[o] = fx(s, 48, 12);
        g[P_][o] = fx(s, 62, 11);
        g[TM_][o] = fx(s, 73, 11);
        DK[o] = fx(s, 84, 11);
        DE[o] = fx(s, 95, 11);
      }
  fclose(ff);

  for (i = 1; i <= LO; i++)
    for (j = 1; j <= MO; j++)
      for (k = 1; k <= NO; k++) {
        int m = (j == 1 || (j == M && !(K->flags & 16))) ? 1 : 0;
        if (i >= L1 && i <= L2 && j >= M1 && j <= M2 && k >= N1 && k <= N2)
          m = 1;
        if (i > L1 && i <= L2 && j > M1 && j <= M2 && k > N1 && k <= N2)
          m = 2;
        MC[ID(i, j, k)] = m;
      }

  nwall = 0;
  wdist = (K->flags & 4) != 0;
  for (i = 2; i <= LT; i++)
    for (j = 2; j <= MT; j++)
      for (k = 2; k <= NT; k++) {
        if (MC[ID(i, j, k)] != 0)
          continue;
        if (MC[ID(i + 1, j, k)] + MC[ID(i - 1, j, k)] + MC[ID(i, j + 1, k)] + MC[ID(i, j - 1, k)] +
                MC[ID(i, j, k + 1)] + MC[ID(i, j, k - 1)] ==
            0)
          continue;
        if (MC[ID(i, j + 1, k)])
          wallj(i, j, k, 1, wdist ? j + 1 : j - 1, wdist ? j - 1 : j - 2);
        if (MC[ID(i, j - 1, k)])
          wallj(i, j, k, 2, wdist ? j - 1 : j + 1, wdist ? j + 1 : j + 2);
        if (MC[ID(i + 1, j, k)])
          walli(i, j, k, 3, wdist ? i + 1 : i - 1, wdist ? i - 1 : i - 2);
        if (MC[ID(i - 1, j, k)] && !((K->flags & 8) && j == M2))
          walli(i, j, k, 4, wdist ? i - 1 : i + 1, wdist ? i + 1 : i + 2);
        if (MC[ID(i, j, k + 1)])
          wallk(i, j, k, 5, wdist ? k + 1 : k - 1, wdist ? k - 1 : k - 2);
        if (MC[ID(i, j, k - 1)])
          wallk(i, j, k, 6, wdist ? k - 1 : k + 1, wdist ? k + 1 : k + 2);
      }

  for (i = 1; i <= L; i++)
    for (j = 1; j <= M; j++)
      for (k = 1; k <= N; k++) {
        int o = ID(i, j, k);
        if (MC[o]) {
          DK[o] = 0;
          DE[o] = 0;
          U[o] = 0;
          V[o] = 0;
          W[o] = 0;
        }
      }

  for (i = 2; i <= L; i++)
    for (j = 2; j <= M; j++)
      for (k = 2; k <= N; k++) {
        float p1 = edge(X, i, j, k, 0), p2 = edge(X, i, j, k, 1), p3 = edge(X, i, j, k, 2);
        float q1 = edge(Y, i, j, k, 0), q2 = edge(Y, i, j, k, 1), q3 = edge(Y, i, j, k, 2);
        float r1 = edge(Z, i, j, k, 0), r2 = edge(Z, i, j, k, 1), r3 = edge(Z, i, j, k, 2);
        float ptr =
            1.0f / (p1 * (q2 * r3 - q3 * r2) - p2 * (q1 * r3 - q3 * r1) + p3 * (q1 * r2 - q2 * r1));
        int o = ID(i, j, k);
        vm[0][o] = ptr * (q2 * r3 - q3 * r2);
        vm[1][o] = -ptr * (p2 * r3 - p3 * r2);
        vm[2][o] = ptr * (p2 * q3 - p3 * q2);
        vm[3][o] = -ptr * (q1 * r3 - q3 * r1);
        vm[4][o] = ptr * (p1 * r3 - p3 * r1);
        vm[5][o] = -ptr * (p1 * q3 - p3 * q1);
        vm[6][o] = ptr * (q1 * r2 - q2 * r1);
        vm[7][o] = -ptr * (p1 * r2 - p2 * r1);
        vm[8][o] = ptr * (p1 * q2 - p2 * q1);
        vm[9][o] = fabsf(1.0f / ptr);
      }
  for (q = 0; q < 10; q++)
    walval(2, L, 2, M, 2, N, vm[q]);
  for (i = 1; i <= L; i++)
    for (j = 1; j <= M; j++)
      for (k = 1; k <= N; k++)
        for (q = 0; q < 10; q++) {
          float *A = vm[q];
          cell[q][ID(i, j, k)] =
              (A[ID(i, j, k)] + A[ID(i, j + 1, k)] + A[ID(i, j, k + 1)] + A[ID(i, j + 1, k + 1)] +
               A[ID(i + 1, j, k)] + A[ID(i + 1, j + 1, k)] + A[ID(i + 1, j, k + 1)] +
               A[ID(i + 1, j + 1, k + 1)]) *
              0.125f;
        }

  for (i = 2; i <= LT; i++)
    for (j = 2; j <= MT; j++)
      for (k = 2; k <= NT; k++) {
        int o = ID(i, j, k), e_ = ID(i + 1, j, k), w_ = ID(i - 1, j, k), n_ = ID(i, j + 1, k),
            s_ = ID(i, j - 1, k), t_ = ID(i, j, k + 1), b_ = ID(i, j, k - 1);
        float cx = CX[o], cy = CY[o], cz = CZ[o], ex = EX[o], ey = EY[o], ez = EZ[o], sx = SX[o],
              sy = SY[o], sz = SZ[o];
        g[TXXE_][o] = WH2(CX, e_) * cx + WH2(CY, e_) * cy + WH2(CZ, e_) * cz;
        g[TXXW_][o] = WH2(CX, w_) * cx + WH2(CY, w_) * cy + WH2(CZ, w_) * cz;
        g[TYYN_][o] = WH2(EX, n_) * ex + WH2(EY, n_) * ey + WH2(EZ, n_) * ez;
        g[TYYS_][o] = WH2(EX, s_) * ex + WH2(EY, s_) * ey + WH2(EZ, s_) * ez;
        g[TZZT_][o] = WH2(SX, t_) * sx + WH2(SY, t_) * sy + WH2(SZ, t_) * sz;
        g[TZZB_][o] = WH2(SX, b_) * sx + WH2(SY, b_) * sy + WH2(SZ, b_) * sz;
        g[TYXE_][o] = (WH2(EX, e_) * cx + WH2(EY, e_) * cy + WH2(EZ, e_) * cz) * 0.25f;
        g[TYXW_][o] = (WH2(EX, w_) * cx + WH2(EY, w_) * cy + WH2(EZ, w_) * cz) * 0.25f;
        g[TYZT_][o] = (WH2(EX, t_) * sx + WH2(EY, t_) * sy + WH2(EZ, t_) * sz) * 0.25f;
        g[TYZB_][o] = (WH2(EX, b_) * sx + WH2(EY, b_) * sy + WH2(EZ, b_) * sz) * 0.25f;
        g[TXYN_][o] = (WH2(CX, n_) * ex + WH2(CY, n_) * ey + WH2(CZ, n_) * ez) * 0.25f;
        g[TXYS_][o] = (WH2(CX, s_) * ex + WH2(CY, s_) * ey + WH2(CZ, s_) * ez) * 0.25f;
        g[TXZT_][o] = (WH2(CX, t_) * sx + WH2(CY, t_) * sy + WH2(CZ, t_) * sz) * 0.25f;
        g[TXZB_][o] = (WH2(CX, b_) * sx + WH2(CY, b_) * sy + WH2(CZ, b_) * sz) * 0.25f;
        g[TZXE_][o] = (WH2(SX, e_) * cx + WH2(SY, e_) * cy + WH2(SZ, e_) * cz) * 0.25f;
        g[TZXW_][o] = (WH2(SX, w_) * cx + WH2(SY, w_) * cy + WH2(SZ, w_) * cz) * 0.25f;
        g[TZYN_][o] = (WH2(SX, n_) * ex + WH2(SY, n_) * ey + WH2(SZ, n_) * ez) * 0.25f;
        g[TZYS_][o] = (WH2(SX, s_) * ex + WH2(SY, s_) * ey + WH2(SZ, s_) * ez) * 0.25f;
      }
  for (q = TXXE_; q <= TZYS_; q++)
    walval(2, LT, 2, MT, 2, NT, g[q]);

  if (K->ig == 2)
    for (i = 1; i <= L; i++)
      for (j = 1; j <= M; j++)
        for (k = 1; k <= N; k++) {
          int o = ID(i, j, k);
          float tv;
          if (DK[o] <= smnum)
            DK[o] = smnum;
          if (DE[o] <= smnum)
            DE[o] = smnum;
          tv = DE[o] <= smnum ? visc : DEN[o] * cmu * DK[o] * DK[o] / DE[o] + visc;
          VISE[o] = VISE[o] + alvis * (tv - VISE[o]);
        }

  for (j = 2; j <= M; j++)
    for (k = 2; k <= N; k++) {
      int a = ID(1, j, k), b = ID(1, j - 1, k), cc = ID(1, j, k - 1), dd = ID(1, j - 1, k - 1);
      float uc = (U[a] + U[b] + U[cc] + U[dd]) * 0.25f,
            denc = (DEN[a] + DEN[b] + DEN[cc] + DEN[dd]) * 0.25f;
      float p1 = (X[a] + X[cc] - X[b] - X[dd]) * 0.5f, p2 = (Y[a] + Y[cc] - Y[b] - Y[dd]) * 0.5f,
            p3 = (Z[a] + Z[cc] - Z[b] - Z[dd]) * 0.5f;
      float q1 = (X[a] + X[b] - X[cc] - X[dd]) * 0.5f, q2 = (Y[a] + Y[b] - Y[cc] - Y[dd]) * 0.5f,
            q3 = (Z[a] + Z[b] - Z[cc] - Z[dd]) * 0.5f;
      flowin = flowin +
               denc * sqrtf(p1 * p1 + p2 * p2 + p3 * p3) * sqrtf(q1 * q1 + q2 * q2 + q3 * q3) * uc;
    }
  K->flowin = flowin;

  {
    size_t off = (size_t)c * NN;
    for (q = 0; q < NSTORED; q++)
      if (stored[q] != MCF_)
        for (n = 0; n < NN; n++)
          host[stored[q]][off + n] = g[stored[q]][n];
    for (n = 0; n < NN; n++)
      host[MCF_][off + n] = MC[n];
    K->walls = wl;
    K->nwall = nwall;
  }
}

static void read_list(const char *path) {
  char ln[2048], extra[2];
  int cap = 64, line = 0, n;
  FILE *f = eopen(path, "r");
  cases = (Case *)emalloc(cap * sizeof(Case));
  hNB = 0;
  while (fgets(ln, sizeof ln, f)) {
    Case *k = &cases[hNB];
    double hyb, c1 = 1.43;
    line++;
    if (strspn(ln, " \t\r\n") == strlen(ln))
      continue;
    memset(k, 0, sizeof *k);
    k->relax = k->relaxp = -1;
    n = sscanf(ln, "%1023s %lf %lf %lf %d %lf %1s", k->dir, &hyb, &k->relax, &k->relaxp, &k->flags,
               &c1, extra);
    if (n < 2 || n > 6)
      die("%s:%d: expected: dir hyb [relax relaxp flags c1]", path, line);
    if (hyb <= 0 || hyb > 1)
      die("%s:%d: hybrid factor %g is not in (0, 1]", path, line, hyb);
    k->hyb = hyb;
    k->c1 = c1;
    if (k->flags < 0 || k->flags > 31)
      die("%s:%d: flags %d is not in [0, 31]", path, line, k->flags);
    if (++hNB == cap) {
      cap *= 2;
      cases = (Case *)erealloc(cases, cap * sizeof(Case));
    }
  }
  if (ferror(f))
    die("%s: %s", path, strerror(errno));
  fclose(f);
  if (hNB == 0)
    die("%s: no cases", path);
}

static R *ral(size_t n) {
  return (R *)emalloc(n * sizeof(R));
}
static int *ial(size_t n) {
  return (int *)emalloc(n * sizeof(int));
}
static float *fal(size_t n) {
  return (float *)emalloc(n * sizeof(float));
}
static void *dalloc(size_t bytes) {
  void *d;
  CK(cudaMalloc(&d, bytes));
  return d;
}
static void *upload(const void *h, size_t bytes) {
  void *d = dalloc(bytes);
  CK(cudaMemcpy(d, h, bytes, cudaMemcpyHostToDevice));
  return d;
}
#define SETI(sym, v)                                                                               \
  do {                                                                                             \
    int t_ = (v);                                                                                  \
    CK(cudaMemcpyToSymbol(sym, &t_, sizeof t_));                                                   \
  } while (0)
#define SETR(sym, v)                                                                               \
  do {                                                                                             \
    R t_ = (v);                                                                                    \
    CK(cudaMemcpyToSymbol(sym, &t_, sizeof t_));                                                   \
  } while (0)

static void setup_all(void) {
  size_t nw_;
  int c, q, n;
  for (q = 0; q < NSTORED; q++) {
    host[stored[q]] = ral((size_t)hNB * hNT);
    g[stored[q]] = fal(hNT);
  }
  for (q = 0; q < 10; q++)
    vm[q] = fal(hNT);
  MC = ial(hNT);
  for (c = 0; c < hNB; c++) {
    Case *k = &cases[c], *k0 = &cases[0];
    setup(c);
    if (k->l != k0->l || k->m != k0->m || k->n != k0->n)
      die("%s: grid %dx%dx%d, %s has %dx%dx%d", k->dir, k->l, k->m, k->n, k0->dir, k0->l, k0->m,
          k0->n);
    if (k->ig != k0->ig || k->nlimt != k0->nlimt || k->nswp != k0->nswp)
      die("%s: model, iteration limit or sweeps differ from %s", k->dir, k0->dir);
  }
  hMAXW = 1;
  for (c = 0; c < hNB; c++)
    if (cases[c].nwall > hMAXW)
      hMAXW = cases[c].nwall;
  nw_ = (size_t)hNB * hMAXW;
  ibc = ial(nw_);
  jbc = ial(nw_);
  kbc = ial(nw_);
  ity = ial(nw_);
  nw = ial(hNB);
  hyn = ral(nw_);
  hyn1 = ral(nw_);
  sxw = ral(nw_);
  syw = ral(nw_);
  szw = ral(nw_);
  for (c = 0; c < hNB; c++) {
    Case *k = &cases[c];
    for (n = 0; n < k->nwall; n++) {
      size_t o = (size_t)c * hMAXW + n;
      Wall *w = &k->walls[n];
      ibc[o] = w->i;
      jbc[o] = w->j;
      kbc[o] = w->k;
      ity[o] = w->t;
      hyn[o] = w->yn;
      hyn1[o] = w->yn1;
      sxw[o] = w->sx;
      syw[o] = w->sy;
      szw[o] = w->sz;
    }
    nw[c] = k->nwall;
    free(k->walls);
    k->walls = NULL;
  }
  hL = cases[0].l;
  hM = cases[0].m;
  hN = cases[0].n;
  hLT = hL - 1;
  hMT = hM - 1;
  hNTT = hN - 1;
  hIG = cases[0].ig;
  nsweep[0] = cases[0].nswp;
  nsweep[1] = nsweep[2] = nsweep[3] = 7;
  nsweep[5] = nsweep[6] = 5;
  sigma[1] = 1.0f;
  sigma[5] = 1.0f;
  sigma[6] = 1.3f;
}

static void upload_all(void) {
  size_t nw_ = (size_t)hNB * hMAXW, nf = sizeof(R) * hNB * hNT;
  R *cfac = ral(hNB), *c1 = ral(hNB), *visc = ral(hNB), *flowin = ral(hNB);
  R *alf = ral((size_t)hNB * 8), *zero = ral(nw_);
  int *flg = ial(hNB), *act = ial(hNB);
  int c, q;
  for (c = 0; c < hNB; c++) {
    Case *k = &cases[c];
    float *a = k->alr;
    R al[8] = {a[3], a[0], a[1], a[2], a[2], a[4], a[5], a[3]};
    if (k->relax > 0)
      for (q = 1; q < 7; q++)
        al[q] = k->relax;
    if (k->relaxp > 0)
      al[7] = k->relaxp;
    for (q = 0; q < 8; q++)
      alf[c * 8 + q] = al[q];
    cfac[c] = k->hyb;
    c1[c] = k->c1;
    visc[c] = k->visc;
    flowin[c] = k->flowin;
    flg[c] = k->flags;
    act[c] = 1;
  }
  for (q = 0; q < NA_; q++) {
    hs.a[q] = (R *)dalloc(nf);
    if (host[q])
      CK(cudaMemcpy(hs.a[q], host[q], nf, cudaMemcpyHostToDevice));
    else
      CK(cudaMemset(hs.a[q], 0, nf));
  }
  hs.ibc = (int *)upload(ibc, sizeof(int) * nw_);
  hs.jbc = (int *)upload(jbc, sizeof(int) * nw_);
  hs.kbc = (int *)upload(kbc, sizeof(int) * nw_);
  hs.ity = (int *)upload(ity, sizeof(int) * nw_);
  hs.nw = (int *)upload(nw, sizeof(int) * hNB);
  hs.yn = (R *)upload(hyn, sizeof(R) * nw_);
  hs.yn1 = (R *)upload(hyn1, sizeof(R) * nw_);
  hs.sx = (R *)upload(sxw, sizeof(R) * nw_);
  hs.sy = (R *)upload(syw, sizeof(R) * nw_);
  hs.sz = (R *)upload(szw, sizeof(R) * nw_);
  hs.anv = (R *)upload(zero, sizeof(R) * nw_);
  hs.anw = (R *)upload(zero, sizeof(R) * nw_);
  hs.ypl = (R *)upload(zero, sizeof(R) * nw_);
  hs.tau = (R *)upload(zero, sizeof(R) * nw_);
  hs.visc = (R *)upload(visc, sizeof(R) * hNB);
  hs.flowin = (R *)upload(flowin, sizeof(R) * hNB);
  hs.cfac = (R *)upload(cfac, sizeof(R) * hNB);
  hs.alf = (R *)upload(alf, sizeof(R) * hNB * 8);
  hs.flg = (int *)upload(flg, sizeof(int) * hNB);
  hs.c1 = (R *)upload(c1, sizeof(R) * hNB);
  hs.act = (int *)upload(act, sizeof(int) * hNB);
  hs.err = (R *)dalloc(sizeof(R) * hNB * 8);
  d_wk = (R *)dalloc(sizeof(R) * hNB * hDI);
  CK(cudaFuncSetAttribute(k_pblock, cudaFuncAttributeMaxDynamicSharedMemorySize, PBLOCK_SHARED));
  d_flow = (R *)dalloc(sizeof(R) * hNB);
  CK(cudaMemcpyToSymbol(G, &hs, sizeof(S)));
  SETI(DI, hDI);
  SETI(DJ, hDJ);
  SETI(NT, hNT);
  SETI(NB, hNB);
  SETI(MAXW, hMAXW);
  SETI(N, hN);
  SETI(LT, hLT);
  SETI(MT, hMT);
  SETI(NTT, hNTT);
  SETR(CMU, 0.09f);
  SETR(C2, 1.92f);
  SETR(CMU1, powf(0.09f, 0.25f));
  SETR(CMU2, powf(0.09f, 0.75f));
  SETR(EE, 9.01069f);
  SETR(CK, 0.4f);
  SETR(HINUM, 1.E30f);
  SETR(SMNUM, 1.E-30f);
  SETR(ALVIS, cases[0].alr[6]);
  free(cfac);
  free(c1);
  free(visc);
  free(flowin);
  free(alf);
  free(zero);
  free(flg);
  free(act);
}

enum { RUNNING, CONVERGED, DIVERGED };

static void run(int nlimt, R erext, int *iters, int *state) {
  cudaGraph_t graph;
  cudaGraphExec_t exec;
  R *err = ral((size_t)hNB * 8);
  int c, it, left = hNB, verb = getenv("VERB") != NULL, *act = ial(hNB);
  for (c = 0; c < hNB; c++)
    act[c] = 1;
  CK(cudaStreamCreate(&stream));
  CK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeThreadLocal));
  iterate();
  CK(cudaStreamEndCapture(stream, &graph));
  CK(cudaGraphInstantiate(&exec, graph, 0));
  for (it = 1; it <= nlimt && left > 0; it++) {
    int changed = 0;
    CK(cudaMemsetAsync(hs.err, 0, sizeof(R) * hNB * 8, stream));
    CK(cudaGraphLaunch(exec, stream));
    if (it % 10 && it != nlimt)
      continue;
    CK(cudaMemcpyAsync(err, hs.err, sizeof(R) * hNB * 8, cudaMemcpyDeviceToHost, stream));
    CK(cudaStreamSynchronize(stream));
    for (c = 0; c < hNB; c++)
      if (act[c]) {
        R *e = &err[c * 8], em = e[0] + e[1] + e[2] + e[3] + e[5] + e[6];
        iters[c] = it;
        if (!isfinite((double)em) || (it >= 20 && em > 1e3))
          state[c] = DIVERGED;
        else if (em <= erext)
          state[c] = CONVERGED;
        if (state[c]) {
          act[c] = 0;
          left--;
          changed = 1;
        }
        if (verb && c == 0)
          fprintf(stderr, "%d %g %g %g %g\n", it, (double)e[1], (double)e[2], (double)e[3],
                  (double)e[0]);
      }
    if (changed)
      CK(cudaMemcpy(hs.act, act, sizeof(int) * hNB, cudaMemcpyHostToDevice));
  }
  CK(cudaDeviceSynchronize());
  CK(cudaGraphExecDestroy(exec));
  CK(cudaGraphDestroy(graph));
  free(err);
  free(act);
}

static int write_case(int c, const char *suffix) {
  char p[1100];
  FILE *f;
  int i, j, k;
  snprintf(p, sizeof p, "%s/field%s", cases[c].dir, suffix);
  if (!(f = fopen(p, "w"))) {
    warn("%s: %s", p, strerror(errno));
    return 1;
  }
  int *b = cases[c].block;
  fprintf(f, "%5d%5d%5d%5d%5d%5d%5d%5d%5d\n", hL, hM, hN, b[0], b[1], b[2], b[3], b[4], b[5]);
  for (k = 1; k <= hN; k++)
    for (i = 1; i <= hL; i++) {
      fputs("\n\n", f);
      for (j = 1; j <= hM; j++) {
        size_t o = (size_t)c * hNT + ID(i, j, k);
        fprintf(f, "%8.4f%8.4f%8.4f%12.4E%12.4E%12.4E  %11.4E%11.4E%11.4E%11.4E%11.4E%11.4E\n",
                (double)host[X_][o], (double)host[Y_][o], (double)host[Z_][o], (double)host[U_][o],
                (double)host[V_][o], (double)host[W_][o], (double)host[P_][o], (double)host[TM_][o],
                (double)host[DK_][o], (double)host[DE_][o], (double)host[VISE_][o],
                (double)host[DEN_][o]);
      }
    }
  if (ferror(f) | fclose(f)) {
    warn("%s: write error", p);
    return 1;
  }
  return 0;
}

static void write_fields(const char *suffix) {
  int out[9] = {U_, V_, W_, P_, TM_, DK_, DE_, VISE_, DEN_};
  int c, q, bad = 0;
  for (q = 0; q < 9; q++)
    CK(cudaMemcpy(host[out[q]], hs.a[out[q]], sizeof(R) * hNB * hNT, cudaMemcpyDeviceToHost));
#pragma omp parallel for schedule(dynamic) reduction(+ : bad)
  for (c = 0; c < hNB; c++)
    bad += write_case(c, suffix);
  if (bad)
    die("%d of %d output files not written", bad, hNB);
}

static int envint(const char *name, int def) {
  char *s = getenv(name), *end;
  long v;
  if (!s)
    return def;
  v = strtol(s, &end, 10);
  if (*s == 0 || *end != 0 || v < 1 || v > 1000000000)
    die("%s=%s: not a positive integer", name, s);
  return (int)v;
}

static double envreal(const char *name, double def) {
  char *s = getenv(name), *end;
  double v;
  if (!s)
    return def;
  v = strtod(s, &end);
  if (*s == 0 || *end != 0 || !(v > 0))
    die("%s=%s: not a positive number", name, s);
  return v;
}

static int arg(const char *s, const char *name) {
  char *end;
  long v = strtol(s, &end, 10);
  if (*s == 0 || *end != 0 || v < 3 || v > 4096)
    die("%s %s: not an integer in [3, 4096]", name, s);
  return (int)v;
}

static void usage(void) {
  fprintf(stderr, "usage: %s di dj dk list suffix\n", argv0);
  exit(1);
}

int main(int argc, char **argv) {
  int c, nlimt, *iters, *state, times = getenv("TIMES") != NULL;
  double t0;
  R erext;
  if (argc != 6)
    usage();
  hDI = arg(argv[1], "di");
  hDJ = arg(argv[2], "dj");
  hDK = arg(argv[3], "dk");
  if (hDK > MAXDK)
    die("dk %d: at most %d", hDK, MAXDK);
  if (PBLOCK_SHARED > 227 * 1024)
    die("di %d: at most %d", hDI, (int)(227 * 1024 / (28 * sizeof(R))));
  hNT = hDI * hDJ * hDK;
  read_list(argv[4]);
  t0 = tnow();
  CK(cudaFree(0));
  if (times)
    fprintf(stderr, "init %.2f s\n", tnow() - t0);
  setup_all();
  if (times)
    fprintf(stderr, "setup %.2f s\n", tnow() - t0);
  upload_all();
  if (times)
    fprintf(stderr, "upload %.2f s\n", tnow() - t0);
  nlimt = envint("ITER", cases[0].nlimt);
  erext = envreal("TOL", 1.E-5f);
  iters = ial(hNB);
  state = ial(hNB);
  run(nlimt, erext, iters, state);
  if (times)
    fprintf(stderr, "iterate %.2f s\n", tnow() - t0);
  write_fields(argv[5]);
  if (times)
    fprintf(stderr, "write %.2f s\n", tnow() - t0);
  for (c = 0; c < hNB; c++)
    printf("%s iterations %d %s\n", cases[c].dir, iters[c],
           state[c] == CONVERGED  ? "done"
           : state[c] == DIVERGED ? "diverged"
                                  : "not converged");
  if (fflush(stdout) != 0)
    die("stdout: %s", strerror(errno));
  _exit(0);
}
