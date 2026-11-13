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
__constant__ int LT, MT;
__constant__ R CMU, C2, CMU1, CMU2, EE, CK, HINUM, SMNUM, ALVIS;

typedef struct {
  int i0, i1, j0, j1;
} Box;

enum {
  U_,
  V_,
  P_,
  PP_,
  DK_,
  DE_,
  TM_,
  VISE_,
  DEN_,
  DU_,
  DV_,
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
  X_,
  Y_,
  TJO_,
  CX_,
  CY_,
  EX_,
  EY_,
  TXXE_,
  TXXW_,
  TYYN_,
  TYYS_,
  TYXE_,
  TYXW_,
  TXYN_,
  TXYS_,
  F1_,
  LP_,
  LQ_,
  RU_,
  RV_,
  RK_,
  RE_,
  RM_,
  ZU_,
  ZV_,
  ZK_,
  ZE_,
  ZM_,
  OU_,
  OV_,
  OP_,
  OK_,
  OE_,
  UT_,
  VT_,
  NA_
};

typedef struct {
  R *a[NA_];
  int *ibc, *jbc, *ity, *nw;
  R *yn, *yn1, *sx, *sy, *anv, *ypl, *tau;
  R *visc, *flowin, *cfac, *err, *alf, *c1, *res;
  int *act, *flg;
} S;
__constant__ S G;

#define IDX(i, j) ((i) - 1 + DI * ((j) - 1))
#define A(n, i, j) G.a[n##_][(size_t)c * NT + IDX(i, j)]
#define CELL(B)                                                                                    \
  int tid = blockIdx.x * blockDim.x + threadIdx.x;                                                 \
  int ni = B.i1 - B.i0 + 1, nj = B.j1 - B.j0 + 1, nn = ni * nj;                                    \
  if (ni <= 0 || nj <= 0 || tid >= nn * NB)                                                        \
    return;                                                                                        \
  int c = tid / nn, r = tid % nn;                                                                  \
  int i = B.i0 + r % ni, j = B.j0 + r / ni
#define DECODE(B)                                                                                  \
  CELL(B);                                                                                         \
  if (!G.act[c])                                                                                   \
  return

__device__ __forceinline__ R *field(int which, int c) {
  return G.a[which] + (size_t)c * NT;
}
#define F(i, j) f[IDX(i, j)]

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
  A(F1, i, j) = A(VISE, i, j) / sigf;
}

__global__ void __launch_bounds__(256) k_coef(Box B, int fv) {
  DECODE(B);
  R *f = field(fv, c);
  R GAE = 0.5 * (A(F1, i + 1, j) + A(F1, i, j)), GAW = 0.5 * (A(F1, i - 1, j) + A(F1, i, j));
  R GAN = 0.5 * (A(F1, i, j + 1) + A(F1, i, j)), GAS = 0.5 * (A(F1, i, j - 1) + A(F1, i, j));
  R DENE = 0.5 * (A(DEN, i + 1, j) + A(DEN, i, j)), DENW = 0.5 * (A(DEN, i - 1, j) + A(DEN, i, j));
  R DENN = 0.5 * (A(DEN, i, j + 1) + A(DEN, i, j)), DENS = 0.5 * (A(DEN, i, j - 1) + A(DEN, i, j));
#define AV(v, di, dj) (0.5 * (A(v, i + di, j + dj) + A(v, i, j)))
  R UE = AV(U, 1, 0), UW = AV(U, -1, 0), UN = AV(U, 0, 1), US = AV(U, 0, -1);
  R VE = AV(V, 1, 0), VW = AV(V, -1, 0), VN = AV(V, 0, 1), VS = AV(V, 0, -1);
  R cx = A(CX, i, j), cy = A(CY, i, j), ex = A(EX, i, j), ey = A(EY, i, j), tj = A(TJO, i, j);
  R CE = DENE * (UE * cx + VE * cy), CW = DENW * (UW * cx + VW * cy);
  R CN = DENN * (UN * ex + VN * ey), CS = DENS * (US * ex + VS * ey);
  R TXXE = A(TXXE, i, j), TXXW = A(TXXW, i, j), TYYN = A(TYYN, i, j), TYYS = A(TYYS, i, j);
  R TYXE = A(TYXE, i, j), TYXW = A(TYXW, i, j), TXYN = A(TXYN, i, j), TXYS = A(TXYS, i, j);
  R DDE = GAE * TXXE + GAN * TXYN - GAS * TXYS;
  R DDW = GAW * TXXW + GAS * TXYS - GAN * TXYN;
  R DDN = GAN * TYYN + GAE * TYXE - GAW * TYXW;
  R DDS = GAS * TYYS + GAW * TYXW - GAE * TYXE;
  R CPO = fabs(CE - CW + CN - CS), cf = G.cfac[c];
  A(AE, i, j) = (fmax(fabs(cf * CE), DDE) - 0.5 * CE) * tj;
  A(AW, i, j) = (fmax(fabs(cf * CW), DDW) + 0.5 * CW) * tj;
  A(AN, i, j) = (fmax(fabs(cf * CN), DDN) - 0.5 * CN) * tj;
  A(AS, i, j) = (fmax(fabs(cf * CS), DDS) + 0.5 * CS) * tj;
  R DDNE = GAE * TYXE + GAN * TXYN, DDSE = -GAE * TYXE - GAS * TXYS,
    DDNW = -GAW * TYXW - GAN * TXYN, DDSW = GAW * TYXW + GAS * TXYS;
#define CROSS(Q)                                                                                   \
  (CPO * Q(i, j) + DDNE * Q(i + 1, j + 1) + DDSE * Q(i + 1, j - 1) + DDNW * Q(i - 1, j + 1) +       \
   DDSW * Q(i - 1, j - 1))
#define QF(a, b) F(a, b)
#define QV(a, b) A(V, a, b)
#define QK(a, b) A(DK, a, b)
  A(SU, i, j) = CROSS(QF);
  A(DV, i, j) = CROSS(QV);
  A(PP, i, j) = CROSS(QK);
  A(SUK, i, j) = CPO;
  A(SPK, i, j) = -CPO;
  R UX = (UE - UW) * cx + (UN - US) * ex, UY = (UE - UW) * cy + (UN - US) * ey;
  R VX = (VE - VW) * cx + (VN - VS) * ex, VY = (VE - VW) * cy + (VN - VS) * ey;
  A(GEN, i, j) = A(VISE, i, j) *
                 ((UY + VX) * (UY + VX) + ((G.flg[c] & 1) ? (R)0 : 2 * (UX * UX + VY * VY)));
}

__global__ void __launch_bounds__(256) k_src_u(Box B) {
  DECODE(B);
#define F1A(a, b) A(F1, a, b)
  R GAE = 0.5 * (F1A(i + 1, j) + F1A(i, j)), GAW = 0.5 * (F1A(i - 1, j) + F1A(i, j));
  R GAN = 0.5 * (F1A(i, j + 1) + F1A(i, j)), GAS = 0.5 * (F1A(i, j - 1) + F1A(i, j));
#define PA(a, b) A(P, a, b)
  R PCXI = (PA(i + 1, j) + PA(i + 1, j + 1) - PA(i, j) - PA(i, j + 1)) * 0.5;
  R PEDA = (PA(i, j + 1) + PA(i + 1, j + 1) - PA(i, j) - PA(i + 1, j)) * 0.5;
#define Q(q, a, b) q##p[(a) + (b) * DI]
#define GRADS(q)                                                                                   \
  const R *q##p = &A(q, i, j);                                                                     \
  R q##CE = Q(q, 1, 0) - Q(q, 0, 0), q##CW = Q(q, 0, 0) - Q(q, -1, 0);                             \
  R q##CN = (Q(q, 1, 1) + Q(q, 1, 0) - Q(q, -1, 1) - Q(q, -1, 0)) * 0.25;                          \
  R q##CS = (Q(q, 1, 0) + Q(q, 1, -1) - Q(q, -1, 0) - Q(q, -1, -1)) * 0.25;                        \
  R q##EE = (Q(q, 1, 1) + Q(q, 0, 1) - Q(q, 1, -1) - Q(q, 0, -1)) * 0.25;                          \
  R q##EW = (Q(q, 0, 1) + Q(q, -1, 1) - Q(q, 0, -1) - Q(q, -1, -1)) * 0.25;                        \
  R q##EN = Q(q, 0, 1) - Q(q, 0, 0), q##ES = Q(q, 0, 0) - Q(q, 0, -1);
  GRADS(U);
  GRADS(V);
  R cx = A(CX, i, j), cy = A(CY, i, j), ex = A(EX, i, j), ey = A(EY, i, j);
  R su = A(SU, i, j) - PCXI * cx - PEDA * ex;
  R dv = A(DV, i, j) - PCXI * cy - PEDA * ey;
#define FACE(c0)                                                                                   \
  R c0##E = (A(c0, i + 1, j) + A(c0, i, j)) * 0.5, c0##W = (A(c0, i - 1, j) + A(c0, i, j)) * 0.5, \
    c0##N = (A(c0, i, j + 1) + A(c0, i, j)) * 0.5, c0##S = (A(c0, i, j) + A(c0, i, j - 1)) * 0.5;
#define SOC(q, Ca, Ea, x, e)                                                                       \
  (x * (GAE * (q##CE * Ca##E + q##EE * Ea##E) - GAW * (q##CW * Ca##W + q##EW * Ea##W)) +           \
   e * (GAN * (q##CN * Ca##N + q##EN * Ea##N) - GAS * (q##CS * Ca##S + q##ES * Ea##S)))
  FACE(CX);
  FACE(EX);
  FACE(CY);
  FACE(EY);
  su += SOC(U, CX, EX, cx, ex) + SOC(V, CX, EX, cy, ey);
  dv += SOC(U, CY, EY, cx, ex) + SOC(V, CY, EY, cy, ey);
  R tj = A(TJO, i, j);
  A(SU, i, j) = su * tj;
  A(DV, i, j) = dv;
  A(SP, i, j) = A(SPK, i, j) * tj;
}
__global__ void k_src_v(Box B) {
  DECODE(B);
  A(SU, i, j) = A(DV, i, j) * A(TJO, i, j);
}
__global__ void k_src_k(Box B) {
  DECODE(B);
  R tj = A(TJO, i, j), d = A(DEN, i, j);
  A(SU, i, j) = (A(GEN, i, j) + A(PP, i, j)) * tj;
  A(SP, i, j) = (A(SPK, i, j) - CMU * d * d * A(DK, i, j) / A(VISE, i, j)) * tj;
}
__global__ void k_src_e(Box B) {
  DECODE(B);
  R tj = A(TJO, i, j), d = A(DEN, i, j);
  A(SU, i, j) =
      (A(SU, i, j) + G.c1[c] * CMU * A(GEN, i, j) * d * d * A(DK, i, j) / A(VISE, i, j)) * tj;
  A(SP, i, j) = (A(SPK, i, j) - C2 * d * A(DE, i, j) / (A(DK, i, j) + SMNUM)) * tj;
}

__device__ void wall_function(int ie, int c, int w, int i, int j, int ii, int jj, int i2, int j2,
                              R *f, R *an) {
  size_t o = (size_t)c * MAXW + w;
  R YP = G.yn[o], YP1 = G.yn1[o], DK1 = A(DK, i, j), DEN1 = A(DEN, i, j), TJ = A(TJO, i, j),
    visc = G.visc[c];
  R *SU1 = &A(SU, i, j), *SP1 = &A(SP, i, j);
  R F1v = f[IDX(ii, jj)];
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
    *an = P1 * *an;
  } else if (ie == 2) {
    R TM = -G.tau[o], PTA = TJ * G.sy[o] * TM / YP1;
    *SP1 -= PTA;
    *SU1 += PTA * F1v;
    *an = G.anv[o];
  } else if (ie == 5) {
    R SQ = sqrt(DK1), YPL = G.ypl[o];
    R DITM =
        YPL <= 11.63 ? DEN1 * CMU2 * SQ * YPL / YP : DEN1 * CMU2 * SQ * log(EE * YPL) / (CK * YP);
    R du = A(U, i, j) - A(U, ii, jj), dv = A(V, i, j) - A(V, ii, jj),
      P1 = (G.flg[c] & 2) ? (G.sx[o] * du) * (G.sx[o] * du) + (G.sy[o] * dv) * (G.sy[o] * dv)
                          : du * du + dv * dv;
    du = A(U, i2, j2) - A(U, ii, jj);
    dv = A(V, i2, j2) - A(V, ii, jj);
    R P2 = du * du + dv * dv;
    R YP2 = 2 * YP1 - YP, G1 = A(GEN, i, j), VI = A(VISE, i, j), GENR = G1 - VI * P2 / YP2 / YP2;
    if (GENR <= 0)
      GENR = 0;
    R SPK1 = A(SPK, i, j), tau = G.tau[o];
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
  if (w0 > 0 && G.ibc[base + w0 - 1] == G.ibc[base + w0] && G.jbc[base + w0 - 1] == G.jbc[base + w0])
    return;
  int lk[5] = {0, AN_, AS_, AE_, AW_};
  R *f = field(fv, c);
  for (int w = w0;
       w < G.nw[c] && G.ibc[base + w] == G.ibc[base + w0] && G.jbc[base + w] == G.jbc[base + w0];
       w++) {
    size_t o = base + w;
    int i = G.ibc[o], j = G.jbc[o], t = G.ity[o];
    int di = (t == 3) - (t == 4), dj = (t == 1) - (t == 2);
    wall_function(ie, c, w, i, j, i + di, j + dj, i - di, j - dj, f,
                  &G.a[lk[t]][(size_t)c * NT + IDX(i, j)]);
  }
}
__global__ void k_outlet_links(Box B) {
  DECODE(B);
  if (G.flg[c] & 16)
    return;
  A(AE, i, j) = 0;
}
__global__ void k_assemble(Box B, int fv, int ie) {
  DECODE(B);
  R *f = field(fv, c);
  R alf = G.alf[c * 8 + ie];
  A(F1, i, j) = F(i, j);
  R ANAB = A(AE, i, j) + A(AW, i, j) + A(AN, i, j) + A(AS, i, j);
  R ap = ANAB - A(SP, i, j), pduv = 1;
  if (A(MCF, i, j) >= 1) {
    ap = alf;
    A(AN, i, j) = A(AS, i, j) = A(AE, i, j) = A(AW, i, j) = 0;
    A(SU, i, j) = F(i, j);
    pduv = 0;
  }
  R P1 = 1.2 * ap;
  ap = ap / alf;
  A(AP, i, j) = ap;
  A(SU, i, j) += pduv * (1 - alf) * ap * F(i, j);
  if (ie == 1)
    A(DU, i, j) = A(TJO, i, j) * pduv / (P1 - ANAB);
}
__device__ void wmaxp(int c, R *a, R v) {
  unsigned m = __activemask();
  int c0 = __shfl_sync(m, c, __ffs(m) - 1);
  int same = __all_sync(m, c == c0);
  if (same) {
    for (int o = 16; o > 0; o >>= 1) {
      R w = __shfl_xor_sync(m, v, o);
      v = fmax(v, w);
    }
    if ((threadIdx.x & 31) == __ffs(m) - 1)
      atomicMaxD(a, v);
  } else
    atomicMaxD(a, v);
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
__global__ void k_jline(Box B, int fv, int color) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x, ni = B.i1 - B.i0 + 1;
  if (tid >= ni * NB)
    return;
  int c = tid / ni, i = B.i0 + tid % ni;
  if (!G.act[c] || (i & 1) != color)
    return;
  R *f = field(fv, c), pp = 0, qq = F(i, B.j0 - 1);
  for (int j = B.j0; j <= B.j1; j++) {
    R as = A(AS, i, j), t = (R)1 / (A(AP, i, j) - as * pp + (R)1e-30);
    pp = A(AN, i, j) * t;
    qq = (A(AE, i, j) * F(i + 1, j) + A(AW, i, j) * F(i - 1, j) + A(SU, i, j) + as * qq) * t;
    A(LP, i, j) = pp;
    A(LQ, i, j) = qq;
  }
  R x = F(i, B.j1 + 1);
  for (int j = B.j1; j >= B.j0; j--) {
    x = A(LP, i, j) * x + A(LQ, i, j);
    F(i, j) = x;
  }
}
__global__ void k_iline(Box B, int fv, int color) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x, nj = B.j1 - B.j0 + 1;
  if (tid >= nj * NB)
    return;
  int c = tid / nj, j = B.j0 + tid % nj;
  if (!G.act[c] || (j & 1) != color)
    return;
  R *f = field(fv, c), pp = 0, qq = F(B.i0 - 1, j);
  for (int i = B.i0; i <= B.i1; i++) {
    R aw = A(AW, i, j), t = (R)1 / (A(AP, i, j) - aw * pp + (R)1e-30);
    pp = A(AE, i, j) * t;
    qq = (A(AN, i, j) * F(i, j + 1) + A(AS, i, j) * F(i, j - 1) + A(SU, i, j) + aw * qq) * t;
    A(LP, i, j) = pp;
    A(LQ, i, j) = qq;
  }
  R x = F(B.i1 + 1, j);
  for (int i = B.i1; i >= B.i0; i--) {
    x = A(LP, i, j) * x + A(LQ, i, j);
    F(i, j) = x;
  }
}
__global__ void k_err(Box B, int fv, int slot) {
  CELL(B);
  R v = 0;
  if (G.act[c]) {
    R *f = field(fv, c);
    v = fabs(F(i, j) - A(F1, i, j));
    if (!isfinite(v))
      v = INFINITY;
  }
  wmax(c, slot, v);
}

__global__ void k_pclear(Box B) {
  DECODE(B);
  A(PP, i, j) = 0;
  A(SUK, i, j) = 0;
  A(SPK, i, j) = 0;
  A(SP, i, j) = 0;
}
__global__ void k_pflux(Box B, int dir) {
  DECODE(B);
  R De, Dw, Dn, Ds, Ue, Uw, Un, Us, Ve, Vw, Vn, Vs, cx, cy, ex, ey, tq;
  if (dir == 0) {
#define SX0(v, x)                                                                                  \
  x##e = A(v, i, j);                                                                               \
  x##w = A(v, i - 1, j);                                                                           \
  x##n = 0.25 * (A(v, i, j) + A(v, i - 1, j) + A(v, i, j + 1) + A(v, i - 1, j + 1));              \
  x##s = 0.25 * (A(v, i, j) + A(v, i - 1, j) + A(v, i, j - 1) + A(v, i - 1, j - 1))
    SX0(DEN, D);
    SX0(UT, U);
    SX0(VT, V);
    cx = 0.5 * (A(CX, i, j) + A(CX, i - 1, j));
    cy = 0.5 * (A(CY, i, j) + A(CY, i - 1, j));
    ex = 0.5 * (A(EX, i, j) + A(EX, i - 1, j));
    ey = 0.5 * (A(EY, i, j) + A(EY, i - 1, j));
    tq = 0.5 * (A(TJO, i - 1, j) + A(TJO, i, j));
  } else if (dir == 1) {
#define SX1(v, x)                                                                                  \
  x##n = A(v, i, j);                                                                               \
  x##s = A(v, i, j - 1);                                                                           \
  x##e = 0.25 * (A(v, i, j) + A(v, i, j - 1) + A(v, i + 1, j) + A(v, i + 1, j - 1));              \
  x##w = 0.25 * (A(v, i, j) + A(v, i, j - 1) + A(v, i - 1, j) + A(v, i - 1, j - 1))
    SX1(DEN, D);
    SX1(UT, U);
    SX1(VT, V);
    cx = 0.5 * (A(CX, i, j) + A(CX, i, j - 1));
    cy = 0.5 * (A(CY, i, j) + A(CY, i, j - 1));
    ex = 0.5 * (A(EX, i, j) + A(EX, i, j - 1));
    ey = 0.5 * (A(EY, i, j) + A(EY, i, j - 1));
    tq = 0.5 * (A(TJO, i, j - 1) + A(TJO, i, j));
  } else {
#define SX2(v, x)                                                                                  \
  x##e = 0.5 * (A(v, i, j) + A(v, i + 1, j));                                                      \
  x##w = 0.5 * (A(v, i, j) + A(v, i - 1, j));                                                      \
  x##n = 0.5 * (A(v, i, j) + A(v, i, j + 1));                                                      \
  x##s = 0.5 * (A(v, i, j) + A(v, i, j - 1))
    SX2(DEN, D);
    SX2(UT, U);
    SX2(VT, V);
    cx = A(CX, i, j);
    cy = A(CY, i, j);
    ex = A(EX, i, j);
    ey = A(EY, i, j);
    tq = A(TJO, i, j);
  }
  R CE = De * (Ue * cx + Ve * cy), CW = Dw * (Uw * cx + Vw * cy);
  R CN = Dn * (Un * ex + Vn * ey), CS = Ds * (Us * ex + Vs * ey);
  R val = -(CE - CW + CN - CS) * tq;
  if (dir == 0)
    A(SUK, i, j) = val;
  else if (dir == 1)
    A(SPK, i, j) = val;
  else
    A(SP, i, j) = val;
}
__global__ void k_rccopy(Box B) {
  DECODE(B);
  A(UT, i, j) = A(U, i, j);
  A(VT, i, j) = A(V, i, j);
}
__global__ void k_rcvel(Box B) {
  DECODE(B);
  if (A(MCF, i, j) >= 1)
    return;
  int II = i + 1, JJ = j + 1;
  R PE = (A(P, II, j) + A(P, II, JJ)) * 0.5, PW = (A(P, i, j) + A(P, i, JJ)) * 0.5;
  R PN = (A(P, i, JJ) + A(P, II, JJ)) * 0.5, PS = (A(P, i, j) + A(P, II, j)) * 0.5;
  R a = PE - PW, e = PN - PS, du = (R)0.2 * A(DU, i, j);
  A(UT, i, j) = A(U, i, j) + du * (a * A(CX, i, j) + e * A(EX, i, j));
  A(VT, i, j) = A(V, i, j) + du * (a * A(CY, i, j) + e * A(EY, i, j));
}
__global__ void k_rcsrc(Box B) {
  DECODE(B);
  if (A(MCF, i, j) >= 2)
    return;
  R p = A(P, i, j);
  A(SU, i, j) -= (R)0.2 * (A(AE, i, j) * (p - A(P, i + 1, j)) + A(AW, i, j) * (p - A(P, i - 1, j)) +
                           A(AN, i, j) * (p - A(P, i, j + 1)) + A(AS, i, j) * (p - A(P, i, j - 1)));
}
__global__ void k_pcoef(Box B) {
  DECODE(B);
  int II = i - 1, JJ = j - 1;
  R GAE = (A(DU, i, j) + A(DU, i, JJ)) * 0.5;
  R GAW = (A(DU, II, j) + A(DU, II, JJ)) * 0.5;
  R GAN = (A(DU, i, j) + A(DU, II, j)) * 0.5;
  R GAS = (A(DU, i, JJ) + A(DU, II, JJ)) * 0.5;
  R GAZ = (A(DU, i, j) + A(DU, II, j) + A(DU, i, JJ) + A(DU, II, JJ)) * 0.25;
#define T2(v) ((A(v, i, j) + A(v, II, JJ)) * 0.5)
  R TQ = T2(TJO);
  A(AE, i, j) = TQ * GAE * T2(TXXE);
  A(AW, i, j) = TQ * GAW * T2(TXXW);
  A(AN, i, j) = TQ * GAN * T2(TYYN);
  A(AS, i, j) = TQ * GAS * T2(TYYS);
  R P1 = GAN != 0, P2 = GAS != 0, P3 = GAE != 0, P4 = GAW != 0, PZ = GAZ != 0;
  R PQ = 0.5 * ((P1 + P2) * (P3 + P4 + 2 * PZ) + (P3 + P4) * (P1 + P2 + 2 * PZ) +
                2 * PZ * (P1 + P2 + P3 + P4));
  if (PQ == 0)
    PQ = 1;
  A(SU, i, j) = (2 * PZ *
                     (P1 * A(SUK, i, j) + P2 * A(SUK, i, j - 1) + P3 * A(SPK, i, j) +
                      P4 * A(SPK, i - 1, j)) +
                 P1 * P3 * A(SP, i, j) + P3 * P2 * A(SP, i, j - 1) +
                 P2 * P4 * A(SP, i - 1, j - 1) + P4 * P1 * A(SP, i - 1, j)) /
                PQ;
}
__global__ void k_pdiag(Box B) {
  DECODE(B);
  A(AP, i, j) = A(AE, i, j) + A(AW, i, j) + A(AN, i, j) + A(AS, i, j);
  if (A(MCF, i, j) >= 2) {
    A(AP, i, j) = 1;
    A(AE, i, j) = A(AW, i, j) = A(AN, i, j) = A(AS, i, j) = 0;
    A(SU, i, j) = 0;
  }
}
__global__ void k_pmean(Box B, R *flow, int acc) {
  int c = blockIdx.x;
  if (c >= NB || !G.act[c] || (!acc && !(G.flg[c] & 16)))
    return;
  int ni = B.i1 - B.i0 + 1, nn = ni * (B.j1 - B.j0 + 1);
  __shared__ double sm[2][512], sn[2][512];
  double s0 = 0, s1 = 0, n0 = 0, n1 = 0;
  for (int t = threadIdx.x; t < nn; t += blockDim.x) {
    int i = B.i0 + t % ni, j = B.j0 + t / ni;
    if (A(MCF, i, j) >= 2)
      continue;
    if ((i + j) & 1) {
      s1 += A(SU, i, j);
      n1 += 1;
    } else {
      s0 += A(SU, i, j);
      n0 += 1;
    }
  }
  sm[0][threadIdx.x] = s0;
  sm[1][threadIdx.x] = s1;
  sn[0][threadIdx.x] = n0;
  sn[1][threadIdx.x] = n1;
  __syncthreads();
  for (int o = blockDim.x / 2; o > 0; o >>= 1) {
    if (threadIdx.x < o)
      for (int q = 0; q < 2; q++) {
        sm[q][threadIdx.x] += sm[q][threadIdx.x + o];
        sn[q][threadIdx.x] += sn[q][threadIdx.x + o];
      }
    __syncthreads();
  }
  R n = sn[0][0] + sn[1][0], m = n > 0 ? (sm[0][0] + sm[1][0]) / n : 0;
  if (acc && threadIdx.x == 0 && !(G.flg[c] & 16))
    flow[2 * NB + c] += (sm[0][0] + sm[1][0]) / flow[NB + c];
  for (int t = threadIdx.x; t < nn; t += blockDim.x) {
    int i = B.i0 + t % ni, j = B.j0 + t / ni;
    if (A(MCF, i, j) < 2)
      A(SU, i, j) -= m;
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
    for (int j = B.j0 + part; j <= B.j1; j += np_) {
      sa -= A(AE, i, j);
      sb -= A(AW, i, j);
      sc += A(SU, i, j);
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
  if (j > B.j0)
    F(i, j) = wk[(size_t)c * DI + i];
}
__global__ void k_pupdate(Box B) {
  CELL(B);
  R v = 0;
  if (G.act[c] && A(MCF, i, j) < 2) {
    R *f = field(PP_, c);
    R p1 = F(i, j) - f[IDX(LT, MT / 2)];
    A(P, i, j) += G.alf[c * 8 + 7] * p1;
    v = fabs(p1);
    if (!isfinite(v))
      v = INFINITY;
  }
  wmax(c, 0, v);
}
__global__ void k_correct(Box B) {
  DECODE(B);
  if (A(MCF, i, j) >= 1)
    return;
  R *f = field(PP_, c);
  int II = i + 1, JJ = j + 1;
  R PE = (F(II, j) + F(II, JJ)) * 0.5, PW = (F(i, j) + F(i, JJ)) * 0.5;
  R PN = (F(i, JJ) + F(II, JJ)) * 0.5, PS = (F(i, j) + F(II, j)) * 0.5;
  R a = PE - PW, e = PN - PS, du = A(DU, i, j);
  A(U, i, j) -= du * (a * A(CX, i, j) + e * A(EX, i, j));
  A(V, i, j) -= du * (a * A(CY, i, j) + e * A(EY, i, j));
}
__global__ void k_positive(Box B) {
  DECODE(B);
  if (A(DK, i, j) <= SMNUM)
    A(DK, i, j) = SMNUM;
}
__global__ void k_viscosity(Box B) {
  DECODE(B);
  R *dk = &A(DK, i, j), *de = &A(DE, i, j);
  if (*dk <= SMNUM)
    *dk = SMNUM;
  if (*de <= SMNUM)
    *de = SMNUM;
  R tv = *de <= SMNUM ? G.visc[c] : A(DEN, i, j) * CMU * *dk * *dk / *de + G.visc[c];
  A(VISE, i, j) += ALVIS * (tv - A(VISE, i, j));
}
__global__ void k_outflow(R *flow, int reset) {
  int c = blockIdx.x;
  if (c >= NB || !G.act[c])
    return;
  int i = LT;
  double fl = 0, ar = 0;
  for (int j = 2 + threadIdx.x; j <= MT + 1; j += blockDim.x) {
    R uc = (A(U, i, j) + A(U, i, j - 1)) * 0.5f, dc = (A(DEN, i, j) + A(DEN, i, j - 1)) * 0.5f;
    R p1 = A(X, i, j) - A(X, i, j - 1), p2 = A(Y, i, j) - A(Y, i, j - 1);
    R area = sqrt(p1 * p1 + p2 * p2);
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
  if (threadIdx.x == 0) {
    if (reset) {
      G.flowin[c] = s1[0];
      flow[2 * NB + c] = 0;
    }
    flow[NB + c] = s2[0];
    flow[c] = (s1[0] - G.flowin[c]) / s2[0] - flow[2 * NB + c];
  }
}
__global__ void k_outlet(Box B, R *uinc) {
  DECODE(B);
  if (G.flg[c] & 16)
    return;
  int I = LT;
  A(U, I + 1, j) = A(U, I, j) - uinc[c];
  A(V, I + 1, j) = A(V, I, j);
  A(DK, I + 1, j) = A(DK, I, j);
  A(DE, I + 1, j) = A(DE, I, j);
}

__global__ void k_resid(Box B, int fv, int rv, int q, int norm) {
  CELL(B);
  R res = 0, den = 0;
  if (G.act[c]) {
    R *f = field(fv, c);
    if (A(MCF, i, j) < 1 && A(SP, i, j) > (R)-1e20) {
      R ae = A(AE, i, j), aw = A(AW, i, j), an = A(AN, i, j), as = A(AS, i, j);
      R ap = ae + aw + an + as - A(SP, i, j);
      res = A(SU, i, j) + ae * F(i + 1, j) + aw * F(i - 1, j) + an * F(i, j + 1) +
            as * F(i, j - 1) - ap * F(i, j);
      den = fabs(ap * F(i, j));
    }
    field(rv, c)[IDX(i, j)] = res;
  }
  if (norm) {
    wmaxp(c, &G.res[c * 16 + 2 * q], isfinite(res) ? fabs(res) : INFINITY);
    wmaxp(c, &G.res[c * 16 + 2 * q + 1], den);
  }
}
__global__ void k_mres(Box B, int norm) {
  CELL(B);
  R res = 0, den = 0;
  if (G.act[c]) {
    res = A(MCF, i, j) < 2 ? A(SU, i, j) : 0;
    A(RM, i, j) = res;
    R u = A(U, i, j), v = A(V, i, j);
    den = (fabs(u * A(CX, i, j) + v * A(CY, i, j)) + fabs(u * A(EX, i, j) + v * A(EY, i, j))) *
          A(TJO, i, j) * A(DEN, i, j);
  }
  if (norm) {
    wmaxp(c, &G.res[c * 16 + 8], isfinite(res) ? fabs(res) : INFINITY);
    wmaxp(c, &G.res[c * 16 + 9], den);
  }
}
__global__ void k_addtau(Box B, int zv, int lim) {
  DECODE(B);
  if (A(MCF, i, j) < lim && (lim == 2 || A(SP, i, j) > (R)-1e20))
    A(SU, i, j) += field(zv, c)[IDX(i, j)];
}
__global__ void k_tau(Box B) {
  DECODE(B);
  A(ZU, i, j) -= A(RU, i, j);
  A(ZV, i, j) -= A(RV, i, j);
  A(ZK, i, j) -= A(RK, i, j);
  A(ZE, i, j) -= A(RE, i, j);
  A(ZM, i, j) -= A(RM, i, j);
}
typedef struct {
  int DI, NT, L, M, si, sj, oi, oj;
} Dim;
#define AT(s, d, n, i, j) s.a[n##_][(size_t)c * d.NT + (i) - 1 + d.DI * ((j) - 1)]
__global__ void k_down(S f, Dim df, S g, Dim dc) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x, n = dc.L * dc.M;
  if (tid >= n * NB)
    return;
  int c = tid / n, r = tid % n, I = 1 + r % dc.L, J = 1 + r / dc.L;
  int li = dc.oi && I == dc.L, lj = dc.oj && J == dc.M;
  int i = li ? df.L : dc.si * (I - 1) + 1, j = lj ? df.M : dc.sj * (J - 1) + 1;
  int na = li ? 3 : dc.si, nb = lj ? 3 : dc.sj;
  if (!f.act[c])
    return;
  AT(g, dc, OU, I, J) = AT(g, dc, U, I, J) = AT(f, df, U, i, j);
  AT(g, dc, OV, I, J) = AT(g, dc, V, I, J) = AT(f, df, V, i, j);
  AT(g, dc, OK, I, J) = AT(g, dc, DK, I, J) = AT(f, df, DK, i, j);
  AT(g, dc, OE, I, J) = AT(g, dc, DE, I, J) = AT(f, df, DE, i, j);
  AT(g, dc, VISE, I, J) = AT(f, df, VISE, i, j);
  if (I >= 2 && J >= 2) {
    R p = 0, m = 0;
    for (int b = 0; b < nb; b++)
      for (int a = 0; a < na; a++) {
        p += AT(f, df, P, i - na + 1 + a, j - nb + 1 + b);
        m += AT(f, df, RM, i - na + 1 + a, j - nb + 1 + b);
      }
    AT(g, dc, OP, I, J) = AT(g, dc, P, I, J) = p / (na * nb);
    AT(g, dc, ZM, I, J) = m;
  }
  if (I >= 2 && I < dc.L && J >= 2 && J < dc.M) {
    R su = 0, sv = 0, sk = 0, se = 0;
    for (int b = 1 - dc.sj; b < dc.sj; b++)
      for (int a = 1 - dc.si; a < dc.si; a++) {
        R w = (a ? (R)0.5 : (R)1) * (b ? (R)0.5 : (R)1);
        su += w * AT(f, df, RU, i + a, j + b);
        sv += w * AT(f, df, RV, i + a, j + b);
        sk += w * AT(f, df, RK, i + a, j + b);
        se += w * AT(f, df, RE, i + a, j + b);
      }
    AT(g, dc, ZU, I, J) = su;
    AT(g, dc, ZV, I, J) = sv;
    AT(g, dc, ZK, I, J) = sk;
    AT(g, dc, ZE, I, J) = se;
  }
}
__global__ void k_up(S f, Dim df, S g, Dim dc, R w, int turb) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x, n = df.L * df.M;
  if (tid >= n * NB)
    return;
  int c = tid / n, r = tid % n, i = 1 + r % df.L, j = 1 + r / df.L;
  if (!f.act[c])
    return;
  if (i >= 2 && i < df.L && j >= 2 && j < df.M && AT(f, df, MCF, i, j) < 1) {
    int I0 = dc.si == 2 ? (i + 1) / 2 : i, J0 = dc.sj == 2 ? (j + 1) / 2 : j;
    int I1 = I0 + (dc.si == 2 && !(i & 1)), J1 = J0 + (dc.sj == 2 && !(j & 1));
    R wa = I1 == I0 ? (R)1 : (R)0.5, wb = J1 == J0 ? (R)1 : (R)0.5;
    if (dc.oi && i > df.L - 3) {
      I0 = dc.L - 1;
      I1 = dc.L;
      wa = (R)(df.L - i) / 3;
    }
    if (dc.oj && j > df.M - 3) {
      J0 = dc.M - 1;
      J1 = dc.M;
      wb = (R)(df.M - j) / 3;
    }
#define DL(q, Q, a, b) (AT(g, dc, q, a, b) - AT(g, dc, Q, a, b))
#define DV(q, Q, a, b) (a > 1 && a < dc.L && b > 1 && b < dc.M ? DL(q, Q, a, b) : (R)0)
#define BIL(q, Q)                                                                                  \
  (wa * wb * DV(q, Q, I0, J0) + (1 - wa) * wb * DV(q, Q, I1, J0) +                                 \
   wa * (1 - wb) * DV(q, Q, I0, J1) + (1 - wa) * (1 - wb) * DV(q, Q, I1, J1))
    AT(f, df, U, i, j) += w * BIL(U, OU);
    AT(f, df, V, i, j) += w * BIL(V, OV);
    if (turb) {
      R k0 = AT(f, df, DK, i, j), e0 = AT(f, df, DE, i, j);
      AT(f, df, DK, i, j) = fmax(k0 + w * BIL(DK, OK), (R)0.5 * k0);
      AT(f, df, DE, i, j) = fmax(e0 + w * BIL(DE, OE), (R)0.5 * e0);
    }
  }
  if (i >= 2 && j >= 2 && AT(f, df, MCF, i, j) < 2) {
    int K = dc.si == 2 ? (i + 2) / 2 : i, L = dc.sj == 2 ? (j + 2) / 2 : j;
    int Kn = dc.si == 1 ? K : (i & 1) ? K + 1 : K - 1;
    int Ln = dc.sj == 1 ? L : (j & 1) ? L + 1 : L - 1;
    R wa = dc.si == 2 ? (R)0.75 : (R)1, wb = dc.sj == 2 ? (R)0.75 : (R)1;
    if (dc.oi && i > df.L - 3)
      K = Kn = dc.L;
    if (dc.oj && j > df.M - 3)
      L = Ln = dc.M;
    if (Kn < 2 || Kn > dc.L)
      Kn = K;
    if (Ln < 2 || Ln > dc.M)
      Ln = L;
#define DP(a, b) (AT(g, dc, MCF, a, b) < 2 ? DL(P, OP, a, b) : d00)
    R d00 = DL(P, OP, K, L), d10 = DP(Kn, L), d01 = DP(K, Ln), d11 = DP(Kn, Ln);
    AT(f, df, P, i, j) +=
        w * (wa * wb * d00 + (1 - wa) * wb * d10 + wa * (1 - wb) * d01 + (1 - wa) * (1 - wb) * d11);
  }
}

static const char *argv0 = "ns2d";

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
  int i, j, t;
  float yn, yn1, sx, sy;
} Wall;

typedef struct {
  char dir[1024];
  R hyb, c1;
  double relax, relaxp;
  int flags, block[6];
  int l, m, n, ig, nlimt, nswp;
  float visc, flowin, alr[8];
  float *z;
  Wall *walls;
  int nwall;
} Case;

static int hDI, hDJ, hNT, hNB, hMAXW;
static int hL, hM, hN, hLT, hMT, hIG, curl;
static Case *cases;
static S hs;
static R *d_wk, *d_flow;
static int pnp = 8;
#define PBLOCK_SHARED (sizeof(R) * (4 + 3 * pnp) * hDI)
static cudaStream_t stream;
static int nsweep[8];
static R sigma[8];

static int nblk(long n) {
  return (int)((n + 255) / 256);
}
static long cnt(Box b) {
  long n = (long)(b.i1 - b.i0 + 1) * (b.j1 - b.j0 + 1);
  return n > 0 ? n * hNB : 0;
}
#define RUN(kern, b, ...)                                                                          \
  do {                                                                                             \
    long n_ = cnt(b);                                                                              \
    if (n_)                                                                                        \
      kern<<<nblk(n_), 256, 0, stream>>>(b, ##__VA_ARGS__);                                        \
  } while (0)

static void sweeps(Box b, int fv, int nsw) {
  int ni = b.i1 - b.i0 + 1, nj = b.j1 - b.j0 + 1;
  for (int s = 0; s < nsw; s++) {
    for (int color = 0; color < 2; color++)
      k_jline<<<nblk((long)ni * hNB), 256, 0, stream>>>(b, fv, color);
    for (int color = 0; color < 2; color++)
      k_iline<<<nblk((long)nj * hNB), 256, 0, stream>>>(b, fv, color);
  }
}

static void prep(int ie, int fv) {
  Box in = {2, hLT, 2, hMT};
  Box all = {1, hL, 1, hM};
  Box outlet = {hLT, hLT, 2, hMT};
  if (ie == 1 || ie == 6) {
    RUN(k_diffusivity, all, sigma[ie]);
    RUN(k_coef, in, fv);
  }
  if (ie == 1)
    RUN(k_src_u, in);
  else if (ie == 2)
    RUN(k_src_v, in);
  else if (ie == 5)
    RUN(k_src_k, in);
  else if (ie == 6)
    RUN(k_src_e, in);
  if (hIG == 2)
    k_wall<<<nblk((long)hNB * hMAXW), 256, 0, stream>>>(ie, fv);
  RUN(k_outlet_links, outlet);
}
static const int zof[7] = {0, ZU_, ZV_, 0, 0, ZK_, ZE_}, rof[7] = {0, RU_, RV_, 0, 0, RK_, RE_};
static void solveq(int ie, int fv, int tau) {
  Box in = {2, hLT, 2, hMT};
  prep(ie, fv);
  if (tau)
    RUN(k_addtau, in, zof[ie], 1);
  RUN(k_assemble, in, fv, ie);
  sweeps(in, fv, nsweep[ie]);
  RUN(k_err, in, fv, ie);
}
static void pmass(int tau) {
  Box in = {2, hL, 2, hM}, all = {1, hL, 1, hM}, cells = {2, hLT, 2, hMT};
  RUN(k_rccopy, all);
  RUN(k_rcvel, cells);
  RUN(k_pclear, in);
  RUN(k_pflux, in, 0);
  RUN(k_pflux, in, 1);
  RUN(k_pflux, in, 2);
  RUN(k_pcoef, in);
  RUN(k_rcsrc, in);
  if (tau)
    RUN(k_addtau, in, ZM_, 2);
}
static void solvep(int tau) {
  Box in = {2, hL, 2, hM};
  Box cells = {2, hLT, 2, hMT};
  pmass(tau);
  RUN(k_pdiag, in);
  k_pmean<<<hNB, 512, 0, stream>>>(in, d_flow, 1);
  k_pblock<<<hNB, 64 * pnp, PBLOCK_SHARED, stream>>>(in, d_wk);
  RUN(k_pblock_set, in, d_wk);
  sweeps(in, PP_, nsweep[0]);
  RUN(k_pupdate, in);
  RUN(k_correct, cells);
}

static void iterate(int tau) {
  Box outlet = {1, 1, 2, hMT};
  Box all = {1, hL, 1, hM};
  k_outflow<<<hNB, 256, 0, stream>>>(d_flow, 0);
  RUN(k_outlet, outlet, d_flow);
  solveq(1, U_, tau);
  solveq(2, V_, tau);
  if (hIG == 2 && !curl) {
    solveq(5, DK_, tau);
    RUN(k_positive, all);
    solveq(6, DE_, tau);
  }
  solvep(tau);
  if (hIG == 2 && !curl)
    RUN(k_viscosity, all);
}

static void resid(int tau, int norm) {
  Box in = {2, hLT, 2, hMT}, pin = {2, hL, 2, hM}, outlet = {1, 1, 2, hMT};
  int ie;
  k_outflow<<<hNB, 256, 0, stream>>>(d_flow, 0);
  RUN(k_outlet, outlet, d_flow);
  for (ie = 1; ie <= 6; ie++) {
    int fv = ie == 1 ? U_ : ie == 2 ? V_ : ie == 5 ? DK_ : DE_;
    if (ie == 3 || ie == 4 || (ie >= 5 && (hIG != 2 || curl)))
      continue;
    prep(ie, fv);
    if (tau)
      RUN(k_addtau, in, zof[ie], 1);
    RUN(k_resid, in, fv, rof[ie], ie == 1 ? 0 : ie == 2 ? 1 : ie == 5 ? 2 : 3, norm);
    if (ie == 1)
      RUN(k_assemble, in, fv, ie);
  }
  pmass(tau);
  k_pmean<<<hNB, 512, 0, stream>>>(pin, d_flow, 0);
  RUN(k_mres, pin, norm);
}

typedef struct {
  int DI, DJ, NT, L, M, LT, MT, MAXW;
  S hs, *dG;
  R *wk, *flow, *host[NA_];
  int *dims, blk[4], si, sj, oi, oj;
} Level;
#define MAXLV 12
static Level lv[MAXLV];
static int nlv = 1, nu1, nu2, nuc, mgmin;
static R mgw;
static void *symp[7];

static Dim dim(int l) {
  Dim d = {lv[l].DI, lv[l].NT, lv[l].L, lv[l].M, lv[l].si, lv[l].sj, lv[l].oi, lv[l].oj};
  return d;
}
static void use(int l) {
  Level *v = &lv[l];
  int k;
  curl = l;
  hDI = v->DI;
  hDJ = v->DJ;
  hNT = v->NT;
  hL = v->L;
  hM = v->M;
  hLT = v->LT;
  hMT = v->MT;
  hMAXW = v->MAXW;
  hs = v->hs;
  d_wk = v->wk;
  d_flow = v->flow;
  for (k = 0; k < 6; k++)
    CK(cudaMemcpyAsync(symp[k], v->dims + k, sizeof(int), cudaMemcpyDeviceToDevice, stream));
  CK(cudaMemcpyAsync(symp[6], v->dG, sizeof(S), cudaMemcpyDeviceToDevice, stream));
}
static void cycle(int l) {
  int k;
  use(l);
  if (l == nlv - 1) {
    for (k = 0; k < nuc; k++)
      iterate(l > 0);
    return;
  }
  for (k = 0; k < nu1; k++)
    iterate(l > 0);
  resid(l > 0, 0);
  Dim df = dim(l), dc = dim(l + 1);
  k_down<<<nblk((long)dc.L * dc.M * hNB), 256, 0, stream>>>(lv[l].hs, df, lv[l + 1].hs, dc);
  use(l + 1);
  k_outflow<<<hNB, 256, 0, stream>>>(d_flow, 1);
  resid(0, 0);
  {
    Box all = {1, hL, 1, hM};
    RUN(k_tau, all);
  }
  cycle(l + 1);
  use(l);
  k_up<<<nblk((long)df.L * df.M * hNB), 256, 0, stream>>>(lv[l].hs, df, lv[l + 1].hs, dc, mgw, 0);
  for (k = 0; k < nu2; k++) {
    if (l == 0 && k == nu2 - 1)
      CK(cudaMemsetAsync(hs.err, 0, sizeof(R) * hNB * 8, stream));
    iterate(l > 0);
  }
}

#define ID(i, j) ((i) - 1 + hDI * ((j) - 1))

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

static const int stored[] = {X_,    Y_,    TJO_,  CX_,   CY_,   EX_,   EY_,  TXXE_,
                             TXXW_, TYYN_, TYYS_, TYXE_, TYXW_, TXYN_, TXYS_, MCF_,
                             U_,    V_,    P_,    DK_,   DE_,   VISE_, DEN_,  TM_};
#define NSTORED ((int)(sizeof stored / sizeof stored[0]))

static R *host[NA_];
static int *ibc, *jbc, *ity, *nw;
static R *hyn, *hyn1, *sxw, *syw;
static Wall *wl;
static int wcap;

static float *g[NA_], *vm[5];
static float *X, *Y;
static int *MC, L1, L2, M1, M2, nwall;

static void walval(int is, int it, int js, int jt, float *F) {
  int i, j, w;
  for (j = js; j <= jt; j++) {
    F[ID(is - 1, j)] = F[ID(is, j)];
    F[ID(it + 1, j)] = F[ID(it, j)];
  }
  for (i = is - 1; i <= it + 1; i++) {
    F[ID(i, js - 1)] = F[ID(i, js)];
    F[ID(i, jt + 1)] = F[ID(i, jt)];
  }
  for (w = 0; w < nwall; w++) {
    int i = wl[w].i, j = wl[w].j, o = ID(i, j);
    switch (wl[w].t) {
    case 1:
      F[ID(i, j + 1)] = F[o];
      break;
    case 2:
      F[ID(i, j - 1)] = F[o];
      break;
    case 3:
      F[ID(i + 1, j)] = F[o];
      break;
    case 4:
      F[ID(i - 1, j)] = F[o];
      break;
    }
  }
}

static void wall(int i, int j, int t, float p1, float p2, int a1, int b1, int a2, int b2) {
  float pq = sqrtf(p1 * p1 + p2 * p2), q1, q2, aa, bb, cc = 1.0f, coth, yn;
  int o = ID(i, j), o1 = ID(a1, b1), o2 = ID(a2, b2), n;
  p1 = p1 / pq;
  p2 = p2 / pq;
  if (nwall == wcap) {
    wcap = wcap ? 2 * wcap : 256;
    wl = (Wall *)erealloc(wl, wcap * sizeof(Wall));
  }
  n = nwall++;
  wl[n].i = i;
  wl[n].j = j;
  wl[n].t = t;
  wl[n].sx = sqrtf(1.f - p1 * p1);
  wl[n].sy = sqrtf(1.f - p2 * p2);
  q1 = X[o] - X[o1];
  q2 = Y[o] - Y[o1];
  aa = sqrtf((q1 - p1) * (q1 - p1) + (q2 - p2) * (q2 - p2));
  bb = sqrtf(q1 * q1 + q2 * q2);
  coth = (bb * bb + cc * cc - aa * aa) / (2 * bb * cc);
  yn = bb * fabsf(coth);
  wl[n].yn = yn;
  q1 = X[o] - X[o2];
  q2 = Y[o] - Y[o2];
  bb = sqrtf(q1 * q1 + q2 * q2);
  aa = sqrtf((q1 - p1) * (q1 - p1) + (q2 - p2) * (q2 - p2));
  coth = (bb * bb + cc * cc - aa * aa) / (2 * bb * cc);
  wl[n].yn1 = (bb * fabsf(coth) + yn) * 0.5f;
}

static void wallj(int i, int j, int t, int j1, int j2) {
  int i1 = i + 1, i2 = i - 1;
  if (t == 1) {
    if (i == L1)
      i2 = i;
    if (i == L2)
      i1 = i;
  }
  wall(i, j, t, -(Y[ID(i1, j)] - Y[ID(i2, j)]), X[ID(i1, j)] - X[ID(i2, j)], i, j1, i, j2);
}
static void walli(int i, int j, int t, int i1, int i2) {
  int j1 = j + 1, j2 = j - 1;
  if (j == M1)
    j2 = j;
  if (j == M2)
    j1 = j;
  wall(i, j, t, -(Y[ID(i, j1)] - Y[ID(i, j2)]), X[ID(i, j1)] - X[ID(i, j2)], i1, j, i2, j);
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

static void geom(int c, int L, int M);
static void setup(int c) {
  Case *K = &cases[c];
  char s[512], p[1100];
  FILE *ff;
  int i, j, k, q, NN = hNT, L, M, N, LO, MO;
  float visc;
  float *U = g[U_], *V = g[V_], *DK = g[DK_], *DE = g[DE_], *VISE = g[VISE_], *DEN = g[DEN_];
  read_control(K);
  visc = K->visc;
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
  for (q = 0; q < 6; q++)
    K->block[q] = ix(s, 15 + 5 * q, 5);
  K->l = L;
  K->m = M;
  K->n = N;
  LO = L + 1;
  MO = M + 1;
  if (L < 3 || M < 3 || N < 2)
    die("%s: grid %dx%dx%d is too small", p, L, M, N);
  if (LO > hDI || MO > hDJ)
    die("%s: grid %dx%d does not fit %dx%d", p, L, M, hDI, hDJ);
  for (q = 0; q < NSTORED; q++)
    memset(g[stored[q]], 0, sizeof(float) * NN);
  for (q = 0; q < 5; q++)
    memset(vm[q], 0, sizeof(float) * NN);
  memset(MC, 0, sizeof(int) * NN);
  K->z = (float *)emalloc(sizeof(float) * (size_t)N * L * M);
  X = g[X_];
  Y = g[Y_];
  for (i = 1; i <= LO; i++)
    for (j = 1; j <= MO; j++) {
      DEN[ID(i, j)] = 1.0f;
      VISE[ID(i, j)] = visc;
    }
  for (k = 1; k <= N; k++)
    for (i = 1; i <= L; i++)
      for (j = 1; j <= M; j++) {
        int o = ID(i, j);
        do {
          egets(s, sizeof s, ff, p);
        } while (strlen(s) < 100);
        K->z[((size_t)(k - 1) * L + (i - 1)) * M + (j - 1)] = fx(s, 16, 8);
        if (k != 2)
          continue;
        X[o] = fx(s, 0, 8);
        Y[o] = fx(s, 8, 8);
        U[o] = fx(s, 24, 12);
        V[o] = fx(s, 36, 12);
        g[P_][o] = fx(s, 62, 11);
        g[TM_][o] = fx(s, 73, 11);
        DK[o] = fx(s, 84, 11);
        DE[o] = fx(s, 95, 11);
      }
  fclose(ff);
  geom(c, L, M);
}

static void geom(int c, int L, int M) {
  Case *K = &cases[c];
  int i, j, q, n, NN = hNT, LO = L + 1, MO = M + 1, LT = L - 1, MT = M - 1, wdist;
  float visc = K->visc, alvis = K->alr[6], cmu = 0.09f, smnum = 1.E-30f, flowin = 0.0f;
  float *U = g[U_], *V = g[V_], *DK = g[DK_], *DE = g[DE_], *VISE = g[VISE_], *DEN = g[DEN_];
  float *cell[5] = {g[CX_], g[CY_], g[EX_], g[EY_], g[TJO_]};
  float *CX = g[CX_], *CY = g[CY_], *EX = g[EX_], *EY = g[EY_];
  X = g[X_];
  Y = g[Y_];
  wl = NULL;
  wcap = 0;
  for (q = 0; q < 5; q++)
    memset(vm[q], 0, sizeof(float) * NN);
  memset(MC, 0, sizeof(int) * NN);
  for (i = 1; i <= LO; i++)
    for (j = 1; j <= MO; j++) {
      int m = (j == 1 || (j == M && !(K->flags & 16))) ? 1 : 0;
      if (i >= L1 && i <= L2 && j >= M1 && j <= M2)
        m = 1;
      if (i > L1 && i <= L2 && j > M1 && j <= M2)
        m = 2;
      MC[ID(i, j)] = m;
    }

  nwall = 0;
  wdist = (K->flags & 4) != 0;
  for (i = 2; i <= LT; i++)
    for (j = 2; j <= MT; j++) {
      if (MC[ID(i, j)] != 0)
        continue;
      if (MC[ID(i + 1, j)] + MC[ID(i - 1, j)] + MC[ID(i, j + 1)] + MC[ID(i, j - 1)] == 0)
        continue;
      if (MC[ID(i, j + 1)])
        wallj(i, j, 1, wdist ? j + 1 : j - 1, wdist ? j - 1 : j - 2);
      if (MC[ID(i, j - 1)])
        wallj(i, j, 2, wdist ? j - 1 : j + 1, wdist ? j + 1 : j + 2);
      if (MC[ID(i + 1, j)])
        walli(i, j, 3, wdist ? i + 1 : i - 1, wdist ? i - 1 : i - 2);
      if (MC[ID(i - 1, j)] && !((K->flags & 8) && j == M2))
        walli(i, j, 4, wdist ? i - 1 : i + 1, wdist ? i + 1 : i + 2);
    }

  for (i = 1; i <= L; i++)
    for (j = 1; j <= M; j++) {
      int o = ID(i, j);
      if (MC[o]) {
        DK[o] = 0;
        DE[o] = 0;
        U[o] = 0;
        V[o] = 0;
      }
    }

  for (i = 2; i <= L; i++)
    for (j = 2; j <= M; j++) {
      int o = ID(i, j), oi = ID(i - 1, j), oj = ID(i, j - 1), oij = ID(i - 1, j - 1);
      float p1 = (X[o] + X[oj] - X[oi] - X[oij]) * 0.5f, p2 = (X[o] + X[oi] - X[oj] - X[oij]) * 0.5f;
      float q1 = (Y[o] + Y[oj] - Y[oi] - Y[oij]) * 0.5f, q2 = (Y[o] + Y[oi] - Y[oj] - Y[oij]) * 0.5f;
      float ptr = 1.0f / (p1 * q2 - p2 * q1);
      vm[0][o] = ptr * q2;
      vm[1][o] = -ptr * p2;
      vm[2][o] = -ptr * q1;
      vm[3][o] = ptr * p1;
      vm[4][o] = fabsf(1.0f / ptr);
    }
  for (q = 0; q < 5; q++)
    walval(2, L, 2, M, vm[q]);
  for (i = 1; i <= L; i++)
    for (j = 1; j <= M; j++)
      for (q = 0; q < 5; q++) {
        float *A = vm[q];
        cell[q][ID(i, j)] =
            (A[ID(i, j)] + A[ID(i, j + 1)] + A[ID(i + 1, j)] + A[ID(i + 1, j + 1)]) * 0.25f;
      }

  for (i = 2; i <= LT; i++)
    for (j = 2; j <= MT; j++) {
      int o = ID(i, j), e_ = ID(i + 1, j), w_ = ID(i - 1, j), n_ = ID(i, j + 1), s_ = ID(i, j - 1);
      float cx = CX[o], cy = CY[o], ex = EX[o], ey = EY[o];
      g[TXXE_][o] = WH2(CX, e_) * cx + WH2(CY, e_) * cy;
      g[TXXW_][o] = WH2(CX, w_) * cx + WH2(CY, w_) * cy;
      g[TYYN_][o] = WH2(EX, n_) * ex + WH2(EY, n_) * ey;
      g[TYYS_][o] = WH2(EX, s_) * ex + WH2(EY, s_) * ey;
      g[TYXE_][o] = (WH2(EX, e_) * cx + WH2(EY, e_) * cy) * 0.25f;
      g[TYXW_][o] = (WH2(EX, w_) * cx + WH2(EY, w_) * cy) * 0.25f;
      g[TXYN_][o] = (WH2(CX, n_) * ex + WH2(CY, n_) * ey) * 0.25f;
      g[TXYS_][o] = (WH2(CX, s_) * ex + WH2(CY, s_) * ey) * 0.25f;
    }
  for (q = TXXE_; q <= TXYS_; q++)
    walval(2, LT, 2, MT, g[q]);

  if (K->ig == 2)
    for (i = 1; i <= L; i++)
      for (j = 1; j <= M; j++) {
        int o = ID(i, j);
        float tv;
        if (DK[o] <= smnum)
          DK[o] = smnum;
        if (DE[o] <= smnum)
          DE[o] = smnum;
        tv = DE[o] <= smnum ? visc : DEN[o] * cmu * DK[o] * DK[o] / DE[o] + visc;
        VISE[o] = VISE[o] + alvis * (tv - VISE[o]);
      }

  for (j = 2; j <= M; j++) {
    int a = ID(1, j), b = ID(1, j - 1);
    float uc = (U[a] + U[b]) * 0.5f, denc = (DEN[a] + DEN[b]) * 0.5f;
    float p1 = X[a] - X[b], p2 = Y[a] - Y[b];
    flowin = flowin + denc * sqrtf(p1 * p1 + p2 * p2) * uc;
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

static void pack_walls(void);
static void setup_all(void) {
  int c, q;
  for (q = 0; q < NSTORED; q++) {
    host[stored[q]] = ral((size_t)hNB * hNT);
    g[stored[q]] = fal(hNT);
  }
  for (q = 0; q < 5; q++)
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
  hL = cases[0].l;
  hM = cases[0].m;
  hN = cases[0].n;
  hLT = hL - 1;
  hMT = hM - 1;
  hIG = cases[0].ig;
  nsweep[0] = cases[0].nswp;
  nsweep[1] = nsweep[2] = 7;
  nsweep[5] = nsweep[6] = 5;
  sigma[1] = 1.0f;
  sigma[5] = 1.0f;
  sigma[6] = 1.3f;
  pack_walls();
}

static void pack_walls(void) {
  size_t nw_;
  int c, n;
  hMAXW = 1;
  for (c = 0; c < hNB; c++)
    if (cases[c].nwall > hMAXW)
      hMAXW = cases[c].nwall;
  nw_ = (size_t)hNB * hMAXW;
  ibc = ial(nw_);
  jbc = ial(nw_);
  ity = ial(nw_);
  nw = ial(hNB);
  hyn = ral(nw_);
  hyn1 = ral(nw_);
  sxw = ral(nw_);
  syw = ral(nw_);
  for (c = 0; c < hNB; c++) {
    Case *k = &cases[c];
    for (n = 0; n < k->nwall; n++) {
      size_t o = (size_t)c * hMAXW + n;
      Wall *w = &k->walls[n];
      ibc[o] = w->i;
      jbc[o] = w->j;
      ity[o] = w->t;
      hyn[o] = w->yn;
      hyn1[o] = w->yn1;
      sxw[o] = w->sx;
      syw[o] = w->sy;
    }
    nw[c] = k->nwall;
    free(k->walls);
    k->walls = NULL;
  }
}

static void save(int l) {
  Level *v = &lv[l];
  int d[6] = {hDI, hDJ, hNT, hMAXW, hLT, hMT}, q;
  v->DI = hDI;
  v->DJ = hDJ;
  v->NT = hNT;
  v->L = hL;
  v->M = hM;
  v->LT = hLT;
  v->MT = hMT;
  v->MAXW = hMAXW;
  v->hs = hs;
  v->wk = d_wk;
  v->flow = d_flow;
  for (q = 0; q < NA_; q++)
    v->host[q] = host[q];
  v->dims = (int *)upload(d, sizeof d);
  v->dG = (S *)upload(&hs, sizeof(S));
}

static int ratio(int n, int b0, int b1, int *odd) {
  *odd = (n - 1) % 2;
  if ((n - 1) / 2 + 1 < mgmin || (b0 && (b0 - 1) % 2) || (b1 && (b1 - 1) % 2) ||
      (b0 && b1 - b0 < 2)) {
    *odd = 0;
    return 1;
  }
  return 2;
}
static int coarsenable(int l) {
  Level *v = &lv[l];
  lv[l + 1].si = ratio(v->L, v->blk[0], v->blk[1], &lv[l + 1].oi);
  lv[l + 1].sj = ratio(v->M, v->blk[2], v->blk[3], &lv[l + 1].oj);
  return lv[l + 1].si * lv[l + 1].sj > 1;
}

static void upload_all(void);
static void build(int l) {
  Level *p = &lv[l - 1];
  int c, q, i, j, I, J;
  static const int inj[] = {X_, Y_, U_, V_, DK_, DE_, TM_};
  int si = lv[l].si, sj = lv[l].sj;
  for (q = 0; q < 4; q++) {
    int r = q < 2 ? si : sj;
    lv[l].blk[q] = p->blk[q] ? (p->blk[q] - 1) / r + 1 : 0;
  }
  hL = (p->L - 1) / si + 1;
  hM = (p->M - 1) / sj + 1;
  hLT = hL - 1;
  hMT = hM - 1;
  hDI = hL + 1;
  hDJ = hM + 1;
  hNT = hDI * hDJ;
  for (q = 0; q < NA_; q++)
    host[q] = NULL;
  for (q = 0; q < NSTORED; q++)
    host[stored[q]] = ral((size_t)hNB * hNT);
  L1 = lv[l].blk[0];
  L2 = lv[l].blk[1];
  M1 = lv[l].blk[2];
  M2 = lv[l].blk[3];
  for (c = 0; c < hNB; c++) {
    size_t po = (size_t)c * p->NT;
    for (q = 0; q < NSTORED; q++)
      memset(g[stored[q]], 0, sizeof(float) * hNT);
    for (I = 1; I <= hL + 1; I++)
      for (J = 1; J <= hM + 1; J++) {
        g[DEN_][ID(I, J)] = 1.0f;
        g[VISE_][ID(I, J)] = cases[c].visc;
      }
    for (I = 1; I <= hL; I++)
      for (J = 1; J <= hM; J++) {
        int li = lv[l].oi && I == hL, lj = lv[l].oj && J == hM;
        int na = li ? 3 : si, nb = lj ? 3 : sj;
        i = li ? p->L : si * (I - 1) + 1;
        j = lj ? p->M : sj * (J - 1) + 1;
        for (q = 0; q < (int)(sizeof inj / sizeof inj[0]); q++)
          g[inj[q]][ID(I, J)] = p->host[inj[q]][po + (i - 1) + p->DI * (j - 1)];
        if (I >= 2 && J >= 2) {
          R s = 0;
          int a, b;
          for (a = 0; a < na; a++)
            for (b = 0; b < nb; b++)
              s += p->host[P_][po + (i - na + a) + p->DI * (j - nb + b)];
          g[P_][ID(I, J)] = s / (na * nb);
        }
      }
    geom(c, hL, hM);
  }
  pack_walls();
  upload_all();
  save(l);
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
  hs.ity = (int *)upload(ity, sizeof(int) * nw_);
  hs.nw = (int *)upload(nw, sizeof(int) * hNB);
  hs.yn = (R *)upload(hyn, sizeof(R) * nw_);
  hs.yn1 = (R *)upload(hyn1, sizeof(R) * nw_);
  hs.sx = (R *)upload(sxw, sizeof(R) * nw_);
  hs.sy = (R *)upload(syw, sizeof(R) * nw_);
  hs.anv = (R *)upload(zero, sizeof(R) * nw_);
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
  hs.res = (R *)dalloc(sizeof(R) * hNB * 16);
  d_wk = (R *)dalloc(sizeof(R) * hNB * hDI);
  CK(cudaFuncSetAttribute(k_pblock, cudaFuncAttributeMaxDynamicSharedMemorySize, PBLOCK_SHARED));
  d_flow = (R *)dalloc(sizeof(R) * 3 * hNB);
  CK(cudaMemset(d_flow, 0, sizeof(R) * 3 * hNB));
  CK(cudaMemcpyToSymbol(G, &hs, sizeof(S)));
  SETI(DI, hDI);
  SETI(DJ, hDJ);
  SETI(NT, hNT);
  SETI(NB, hNB);
  SETI(MAXW, hMAXW);
  SETI(LT, hLT);
  SETI(MT, hMT);
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
  R *err = ral((size_t)hNB * 8), *res = ral((size_t)hNB * 16);
  int k;
  int c, it, left = hNB, verb = getenv("VERB") != NULL, *act = ial(hNB);
  for (c = 0; c < hNB; c++)
    act[c] = 1;
  CK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeThreadLocal));
  if (nlv > 1)
    cycle(0);
  else
    for (k = 0; k < 10; k++) {
      if (k == 9)
        CK(cudaMemsetAsync(hs.err, 0, sizeof(R) * hNB * 8, stream));
      iterate(0);
    }
  CK(cudaMemsetAsync(hs.res, 0, sizeof(R) * hNB * 16, stream));
  resid(0, 1);
  CK(cudaStreamEndCapture(stream, &graph));
  CK(cudaGraphInstantiate(&exec, graph, 0));
  int per = nlv > 1 ? nu1 + nu2 : 10, nc = (nlimt + per - 1) / per;
  for (it = 1; it <= nc && left > 0; it++) {
    int changed = 0;
    CK(cudaMemsetAsync(lv[0].hs.err, 0, sizeof(R) * hNB * 8, stream));
    CK(cudaGraphLaunch(exec, stream));
    CK(cudaMemcpyAsync(err, lv[0].hs.err, sizeof(R) * hNB * 8, cudaMemcpyDeviceToHost, stream));
    CK(cudaMemcpyAsync(res, lv[0].hs.res, sizeof(R) * hNB * 16, cudaMemcpyDeviceToHost, stream));
    CK(cudaStreamSynchronize(stream));
    for (c = 0; c < hNB; c++)
      if (act[c]) {
        R *e = &err[c * 8], em = e[0] + e[1] + e[2] + e[5] + e[6], rm = 0, *r = &res[c * 16];
        for (k = 0; k < 5; k++) {
          R q = r[2 * k + 1] > 0 ? r[2 * k] / r[2 * k + 1] : r[2 * k] > 0 ? INFINITY : 0;
          if (!(q <= rm))
            rm = q;
        }
        iters[c] = it * per;
        if (!isfinite((double)em) || (it >= 20 && em > 1e3))
          state[c] = DIVERGED;
        else if (rm <= erext)
          state[c] = CONVERGED;
        if (state[c]) {
          act[c] = 0;
          left--;
          changed = 1;
        }
        if (verb && c == 0)
          fprintf(stderr, "%d %g %g %g %g %g %g\n", it * per, (double)e[1], (double)e[2],
                  (double)e[5], (double)e[6], (double)e[0], (double)rm);
      }
    if (changed)
      for (int l = 0; l < nlv; l++)
        CK(cudaMemcpy(lv[l].hs.act, act, sizeof(int) * hNB, cudaMemcpyHostToDevice));
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
        size_t o = (size_t)c * hNT + ID(i, j);
        fprintf(f, "%8.4f%8.4f%8.4f%12.4E%12.4E%12.4E  %11.4E%11.4E%11.4E%11.4E%11.4E%11.4E\n",
                (double)host[X_][o], (double)host[Y_][o],
                (double)cases[c].z[((size_t)(k - 1) * hL + (i - 1)) * hM + (j - 1)],
                (double)host[U_][o], (double)host[V_][o], 0.0, (double)host[P_][o],
                (double)host[TM_][o], (double)host[DK_][o], (double)host[DE_][o],
                (double)host[VISE_][o], (double)host[DEN_][o]);
      }
    }
  if (ferror(f) | fclose(f)) {
    warn("%s: write error", p);
    return 1;
  }
  return 0;
}

static void write_fields(const char *suffix) {
  int out[8] = {U_, V_, P_, TM_, DK_, DE_, VISE_, DEN_};
  int c, q, bad = 0;
  for (q = 0; q < 8; q++)
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
  fprintf(stderr, "usage: %s di dj list suffix\n", argv0);
  exit(1);
}

int main(int argc, char **argv) {
  int c, nlimt, *iters, *state, times = getenv("TIMES") != NULL;
  double t0;
  R erext;
  if (argc != 5)
    usage();
  hDI = arg(argv[1], "di");
  hDJ = arg(argv[2], "dj");
  pnp = (int)((227 * 1024 / (sizeof(R) * hDI) - 4) / 3);
  if (pnp > 8)
    pnp = 8;
  if (pnp < 1)
    die("di %d: at most %d", hDI, (int)(227 * 1024 / (7 * sizeof(R))));
  hNT = hDI * hDJ;
  read_list(argv[3]);
  t0 = tnow();
  CK(cudaFree(0));
  if (times)
    fprintf(stderr, "init %.2f s\n", tnow() - t0);
  setup_all();
  if (times)
    fprintf(stderr, "setup %.2f s\n", tnow() - t0);
  upload_all();
  save(0);
  for (c = 0; c < 4; c++)
    lv[0].blk[c] = cases[0].block[c];
  for (c = 1; c < hNB; c++)
    for (int q = 0; q < 4; q++)
      if (cases[c].block[q] != cases[0].block[q])
        lv[0].blk[0] = -1;
  {
    int mg = envint("MG", MAXLV);
    mgmin = envint("MGMIN", 33);
    nu1 = envint("NU1", 3);
    nu2 = envint("NU2", 3);
    nuc = envint("NUC", 10);
    mgw = envreal("MGW", 0.7);
    while (nlv < mg && nlv < MAXLV && lv[0].blk[0] >= 0 && coarsenable(nlv - 1)) {
      build(nlv);
      nlv++;
    }
  }
  CK(cudaFuncSetAttribute(k_pblock, cudaFuncAttributeMaxDynamicSharedMemorySize,
                          (int)(sizeof(R) * (4 + 3 * pnp) * lv[0].DI)));
  CK(cudaGetSymbolAddress(&symp[0], DI));
  CK(cudaGetSymbolAddress(&symp[1], DJ));
  CK(cudaGetSymbolAddress(&symp[2], NT));
  CK(cudaGetSymbolAddress(&symp[3], MAXW));
  CK(cudaGetSymbolAddress(&symp[4], LT));
  CK(cudaGetSymbolAddress(&symp[5], MT));
  CK(cudaGetSymbolAddress(&symp[6], G));
  CK(cudaStreamCreate(&stream));
  use(0);
  CK(cudaStreamSynchronize(stream));
  for (c = 0; c < NA_; c++)
    host[c] = lv[0].host[c];
  if (times || getenv("VERB"))
    for (c = 0; c < nlv; c++)
      fprintf(stderr, "level %d: %dx%d\n", c, lv[c].L, lv[c].M);
  if (times)
    fprintf(stderr, "upload %.2f s\n", tnow() - t0);
  nlimt = envint("ITER", cases[0].nlimt);
  erext = envreal("TOL", 1.E-5f);
  iters = ial(hNB);
  state = ial(hNB);
  run(nlimt, erext, iters, state);
  if (times)
    fprintf(stderr, "iterate %.2f s\n", tnow() - t0);
  write_fields(argv[4]);
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
