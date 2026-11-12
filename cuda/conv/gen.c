#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

enum { MAXN = 4096 };

static void die(const char *m) {
  fprintf(stderr, "gen: %s\n", m);
  exit(1);
}

static double sgn(double a) { return (a > 0) - (a < 0); }

static int refine(const double *q, int n, int r, double *o) {
  double m[MAXN], d[MAXN];
  int k, i, no = (n - 1) * r + 1;
  if (no > MAXN)
    die("grid too large");
  for (k = 0; k < n - 1; k++)
    m[k] = q[k + 1] - q[k];
  for (k = 1; k < n - 1; k++)
    d[k] = (sgn(m[k - 1]) != sgn(m[k]) || m[k - 1] == 0 || m[k] == 0) ? 0 : 6 / (3 / m[k - 1] + 3 / m[k]);
  d[0] = (3 * m[0] - m[1]) / 2;
  if (sgn(d[0]) != sgn(m[0]))
    d[0] = 0;
  else if (sgn(m[0]) != sgn(m[1]) && fabs(d[0]) > 3 * fabs(m[0]))
    d[0] = 3 * m[0];
  d[n - 1] = (3 * m[n - 2] - m[n - 3]) / 2;
  if (sgn(d[n - 1]) != sgn(m[n - 2]))
    d[n - 1] = 0;
  else if (sgn(m[n - 2]) != sgn(m[n - 3]) && fabs(d[n - 1]) > 3 * fabs(m[n - 2]))
    d[n - 1] = 3 * m[n - 2];
  for (i = 0; i < no; i++) {
    k = i / r;
    if (k == n - 1) {
      o[i] = q[n - 1];
      continue;
    }
    double t = (double)(i % r) / r, t2 = t * t, t3 = t2 * t;
    o[i] = (2 * t3 - 3 * t2 + 1) * q[k] + (t3 - 2 * t2 + t) * d[k] + (-2 * t3 + 3 * t2) * q[k + 1] + (t3 - t2) * d[k + 1];
  }
  return no;
}

static int refine_keep(const double *q, int n, int r, int lo, int hi, double *o) {
  int a = lo ? 1 : 0, b = hi ? n - 2 : n - 1, no = 0, i;
  for (i = 0; i < a; i++)
    o[no++] = q[i];
  no += refine(q + a, b - a + 1, r, o + no);
  for (i = b + 1; i < n; i++)
    o[no++] = q[i];
  return no;
}

static int half(int r, double *o) {
  double y1 = 0.024, q = pow(1.2, 1.0 / r), lo = 1, hi, t;
  int n = (int)ceil(log(1 + 0.5 * (q - 1) / y1) / log(q)), i, it;
  hi = q;
  for (it = 0; it < 200; it++) {
    t = (lo + hi) / 2;
    if (y1 * (pow(t, n) - 1) / (t - 1) > 0.5)
      hi = t;
    else
      lo = t;
  }
  o[0] = 0;
  for (i = 1; i <= n; i++)
    o[i] = o[i - 1] + y1 * pow(lo, i - 1);
  o[n] = 0.5;
  return n;
}

static void section(int r, double *rr, int *M, double *z, int *N) {
  double h[MAXN];
  int n = half(r, h), j;
  for (j = 0; j <= n; j++) {
    rr[j] = 2.8 - h[j];
    rr[2 * n - j] = 1.8 + h[j];
  }
  *M = 2 * n + 1;
  for (j = 0; j <= n; j++)
    z[j] = 0.5 - h[n - j];
  *N = n + 1;
}

static FILE *out;
static int RI, RJ, RK;

static void start(const char *d, const int *h, const char *ctl) {
  char p[1024];
  FILE *f;
  int q;
  mkdir(d, 0755);
  snprintf(p, sizeof p, "%s/control", d);
  if (!(f = fopen(p, "w")))
    die(p);
  fputs(ctl, f);
  fclose(f);
  snprintf(p, sizeof p, "%s/field", d);
  if (!(out = fopen(p, "w")))
    die(p);
  for (q = 0; q < 9; q++)
    fprintf(out, "%5d", h[q]);
  fputc('\n', out);
}

static void row(double x, double y, double z, double u, double v, double w, double tk, double te) {
  fprintf(out, "%8.5f%8.5f%8.5f%12.4E%12.4E%12.4E  %11.4E%11.4E%11.4E%11.4E%11.4E%11.4E\n", x, y, z, u, v, w, 0.0, 0.0, tk, te,
          0.0, 0.0);
}

static void finish(void) {
  if (ferror(out) | fclose(out))
    die("write error");
}

static double prof(double a, double b) {
  double y = a - 0.5, v = 0;
  int n;
  for (n = 1; n < 201; n += 2)
    v += (((n - 1) / 2) % 2 ? -1 : 1) * (1 - cosh(n * M_PI * b) / cosh(n * M_PI * 0.5)) * cos(n * M_PI * y) / ((double)n * n * n);
  return v;
}

static void bend(const char *d, int turb, int keep) {
  static const double up[] = {-2.7447, -1.8373, -1.031, -0.3901}, dn[] = {0.1865, 0.4345, 0.739, 1.0906, 1.4892, 1.9318},
                      rb[] = {0.0, 0.0308, 0.0725, 0.1259, 0.1871, 0.2614, 0.3363, 0.4163, 0.4978,
                              0.5819, 0.6614, 0.7411, 0.8118, 0.8757, 0.9273, 0.9683, 1.0},
                      zb[] = {0.0, 0.0588, 0.183, 0.3268, 0.4902, 0.6601, 0.8105, 0.9346, 1.0};
  double rc = 2.3, sb[20], q[17], s[MAXN], rr[MAXN], z[MAXN], bulk = 0;
  int L, M, N, i, j, k, n = 0;
  for (i = 0; i < 4; i++)
    sb[n++] = up[i];
  for (i = 0; i <= 9; i++)
    sb[n++] = rc * (10.0 * i * M_PI / 180);
  for (i = 0; i < 6; i++)
    sb[n++] = rc * M_PI / 2 + dn[i];
  L = refine(sb, 20, RI, s);
  for (j = 0; j < 17; j++)
    q[j] = rb[16 - j];
  M = keep ? refine_keep(q, 17, RJ, 1, 1, rr) : refine(q, 17, RJ, rr);
  for (j = 0; j < M; j++)
    rr[j] += 1.8;
  for (k = 0; k < 9; k++)
    q[k] = zb[k] * 0.5;
  N = keep ? refine_keep(q, 9, RK, 0, 1, z) : refine(q, 9, RK, z);
  for (i = 0; i < 401; i++)
    for (j = 0; j < 401; j++)
      bulk += prof(j / 400.0, -0.5 + i / 400.0);
  bulk /= 401.0 * 401.0;
  int h[9] = {L, M, N, 0, 0, 0, 0, 0, 0};
  start(d, h,
        turb ? " 3000    2   10    1\n 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.7000 0.5000\n40000.0 0.0000\n"
             : " 3000    1   10    1\n 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000\n  790.0 0.0000\n");
  for (k = 0; k < N; k++)
    for (i = 0; i < L; i++) {
      fputs("\n\n", out);
      for (j = 0; j < M; j++) {
        double v = (j > 0 && j < M - 1 && k < N - 1) ? prof(2.8 - rr[j], z[k]) / bulk : 0, x, y, tx, ty;
        if (s[i] <= 0)
          x = s[i], y = -rr[j], tx = 1, ty = 0;
        else if (s[i] <= rc * M_PI / 2 + 1e-9) {
          double t = s[i] / rc;
          x = rr[j] * sin(t), y = -rr[j] * cos(t), tx = cos(t), ty = sin(t);
        } else
          x = rr[j], y = s[i] - rc * M_PI / 2, tx = 0, ty = 1;
        row(x, y, z[k], v * tx, v * ty, 0, turb ? 0.005 : 0, turb ? 0.001162 : 0);
      }
    }
  finish();
}

static const char *TCTL = " 3000    2   10    1\n 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.7000 0.5000\n40000.0 0.0000\n";

static void tduct(const char *d, int r) {
  double rr[MAXN], z[MAXN];
  int M, N, L = 121, i, j, k;
  section(r, rr, &M, z, &N);
  int h[9] = {L, M, N, 0, 0, 0, 0, 0, 0};
  start(d, h, TCTL);
  for (k = 0; k < N; k++)
    for (i = 0; i < L; i++) {
      fputs("\n\n", out);
      for (j = 0; j < M; j++) {
        int wall = j == 0 || j == M - 1 || k == N - 1;
        double c = rr[j] - 2.3, q = wall ? 0 : M_SQRT1_2;
        row((0.5 * i + c) * M_SQRT1_2, (0.5 * i - c) * M_SQRT1_2, z[k], q, q, 0, 0.005, 0.001162);
      }
    }
  finish();
}

static double inl[4][MAXN];

static void inlet(const char *f, int M, int N) {
  char s[512];
  int L, m, n, i, j, k;
  FILE *p = fopen(f, "r");
  if (!p || !fgets(s, sizeof s, p) || sscanf(s, "%5d%5d%5d", &L, &m, &n) != 3)
    die(f);
  if (m != M || n != N)
    die("inlet grid differs");
  for (k = 0; k < n; k++)
    for (i = 0; i < L; i++)
      for (j = 0; j < m;) {
        if (!fgets(s, sizeof s, p))
          die(f);
        if (strlen(s) <= 100)
          continue;
        if (i == L - 2) {
          char b[16];
          double v[4];
          int c, at[4] = {24, 36, 84, 95}, w[4] = {12, 12, 11, 11};
          for (c = 0; c < 4; c++) {
            memcpy(b, s + at[c], w[c]);
            b[w[c]] = 0;
            v[c] = atof(b);
          }
          inl[0][k * m + j] = (v[0] + v[1]) * M_SQRT1_2;
          inl[1][k * m + j] = (v[1] - v[0]) * M_SQRT1_2;
          inl[2][k * m + j] = v[2];
          inl[3][k * m + j] = v[3];
        }
        j++;
      }
  fclose(p);
}

static int leg(double ds, double len, double *o) {
  double q = 0, h = ds;
  int n = 0;
  o[n++] = 0;
  while (q + h < len) {
    q += h;
    o[n++] = q;
    h = h * 1.1 < 0.3 ? h * 1.1 : 0.3;
  }
  if (len - q < 0.5 * h && n > 1)
    n--;
  o[n++] = len;
  return n;
}

static int designed(int nb, double rc, double *s) {
  double ds = rc * M_PI / 2 / (9 * nb), u[MAXN], w[MAXN];
  int nu = leg(ds, 2.7447, u), nw = leg(ds, 1.9318, w), L = 0, i;
  if (nu + nw + 9 * nb > MAXN)
    die("grid too large");
  for (i = nu - 1; i > 0; i--)
    s[L++] = -u[i];
  for (i = 0; i <= 9 * nb; i++)
    s[L++] = rc * M_PI / 2 * i / (9 * nb);
  for (i = 1; i < nw; i++)
    s[L++] = rc * M_PI / 2 + w[i];
  return L;
}

static void tbendg(const char *d, int r, int rs, int nb, const char *f) {
  static const double up[] = {-2.7447, -1.8373, -1.031, -0.3901}, dn[] = {0.1865, 0.4345, 0.739, 1.0906, 1.4892, 1.9318};
  double rc = 2.3, sb[20], s[MAXN], rr[MAXN], z[MAXN];
  int L, M, N, i, j, k, n = 0;
  for (i = 0; i < 4; i++)
    sb[n++] = up[i];
  for (i = 0; i <= 9; i++)
    sb[n++] = rc * (10.0 * i * M_PI / 180);
  for (i = 0; i < 6; i++)
    sb[n++] = rc * M_PI / 2 + dn[i];
  L = nb > 0 ? designed(nb, rc, s) : refine(sb, 20, rs, s);
  section(r, rr, &M, z, &N);
  if (M * N > MAXN)
    die("grid too large");
  inlet(f, M, N);
  int h[9] = {L, M, N, 0, 0, 0, 0, 0, 0};
  start(d, h, TCTL);
  for (k = 0; k < N; k++)
    for (i = 0; i < L; i++) {
      fputs("\n\n", out);
      for (j = 0; j < M; j++) {
        double v = inl[0][k * M + j], x, y, tx, ty;
        if (s[i] <= 0)
          x = s[i], y = -rr[j], tx = 1, ty = 0;
        else if (s[i] <= rc * M_PI / 2 + 1e-9) {
          double t = s[i] / rc;
          x = rr[j] * sin(t), y = -rr[j] * cos(t), tx = cos(t), ty = sin(t);
        } else
          x = rr[j], y = s[i] - rc * M_PI / 2, tx = 0, ty = 1;
        row(x, y, z[k], v * tx, v * ty, 0, inl[2][k * M + j], inl[3][k * M + j]);
      }
    }
  finish();
}

static void cavity(const char *d, const char *nodes) {
  int n0 = 31, n = (n0 - 1) * RI + 1, i, j, k;
  double *X0 = malloc(sizeof(double) * n0 * n0), *Y0 = malloc(sizeof(double) * n0 * n0);
  double *Xa = malloc(sizeof(double) * n * n0), *Ya = malloc(sizeof(double) * n * n0);
  double *X = malloc(sizeof(double) * n * n), *Y = malloc(sizeof(double) * n * n);
  double a[MAXN], b[MAXN], dz = 1.0 / (n0 - 1);
  if (!X0 || !Y0 || !Xa || !Ya || !X || !Y)
    die("out of memory");
  for (i = 0; i < n0; i++)
    for (j = 0; j < n0; j++) {
      X0[i * n0 + j] = (double)i / (n0 - 1);
      Y0[i * n0 + j] = (double)j / (n0 - 1);
    }
  if (nodes) {
    FILE *f = fopen(nodes, "r");
    char ln[256];
    double px, py, x, y;
    int ii, jj;
    if (!f)
      die(nodes);
    if (!fgets(ln, sizeof ln, f))
      die(nodes);
    while (fgets(ln, sizeof ln, f))
      if (sscanf(ln, "%d,%d,%lf,%lf,%lf,%lf", &ii, &jj, &px, &py, &x, &y) == 6 && ii >= 1 && ii <= n0 && jj >= 1 && jj <= n0) {
        X0[(ii - 1) * n0 + jj - 1] = x;
        Y0[(ii - 1) * n0 + jj - 1] = y;
      }
    fclose(f);
    if (Y0[0] > Y0[n0 - 1])
      for (i = 0; i < n0; i++)
        for (j = 0; j < n0 / 2; j++) {
          double t = X0[i * n0 + j];
          X0[i * n0 + j] = X0[i * n0 + n0 - 1 - j];
          X0[i * n0 + n0 - 1 - j] = t;
          t = Y0[i * n0 + j];
          Y0[i * n0 + j] = Y0[i * n0 + n0 - 1 - j];
          Y0[i * n0 + n0 - 1 - j] = t;
        }
  }
  for (j = 0; j < n0; j++) {
    for (i = 0; i < n0; i++) {
      a[i] = X0[i * n0 + j];
      b[i] = Y0[i * n0 + j];
    }
    double ra[MAXN], rb[MAXN];
    refine(a, n0, RI, ra);
    refine(b, n0, RI, rb);
    for (i = 0; i < n; i++) {
      Xa[i * n0 + j] = ra[i];
      Ya[i * n0 + j] = rb[i];
    }
  }
  for (i = 0; i < n; i++) {
    double ra[MAXN], rb[MAXN];
    refine(Xa + i * n0, n0, RI, ra);
    refine(Ya + i * n0, n0, RI, rb);
    for (j = 0; j < n; j++) {
      X[i * n + j] = ra[j];
      Y[i * n + j] = rb[j];
    }
  }
  for (k = 0; k < n; k++) {
    X[k] = 0;
    X[(n - 1) * n + k] = 1;
    Y[k * n] = 0;
    Y[k * n + n - 1] = 1;
  }
  int h[9] = {n, n, 4, 0, 0, 0, 0, 0, 0};
  start(d, h, " 3000    1   10    1\n 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000\n  400.0 0.0000\n");
  for (k = 0; k < 4; k++)
    for (i = 0; i < n; i++) {
      fputs("\n\n", out);
      for (j = 0; j < n; j++)
        row(X[i * n + j], Y[i * n + j], (k - 0.5) * dz, j == n - 1 && i > 0 && i < n - 1, 0, 0, 0, 0);
    }
  finish();
  free(X0);
  free(Y0);
  free(Xa);
  free(Ya);
  free(X);
  free(Y);
}
static void lstep(const char *d) {
  int L = 44 * RI + 1, M = 44 * RJ + 1, i, j, k;
  double b = 0.5, LX = 30, dz = 2.0 / 44;
  int h[9] = {L, M, 4, 0, 0, 0, 0, 0, 0};
  start(d, h, " 3000    1   10    1\n 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000\n  300.0 0.0000\n");
  for (k = 0; k < 4; k++)
    for (i = 0; i < L; i++) {
      fputs("\n\n", out);
      double x = LX * (exp(b * i / (L - 1.0)) - 1) / (exp(b) - 1);
      for (j = 0; j < M; j++) {
        double y = 1 + (-1 + 2.0 * j / (M - 1)), e = y - 1;
        row(x, y, (k - 0.5) * dz, (e > 0 && e < 1) ? 6 * e * (1 - e) : 0, 0, 0, 0, 0);
      }
    }
  finish();
}

static void seg(double y0, double y1, int n, double qa, double qb, double *o) {
  double w[MAXN], c = 0;
  int i;
  for (i = 0; i < n; i++) {
    double t = (i + 0.5) / n;
    c += w[i] = 1 / (1 + qa * exp(-8 * t) + qb * exp(-8 * (1 - t)));
  }
  for (i = 0, o[0] = 0; i < n; i++)
    o[i + 1] = o[i] + w[i];
  for (i = 0; i <= n; i++)
    o[i] = y0 + (y1 - y0) * o[i] / c;
}

static void tstep(const char *d, int keep, int fit) {
  int M0 = fit ? 42 : 43, L0 = 45, NU = fit ? 3 : 9, Js0 = 15, i, j, k, n = 0;
  double b = fit ? 2.698 : 2.5, LD = fit ? 15.369 : 30, tu = fit ? 0.0297 : 0.04, ls = fit ? 3.741 : 0.14, bu = 2.743;
  double x0[45], y0[MAXN], x[MAXN], y[MAXN];
  for (i = 0; i < NU; i++) {
    double su = i / (NU - 1.0);
    x0[n++] = fit ? -4 * (exp(bu * (1 - su)) - 1) / (exp(bu) - 1) : -4 + 4.0 * su;
  }
  for (i = 1; i <= L0 - NU; i++)
    x0[n++] = LD * (exp(b * i / (double)(L0 - NU)) - 1) / (exp(b) - 1);
  if (fit) {
    seg(0, 1, Js0 - 1, 8.356, 15.903, y0);
    seg(1, 3, M0 - Js0, 15.903, 8.356, y0 + Js0 - 1);
  } else {
    for (j = 0; j < Js0; j++)
      y0[j] = (double)j / (Js0 - 1);
    for (j = 1; j <= M0 - Js0; j++)
      y0[Js0 - 1 + j] = 1 + 2.0 * j / (M0 - Js0);
  }
  int L, M, Is, Js;
  if (keep) {
    L = refine(x0, NU, RI, x);
    Is = L;
    L += refine_keep(x0 + NU - 1, L0 - NU + 1, RI, 1, 0, x + L - 1) - 1;
    M = refine_keep(y0, Js0, RJ, 1, 1, y);
    Js = M;
    M += refine_keep(y0 + Js0 - 1, M0 - Js0 + 1, RJ, 1, 1, y + M - 1) - 1;
  } else {
    L = refine(x0, L0, RI, x);
    M = refine(y0, M0, RJ, y);
    Is = (NU - 1) * RI + 1;
    Js = (Js0 - 1) * RJ + 1;
  }
  double k0 = 1.5 * tu * tu, e0 = pow(0.09, 0.75) * pow(k0, 1.5) / ls, dz = 3.0 / (M0 - 1);
  int h[9] = {L, M, 4, 1, Is, 1, Js, 1, 4};
  start(d, h, " 3000    2   10    1\n 0.5000 0.5000 0.5000 0.3000 0.5000 0.5000 0.7000 0.5000\n45000.0 0.0000\n");
  for (k = 0; k < 4; k++)
    for (i = 0; i < L; i++) {
      fputs("\n\n", out);
      for (j = 0; j < M; j++) {
        int wall = j == 0 || j == M - 1 || (i < Is && j < Js);
        row(x[i], y[j], (k - 0.5) * dz, wall ? 0 : x[i] > 0 ? 2.0 / 3 : fit == 2 ? 8.0 / 7 * pow(fmin(y[j] - 1, 3 - y[j]), 1.0 / 7) : 1, 0, 0, k0, e0);
      }
    }
  finish();
}

static int factors(const char *s) {
  int n = sscanf(s, "%d,%d,%d", &RI, &RJ, &RK);
  if (n < 1 || RI < 1)
    die("bad factor");
  if (n < 2)
    RJ = RI;
  if (n < 3)
    RK = RJ;
  if (RJ < 1 || RK < 1)
    die("bad factor");
  return RI;
}

int main(int argc, char **argv) {
  int r;
  if ((argc == 5 || argc == 6) && !strcmp(argv[1], "tbendg") && (r = atoi(argv[2])) >= 1) {
    int rs = argc == 6 ? atoi(argv[5]) : r;
    if (rs < 1)
      die("bad streamwise factor");
    tbendg(argv[3], r, rs, 0, argv[4]);
    return 0;
  }
  if (argc == 6 && !strcmp(argv[1], "tbendd") && (r = atoi(argv[2])) >= 1) {
    int nb = atoi(argv[5]);
    if (nb < 1)
      die("bad bend factor");
    tbendg(argv[3], r, 1, nb, argv[4]);
    return 0;
  }
  if (argc == 5 && !strcmp(argv[1], "cavityb") && (r = factors(argv[2])) >= 1) {
    cavity(argv[3], argv[4]);
    return 0;
  }
  if (argc != 4 || (r = factors(argv[2])) < 1)
    die("usage: gen cavity|bend|tbend|tbendw|lstep|tstep|tstepw|tstepf|tstepfw|tstepfp|tduct factor[,factor[,factor]] dir | gen tbendg factor dir inlet [streamwise-factor] | gen tbendd factor dir inlet bend-factor | gen cavityb factor dir nodes");
  if (!strcmp(argv[1], "bend"))
    bend(argv[3], 0, 0);
  else if (!strcmp(argv[1], "tbend"))
    bend(argv[3], 1, 0);
  else if (!strcmp(argv[1], "tbendw"))
    bend(argv[3], 1, 1);
  else if (!strcmp(argv[1], "cavity"))
    cavity(argv[3], NULL);
  else if (!strcmp(argv[1], "lstep"))
    lstep(argv[3]);
  else if (!strcmp(argv[1], "tstep"))
    tstep(argv[3], 0, 0);
  else if (!strcmp(argv[1], "tstepw"))
    tstep(argv[3], 1, 0);
  else if (!strcmp(argv[1], "tstepf"))
    tstep(argv[3], 0, 1);
  else if (!strcmp(argv[1], "tstepfw"))
    tstep(argv[3], 1, 1);
  else if (!strcmp(argv[1], "tstepfp"))
    tstep(argv[3], 0, 2);
  else if (!strcmp(argv[1], "tduct"))
    tduct(argv[3], r);
  else
    die("unknown case");
  return 0;
}
