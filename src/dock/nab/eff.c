/* 
 * eff.c: implement energy subroutines for 3 or 4 cartesian coordinates.
 *
 * Parallelization via OpenMP, creation of pair lists using a kd tree, and
 * optional calculation in 4D were added by Russ Brown (russ.brown@sun.com)
 */

#if defined(MPI) || defined(SCALAPACK)
#include "mpi.h"
#endif

/* Here is the structure for a kd tree node. */

typedef struct kdnode {
  INT_T n;
  struct kdnode *lo, *hi;
} KDNODE_T;


/***********************************************************************
                            ECONS()
************************************************************************/
 
/* Calculate the constrained energy and first derivatives. */

static
REAL_T econs( REAL_T *x, REAL_T *f )
{
  int i, foff, threadnum, numthreads;
  REAL_T e_cons, rx, ry, rz, rw;

  e_cons = 0.0;

#pragma omp parallel reduction (+: e_cons) private (i, rx, ry, rz, rw)
  {
    /*
     * Get the thread number and the number of threads for multi-threaded
     * execution under OpenMP.  For all other cases, including ScaLAPACK,
     * MPI and single-threaded execution, use the values that have been
     * stored in mytaskid and numtasks, respectively.
     */

#if defined(OPENMP)
    threadnum = omp_get_thread_num();
    numthreads = omp_get_num_threads();
#else
    threadnum = mytaskid;
    numthreads = numtasks;
#endif

    /*
     * Compute an offset into the gradient array for this thread,
     * but only if OPENMP is defined.
     */

#ifdef OPENMP
    foff = dim * prm->Natom * threadnum;
#else
    foff = 0;
#endif

    /*
     * Loop over all atoms.  Map loop indices onto OpenMP threads
     * or MPI tasks using (static, 1) scheduling.
     */

    for (i = threadnum; i < prm->Natom; i += numthreads) {
      if (constrained[i]) {
        rx = x[dim * i    ] - x0[dim * i    ];
        ry = x[dim * i + 1] - x0[dim * i + 1];
        rz = x[dim * i + 2] - x0[dim * i + 2];

        e_cons += wcons * (rx * rx + ry * ry + rz * rz);

        f[foff + dim * i    ] += 2. * wcons * rx;
        f[foff + dim * i + 1] += 2. * wcons * ry;
        f[foff + dim * i + 2] += 2. * wcons * rz;

        if (dim == 4) {
          rw = x[dim * i + 3] - x0[dim * i + 3];
          e_cons += wcons * rw * rw;
          f[foff + dim * i + 3] += 2. * wcons * rw;
        }
      }
    }
  }

  return (e_cons);
}

/***********************************************************************
                            EBOND()
************************************************************************/
 
/* Calculate the bond stretching energy and first derivatives.*/

static
REAL_T        ebond( int nbond, int *a1, int *a2, int *atype,
               REAL_T *Rk, REAL_T *Req, REAL_T *x, REAL_T *f )
{
  int i, at1, at2, atyp, foff, threadnum, numthreads;
  REAL_T e_bond, r, rx, ry, rz, rw, r2, s, db, df, e;

  e_bond = 0.0;

#pragma omp parallel reduction (+: e_bond) \
  private (i, foff, at1, at2, atyp, threadnum, numthreads, \
           rx, ry, rz, rw, r2, s, r, db, df, e )
  {
    /*
     * Get the thread number and the number of threads for multi-threaded
     * execution under OpenMP.  For all other cases, including ScaLAPACK,
     * MPI and single-threaded execution, use the values that have been
     * stored in mytaskid and numtasks, respectively.
     */

#if defined(OPENMP)
    threadnum = omp_get_thread_num();
    numthreads = omp_get_num_threads();
#else
    threadnum = mytaskid;
    numthreads = numtasks;
#endif

    /*
     * Compute an offset into the gradient array for this thread,
     * but only if OPENMP is defined.
     */

#ifdef OPENMP
    foff = dim * prm->Natom * threadnum;
#else
    foff = 0;
#endif

    /*
     * Loop over all 1-2 bonds.  Map loop indices onto OpenMP threads
     * or MPI tasks using (static, 1) scheduling.
     */

    for (i = threadnum; i < nbond; i += numthreads) {

      at1 = dim * a1[i] / 3;
      at2 = dim * a2[i] / 3;
      atyp = atype[i] - 1;

      rx = x[at1    ] - x[at2    ];
      ry = x[at1 + 1] - x[at2 + 1];
      rz = x[at1 + 2] - x[at2 + 2];
      r2 = rx * rx + ry * ry + rz * rz;

      if (dim == 4) {
        rw = x[at1 + 3] - x[at2 + 3];
        r2 += rw * rw;
      }

      s = sqrt(r2);
      r = 2.0 / s;
      db = s - Req[atyp];
      df = Rk[atyp] * db;
      e = df * db;
      e_bond += e;
      df *= r;

      f[foff + at1 + 0] += rx * df;
      f[foff + at1 + 1] += ry * df;
      f[foff + at1 + 2] += rz * df;

      f[foff + at2 + 0] -= rx * df;
      f[foff + at2 + 1] -= ry * df;
      f[foff + at2 + 2] -= rz * df;

      if (dim == 4) {
        f[foff + at1 + 3] += rw * df;
        f[foff + at2 + 3] -= rw * df;
      }
    }
  }

  return (e_bond);
}

/***********************************************************************
                            EANGL()
************************************************************************/
 
/* Calculate the bond bending energy and first derivatives. */

static
REAL_T        eangl( int nang, int *a1, int *a2, int *a3, int *atype,
               REAL_T *Tk, REAL_T *Teq, REAL_T *x, REAL_T *f )
{
  int i, atyp, at1, at2, at3, foff, threadnum, numthreads;
  REAL_T dxi, dyi, dzi, dwi, dxj, dyj, dzj, dwj, ri2, rj2, ri, rj, rir,
    rjr;
  REAL_T dxir, dyir, dzir, dwir, dxjr, dyjr, dzjr, dwjr, cst, at, da, df,
    e, e_theta;
  REAL_T xtmp, dxtmp, ytmp, wtmp, dytmp, ztmp, dztmp, dwtmp;

  e_theta = 0.0;

#pragma omp parallel reduction (+: e_theta) \
  private (i, foff, at1, at2, at3, atyp, threadnum, numthreads, \
           dxi, dyi, dzi, dwi, dxj, dyj, dzj, dwj, ri2, rj2, ri, rj, rir, rjr, \
           dxir, dyir, dzir, dwir, dxjr, dyjr, dzjr, dwjr, cst, at, da, df, e, \
           xtmp, dxtmp, ytmp, dytmp, ztmp, dztmp, wtmp, dwtmp)
  {
    /*
     * Get the thread number and the number of threads for multi-threaded
     * execution under OpenMP.  For all other cases, including ScaLAPACK,
     * MPI and single-threaded execution, use the values that have been
     * stored in mytaskid and numtasks, respectively.
     */

#if defined(OPENMP)
    threadnum = omp_get_thread_num();
    numthreads = omp_get_num_threads();
#else
    threadnum = mytaskid;
    numthreads = numtasks;
#endif

    /*
     * Compute an offset into the gradient array for this thread,
     * but only if OPENMP is defined.
     */

#ifdef OPENMP
    foff = dim * prm->Natom * threadnum;
#else
    foff = 0;
#endif

    /*
     * Loop over all 1-3 bonds.  Map loop indices onto OpenMP threads
     * or MPI tasks using (static, 1) scheduling.
     */

    for (i = threadnum; i < nang; i += numthreads) {

      at1 = dim * a1[i] / 3;
      at2 = dim * a2[i] / 3;
      at3 = dim * a3[i] / 3;
      atyp = atype[i] - 1;

      dxi = x[at1    ] - x[at2    ];
      dyi = x[at1 + 1] - x[at2 + 1];
      dzi = x[at1 + 2] - x[at2 + 2];

      dxj = x[at3    ] - x[at2    ];
      dyj = x[at3 + 1] - x[at2 + 1];
      dzj = x[at3 + 2] - x[at2 + 2];

      ri2 = dxi * dxi + dyi * dyi + dzi * dzi;
      rj2 = dxj * dxj + dyj * dyj + dzj * dzj;

      if (dim == 4) {
        dwi = x[at1 + 3] - x[at2 + 3];
        dwj = x[at3 + 3] - x[at2 + 3];
        ri2 += dwi * dwi;
        rj2 += dwj * dwj;
      }

      ri = sqrt(ri2);
      rj = sqrt(rj2);
      rir = 1. / ri;
      rjr = 1. / rj;

      dxir = dxi * rir;
      dyir = dyi * rir;
      dzir = dzi * rir;

      dxjr = dxj * rjr;
      dyjr = dyj * rjr;
      dzjr = dzj * rjr;

      cst = dxir * dxjr + dyir * dyjr + dzir * dzjr;

      if (dim == 4) {
        dwir = dwi * rir;
        dwjr = dwj * rjr;
        cst += dwir * dwjr;
      }

      if (cst > 1.0)
        cst = 1.0;
      if (cst < -1.0)
        cst = -1.0;

      at = acos(cst);
      da = at - Teq[atyp];
      df = da * Tk[atyp];
      e = df * da;
      e_theta = e_theta + e;
      df = df + df;
      at = sin(at);
      if (at > 0 && at < 1.e-3)
        at = 1.e-3;
      else if (at < 0 && at > -1.e-3)
        at = -1.e-3;
      df = -df / at;

      xtmp = df * rir * (dxjr - cst * dxir);
      dxtmp = df * rjr * (dxir - cst * dxjr);

      ytmp = df * rir * (dyjr - cst * dyir);
      dytmp = df * rjr * (dyir - cst * dyjr);

      ztmp = df * rir * (dzjr - cst * dzir);
      dztmp = df * rjr * (dzir - cst * dzjr);

      f[foff + at1 + 0] += xtmp;
      f[foff + at3 + 0] += dxtmp;
      f[foff + at2 + 0] -= xtmp + dxtmp;

      f[foff + at1 + 1] += ytmp;
      f[foff + at3 + 1] += dytmp;
      f[foff + at2 + 1] -= ytmp + dytmp;

      f[foff + at1 + 2] += ztmp;
      f[foff + at3 + 2] += dztmp;
      f[foff + at2 + 2] -= ztmp + dztmp;

      if (dim == 4) {
        wtmp = df * rir * (dwjr - cst * dwir);
        dwtmp = df * rjr * (dwir - cst * dwjr);
        f[foff + at1 + 3] += wtmp;
        f[foff + at3 + 3] += dwtmp;
        f[foff + at2 + 3] -= wtmp + dwtmp;
      }
    }
  }

  return (e_theta);
}

/***********************************************************************
                            EPHI()
************************************************************************/
 
/* Calculate the dihedral torsion energy and first derivatives. */

static
REAL_T        ephi( int nphi, int *a1, int *a2, int *a3, int *a4, int *atype,
              REAL_T *Pk, REAL_T *Pn, REAL_T *Phase, REAL_T *x, REAL_T *f )
{
  REAL_T e, co, den, co1, uu, vv, uv, ax, bx, cx, ay, by, cy, az, bz, cz,
    aw, bw, cw;
  REAL_T a0x, b0x, c0x, a0y, b0y, c0y, a0z, b0z, c0z, a0w, b0w, c0w, a1x,
    b1x;
  REAL_T a1y, b1y, a1z, b1z, a1w, b1w, a2x, b2x, a2y, b2y, a2z, b2z, a2w,
    b2w;
  REAL_T dd1x, dd2x, dd3x, dd4x, dd1y, dd2y, dd3y, dd4y, dd1z, dd2z,
    dd3z, dd4z;
  REAL_T dd1w, dd2w, dd3w, dd4w;
  REAL_T df, aa, bb, cc, ab, bc, ac, cosq;
  REAL_T ktors, phase, e_tors;
  int i, at1, at2, at3, at4, atyp, foff, threadnum, numthreads;
  REAL_T ux, uy, uz, vx, vy, vz, delta, phi, dx1, dy1, dz1, yy, pi;
  pi = 3.1415927;

  e_tors = 0.0;

#pragma omp parallel reduction (+: e_tors) \
  private (i, at1, at2, at3, at4, atyp, ax, ay, az, aw, bx, by, bz, bw, \
           cx, cy, cz, cw, ab, bc, ac, aa, bb, cc, uu, vv, uv, den, co, co1, \
           a0x, a0y, a0z, a0w, b0x, b0y, b0z, b0w, c0x, c0y, c0z, c0w, \
           a1x, a1y, a1z, a1w, b1x, b1y, b1z, b1w, a2x, a2y, a2z, a2w, \
           b2x, b2y, b2z, b2w, dd1x, dd1y, dd1z, dd1w, dd2x, dd2y, dd2z, dd2w, \
           dd3x, dd3y, dd3z, dd3w, dd4x, dd4y, dd4z, dd4w, phi, \
           ux, uy, uz, vx, vy, vz, dx1, dy1, dz1, delta, df, e, yy, phase, \
           ktors, cosq, threadnum, numthreads, foff)
  {
    /*
     * Get the thread number and the number of threads for multi-threaded
     * execution under OpenMP.  For all other cases, including ScaLAPACK,
     * MPI and single-threaded execution, use the values that have been
     * stored in mytaskid and numtasks, respectively.
     */

#if defined(OPENMP)
    threadnum = omp_get_thread_num();
    numthreads = omp_get_num_threads();
#else
    threadnum = mytaskid;
    numthreads = numtasks;
#endif

    /*
     * Compute an offset into the gradient array for this thread,
     * but only if OPENMP is defined.
     */

#ifdef OPENMP
    foff = dim * prm->Natom * threadnum;
#else
    foff = 0;
#endif

    /*
     * Loop over all 1-4 bonds.  Map loop indices onto OpenMP threads
     * or MPI tasks using (static, 1) scheduling.
     */

    for (i = threadnum; i < nphi; i += numthreads) {

      at1 = dim * a1[i] / 3;
      at2 = dim * a2[i] / 3;
      at3 = dim * abs(a3[i]) / 3;
      at4 = dim * abs(a4[i]) / 3;
      atyp = atype[i] - 1;

      ax = x[at2 + 0] - x[at1 + 0];
      ay = x[at2 + 1] - x[at1 + 1];
      az = x[at2 + 2] - x[at1 + 2];

      bx = x[at3 + 0] - x[at2 + 0];
      by = x[at3 + 1] - x[at2 + 1];
      bz = x[at3 + 2] - x[at2 + 2];

      cx = x[at4 + 0] - x[at3 + 0];
      cy = x[at4 + 1] - x[at3 + 1];
      cz = x[at4 + 2] - x[at3 + 2];

      if (dim == 4) {
        aw = x[at2 + 3] - x[at1 + 3];
        bw = x[at3 + 3] - x[at2 + 3];
        cw = x[at4 + 3] - x[at3 + 3];

#     define DOT4(a,b,c,d,e,f,g,h) a*e + b*f + c*g + d*h

        ab = DOT4(ax, ay, az, aw, bx, by, bz, bw);
        bc = DOT4(bx, by, bz, bw, cx, cy, cz, cw);
        ac = DOT4(ax, ay, az, aw, cx, cy, cz, cw);
        aa = DOT4(ax, ay, az, aw, ax, ay, az, aw);
        bb = DOT4(bx, by, bz, bw, bx, by, bz, bw);
        cc = DOT4(cx, cy, cz, cw, cx, cy, cz, cw);
      } else {

#     define DOT3(a,b,c,d,e,f) a*d + b*e + c*f

        ab = DOT3(ax, ay, az, bx, by, bz);
        bc = DOT3(bx, by, bz, cx, cy, cz);
        ac = DOT3(ax, ay, az, cx, cy, cz);
        aa = DOT3(ax, ay, az, ax, ay, az);
        bb = DOT3(bx, by, bz, bx, by, bz);
        cc = DOT3(cx, cy, cz, cx, cy, cz);
      }

      uu = (aa * bb) - (ab * ab);
      vv = (bb * cc) - (bc * bc);
      uv = (ab * bc) - (ac * bb);
      den = 1.0 / sqrt(uu * vv);
      co = uv * den;
      co1 = 0.5 * co * den;

      a0x = -bc * bx + bb * cx;
      a0y = -bc * by + bb * cy;
      a0z = -bc * bz + bb * cz;

      b0x = ab * cx + bc * ax - 2. * ac * bx;
      b0y = ab * cy + bc * ay - 2. * ac * by;
      b0z = ab * cz + bc * az - 2. * ac * bz;

      c0x = ab * bx - bb * ax;
      c0y = ab * by - bb * ay;
      c0z = ab * bz - bb * az;

      a1x = 2. * uu * (-cc * bx + bc * cx);
      a1y = 2. * uu * (-cc * by + bc * cy);
      a1z = 2. * uu * (-cc * bz + bc * cz);

      b1x = 2. * uu * (bb * cx - bc * bx);
      b1y = 2. * uu * (bb * cy - bc * by);
      b1z = 2. * uu * (bb * cz - bc * bz);

      a2x = -2. * vv * (bb * ax - ab * bx);
      a2y = -2. * vv * (bb * ay - ab * by);
      a2z = -2. * vv * (bb * az - ab * bz);

      b2x = 2. * vv * (aa * bx - ab * ax);
      b2y = 2. * vv * (aa * by - ab * ay);
      b2z = 2. * vv * (aa * bz - ab * az);

      dd1x = (a0x - a2x * co1) * den;
      dd1y = (a0y - a2y * co1) * den;
      dd1z = (a0z - a2z * co1) * den;

      dd2x = (-a0x - b0x - (a1x - a2x - b2x) * co1) * den;
      dd2y = (-a0y - b0y - (a1y - a2y - b2y) * co1) * den;
      dd2z = (-a0z - b0z - (a1z - a2z - b2z) * co1) * den;

      dd3x = (b0x - c0x - (-a1x - b1x + b2x) * co1) * den;
      dd3y = (b0y - c0y - (-a1y - b1y + b2y) * co1) * den;
      dd3z = (b0z - c0z - (-a1z - b1z + b2z) * co1) * den;

      dd4x = (c0x - b1x * co1) * den;
      dd4y = (c0y - b1y * co1) * den;
      dd4z = (c0z - b1z * co1) * den;

      if (dim == 4) {
        a0w = -bc * bw + bb * cw;
        b0w = ab * cw + bc * aw - 2. * ac * bw;
        c0w = ab * bw - bb * aw;
        a1w = 2. * uu * (-cc * bw + bc * cw);
        b1w = 2. * uu * (bb * cw - bc * bw);
        a2w = -2. * vv * (bb * aw - ab * bw);
        b2w = 2. * vv * (aa * bw - ab * aw);
        dd1w = (a0w - a2w * co1) * den;
        dd2w = (-a0w - b0w - (a1w - a2w - b2w) * co1) * den;
        dd3w = (b0w - c0w - (-a1w - b1w + b2w) * co1) * den;
        dd4w = (c0w - b1w * co1) * den;
      }

      if (prm->Nhparm && a3[i] < 0) {

        /*   here we will use a quadratic form for the improper torsion  */
        /*     we are using the NHPARM variable in prmtop to trigger this   */
        /*
          WARNING: phi itself is here calculated from the first three coords--
          --- may fail!
        */

        /* Note: The following improper torsion code does not support 4D! */

        co = co > 1.0 ? 1.0 : co;
        co = co < -1.0 ? -1.0 : co;

        phi = acos(co);

        /*
          now calculate sin(phi) because cos(phi) is symmetric, so
          we can decide between +-phi.
        */

        ux = ay * bz - az * by;
        uy = az * bx - ax * bz;
        uz = ax * by - ay * bx;

        vx = by * cz - bz * cy;
        vy = bz * cx - bx * cz;
        vz = bx * cy - by * cx;

        dx1 = uy * vz - uz * vy;
        dy1 = uz * vx - ux * vz;
        dz1 = ux * vy - uy * vx;

        dx1 = DOT3(dx1, dy1, dz1, bx, by, bz);
        if (dx1 < 0.0)
          phi = -phi;

        delta = phi - Phase[atyp];
        delta = delta > pi ? pi : delta;
        delta = delta < -pi ? -pi : delta;

        df = Pk[atyp] * delta;
        e = df * delta;
        e_tors += e;
        yy = sin(phi);

        /*
          Decide what expansion to use
          Check first for the "normal" expression, since it will be
          the most used

          the 0.001 value could be lowered for increased precision.
          This insures ~1e-05% error for sin(phi)=0.001
        */

        if (fabs(yy) > 0.001) {
          df = -2.0 * df / yy;
        } else {
          if (fabs(delta) < 0.10) {
            if (Phase[atyp] == 0.0) {
              df = -2.0 * Pk[atyp] * (1 + phi * phi / 6.0);
            } else {
              if (fabs(Phase[atyp]) == pi) {
                df = 2.0 * Pk[atyp] * (1 +
                                       delta * delta / 6.0);
              }
            }
          } else {
            if ((phi > 0.0 && phi < (pi / 2.0)) ||
                (phi < 0.0 && phi > -pi / 2.0))
              df = df * 1000.;
            else
              df = -df * 1000.;
          }
        }
      } else {
      multi_term:

        if (fabs(Phase[atyp] - 3.142) < 0.01)
          phase = -1.0;
        else
          phase = 1.0;

        ktors = Pk[atyp];
        switch ((int) fabs(Pn[atyp])) {

        case 1:
          e = ktors * (1.0 + phase * co);
          df = phase * ktors;
          break;
        case 2:
          e = ktors * (1.0 + phase * (2.*co*co -1.));
          df = phase*ktors*4.*co;
          break;
        case 3:
          cosq = co*co;
          e = ktors * (1.0 + phase * co*(4.*cosq -3.));
          df = phase*ktors*(12.*cosq - 3.);
          break;
        case 4:
          cosq = co*co;
          e = ktors * (1.0 + phase * (8.*cosq*(cosq - 1.) + 1.));
          df = phase*ktors*co*( 32.*cosq -16.);
          break;
        case 6:
          cosq = co*co;
          e = ktors * (1.0 + phase * (32.*cosq*cosq*cosq -
                                      48.*cosq*cosq + 18.*cosq - 1.));
          df = phase*ktors*co*( 192.*cosq*cosq - 192.*cosq + 36. );
          break;
        default:
          fprintf( stderr, "bad value for Pn: %d %d %d %d %8.3f\n",
                  at1, at2, at3, at4, Pn[atyp]);
          exit(1);
        }
        e_tors += e;

      }

      f[foff + at1 + 0] += df * dd1x;
      f[foff + at1 + 1] += df * dd1y;
      f[foff + at1 + 2] += df * dd1z;

      f[foff + at2 + 0] += df * dd2x;
      f[foff + at2 + 1] += df * dd2y;
      f[foff + at2 + 2] += df * dd2z;

      f[foff + at3 + 0] += df * dd3x;
      f[foff + at3 + 1] += df * dd3y;
      f[foff + at3 + 2] += df * dd3z;

      f[foff + at4 + 0] += df * dd4x;
      f[foff + at4 + 1] += df * dd4y;
      f[foff + at4 + 2] += df * dd4z;

      if (dim == 4) {
        f[foff + at1 + 3] += df * dd1w;
        f[foff + at2 + 3] += df * dd2w;
        f[foff + at3 + 3] += df * dd3w;
        f[foff + at4 + 3] += df * dd4w;
      }

#ifdef PRINT_EPHI
      fprintf( nabout,"%4d%4d%4d%4d%4d%8.3f\n", i + 1, at1, at2, at3, at4, e);
      fprintf( nabout,"%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
             -df * dd1x, -df * dd1y, -df * dd1z, -df * dd2x, -df * dd2y,
             -df * dd2z, -df * dd3x, -df * dd3y);
      fprintf( nabout,"%10.5f%10.5f%10.5f%10.5f\n", -df * dd3z, -df * dd4x,
             -df * dd4y, -df * dd4z);
#endif
      if (Pn[atyp] < 0.0) {
        atyp++;
        goto multi_term;
      }
    }
  }

  return (e_tors);
}

#if 0
/***********************************************************************
                            EIMPHI()
************************************************************************/
 
/* 
 * This function for calculation of improper torsion is included because
 * it has historically been included in eff.c.  However, it has not been
 * modified to support 4D, and it is not called by any other function.
 */

static
REAL_T eimphi( int nphi, int *a1, int *a2, int *a3, int *a4, int *atype,
               REAL_T *Pk, REAL_T *Peq, REAL_T *x, REAL_T *f )
{

  /*
    calculate improper dihedral energies and forces according
    to T.Schlick, J.Comput.Chem, vol 10, 7 (1989) 
    with local variant (no delta funcion involved)

    If |sin(phi)|.gt.1e-03 use the angle directly
    V = kimp ( phi -impeq ) ^2
    dV/dcos(phi) = -2*k/sin(phi)

    If |sin(phi)|.lt.1e-03 and (phi-Peq).lt.0.1
    and cnseq(phi)=0.
    V = kimp (phi - impeq)^2
    dV/dcos(phi) = -2*k*(1+phi*phi/6.)

    If |sin(phi)|.lt.1e-03 and (phi-Peq).lt.0.1
    and cnseq(phi)=pi
    V = kimp (phi - impeq)^2
    dV/dcos(phi) = 2*k*(1+phi*phi/6.)

    Now, if |sin(phi)|.lt.1e-03 but far from minimum,
    set sin(phi) to 1e-03 arbitrarily

  */

  REAL_T e;
  REAL_T co, phi, den, co1;
  REAL_T uu, vv, uv;
  REAL_T ax, bx, cx;
  REAL_T ay, by, cy;
  REAL_T az, bz, cz;
  REAL_T ux, uy, uz;
  REAL_T vx, vy, vz;
  REAL_T dx1, dy1, dz1;
  REAL_T a0x, b0x, c0x;
  REAL_T a0y, b0y, c0y;
  REAL_T a0z, b0z, c0z;
  REAL_T a1x, b1x;
  REAL_T a1y, b1y;
  REAL_T a1z, b1z;
  REAL_T a2x, b2x;
  REAL_T a2y, b2y;
  REAL_T a2z, b2z;
  REAL_T dd1x, dd2x, dd3x, dd4x;
  REAL_T dd1y, dd2y, dd3y, dd4y;
  REAL_T dd1z, dd2z, dd3z, dd4z;
  REAL_T df, delta, pi, yy;
  REAL_T aa, bb, cc, ab, bc, ac;
  REAL_T e_imp;
  int i;
  int at1, at2, at3, at4, atyp;

  pi = 3.1415927;

  e_imp = 0.0;
  for (i = 0; i < nphi; i++) {

    at1 = a1[i];
    at2 = a2[i];
    at3 = abs(a3[i]);
    at4 = abs(a4[i]);
    atyp = atype[i] - 1;


    /*        calculate distances  */

    ax = x[at2 + 0] - x[at1 + 0];
    ay = x[at2 + 1] - x[at1 + 1];
    az = x[at2 + 2] - x[at1 + 2];
    bx = x[at3 + 0] - x[at2 + 0];
    by = x[at3 + 1] - x[at2 + 1];
    bz = x[at3 + 2] - x[at2 + 2];
    cx = x[at4 + 0] - x[at3 + 0];
    cy = x[at4 + 1] - x[at3 + 1];
    cz = x[at4 + 2] - x[at3 + 2];

#                define DOT(a,b,c,d,e,f) a*d + b*e + c*f
    ab = DOT(ax, ay, az, bx, by, bz);
    bc = DOT(bx, by, bz, cx, cy, cz);
    ac = DOT(ax, ay, az, cx, cy, cz);
    aa = DOT(ax, ay, az, ax, ay, az);
    bb = DOT(bx, by, bz, bx, by, bz);
    cc = DOT(cx, cy, cz, cx, cy, cz);


    uu = (aa * bb) - (ab * ab);
    vv = (bb * cc) - (bc * bc);
    uv = (ab * bc) - (ac * bb);

    den = 1.0 / sqrt(fabs(uu * vv));
    co = uv * den;
    co = co > 1.0 ? 1.0 : co;
    co = co < -1.0 ? -1.0 : co;

    phi = acos(co);

    /*
      now calculate sin(phi) because cos(phi) is symmetric, so
      we can decide between +-phi.
    */

    ux = ay * bz - az * by;
    uy = az * bx - ax * bz;
    uz = ax * by - ay * bx;

    vx = by * cz - bz * cy;
    vy = bz * cx - bx * cz;
    vz = bx * cy - by * cx;

    dx1 = uy * vz - uz * vy;
    dy1 = uz * vx - ux * vz;
    dz1 = ux * vy - uy * vx;

    dx1 = DOT(dx1, dy1, dz1, bx, by, bz);
    if (dx1 < 0.0)
      phi = -phi;

    delta = phi - Peq[atyp];
    delta = delta > pi ? pi : delta;
    delta = delta < -pi ? -pi : delta;

    df = Pk[atyp] * delta;
    e = df * delta;
    e_imp += e;
    yy = sin(phi);

    /*
      Decide what expansion to use
      Check first for the "normal" expression, since it will be
      the most used

      the 0.001 value could be lowered for increased precision.
      This insures ~1e-05% error for sin(phi)=0.001
    */

    if (fabs(yy) > 0.001) {
      df = -2.0 * df / yy;
    } else {
      if (fabs(delta) < 0.10) {
        if (Peq[atyp] == 0.0) {
          df = -2.0 * Pk[atyp] * (1 + phi * phi / 6.0);
        } else {
          if (fabs(Peq[atyp]) == pi)
            df = 2.0 * Pk[atyp] * (1 + delta * delta / 6.0);
        }
      } else {
        if ((phi > 0.0 && phi < (pi / 2.0)) ||
            (phi < 0.0 && phi > -pi / 2.0))
          df = df * 1000.;
        else
          df = -df * 1000.;
      }
    }

    co1 = 0.5 * co * den;

    a0x = -bc * bx + bb * cx;
    a0y = -bc * by + bb * cy;
    a0z = -bc * bz + bb * cz;

    b0x = ab * cx + bc * ax - 2. * ac * bx;
    b0y = ab * cy + bc * ay - 2. * ac * by;
    b0z = ab * cz + bc * az - 2. * ac * bz;

    c0x = ab * bx - bb * ax;
    c0y = ab * by - bb * ay;
    c0z = ab * bz - bb * az;

    a1x = 2. * uu * (-cc * bx + bc * cx);
    a1y = 2. * uu * (-cc * by + bc * cy);
    a1z = 2. * uu * (-cc * bz + bc * cz);

    b1x = 2. * uu * (bb * cx - bc * bx);
    b1y = 2. * uu * (bb * cy - bc * by);
    b1z = 2. * uu * (bb * cz - bc * bz);

    a2x = -2. * vv * (bb * ax - ab * bx);
    a2y = -2. * vv * (bb * ay - ab * by);
    a2z = -2. * vv * (bb * az - ab * bz);

    b2x = 2. * vv * (aa * bx - ab * ax);
    b2y = 2. * vv * (aa * by - ab * ay);
    b2z = 2. * vv * (aa * bz - ab * az);

    dd1x = (a0x - a2x * co1) * den;
    dd1y = (a0y - a2y * co1) * den;
    dd1z = (a0z - a2z * co1) * den;

    dd2x = (-a0x - b0x - (a1x - a2x - b2x) * co1) * den;
    dd2y = (-a0y - b0y - (a1y - a2y - b2y) * co1) * den;
    dd2z = (-a0z - b0z - (a1z - a2z - b2z) * co1) * den;

    dd3x = (b0x - c0x - (-a1x - b1x + b2x) * co1) * den;
    dd3y = (b0y - c0y - (-a1y - b1y + b2y) * co1) * den;
    dd3z = (b0z - c0z - (-a1z - b1z + b2z) * co1) * den;

    dd4x = (c0x - b1x * co1) * den;
    dd4y = (c0y - b1y * co1) * den;
    dd4z = (c0z - b1z * co1) * den;

    f[at1 + 0] += df * dd1x;
    f[at1 + 1] += df * dd1y;
    f[at1 + 2] += df * dd1z;

    f[at2 + 0] += df * dd2x;
    f[at2 + 1] += df * dd2y;
    f[at2 + 2] += df * dd2z;

    f[at3 + 0] += df * dd3x;
    f[at3 + 1] += df * dd3y;
    f[at3 + 2] += df * dd3z;

    f[at4 + 0] += df * dd4x;
    f[at4 + 1] += df * dd4y;
    f[at4 + 2] += df * dd4z;

  }

  return (e_imp);
}
#endif

/***********************************************************************
                            DOWNHEAP_PAIRS()
************************************************************************/
 
/*
 * The downheap function from Robert Sedgewick's "Algorithms in C++" p. 152,
 * corrected for the fact that Sedgewick indexes the heap array from 1 to n
 * whereas Java indexes the heap array from 0 to n-1. Note, however, that
 * the heap should be indexed conceptually from 1 to n in order that for
 * any node k the two children are found at nodes 2*k and 2*k+1. Move down
 * the heap, exchanging the node at position k with the larger of its two
 * children if necessary and stopping when the node at k is larger than both
 * children or the bottom of the heap is reached. Note that it is possible
 * for the node at k to have only one child: this case is treated properly.
 * A full exchange is not necessary because the variable 'v' is involved in
 * the exchanges.  The 'while' loop has two exits: one for the case that
 * the bottom of the heap is hit, and another for the case that the heap
 * condition (the parent is greater than or equal to both children) is
 * satisfied somewhere in the interior of the heap.
 *
 * Used by the heapsort_pairs function which sorts the pair list arrays.
 *
 * Calling parameters are as follows:
 *
 * a - array of indices into the atomic coordinate array x
 * n - the number of items to be sorted
 * k - the exchange node (or element)
 */

static
void downheap_pairs(int *a, int n, int k) {

  int j, v;

  v = a[k-1];
  while (k <= n/2) {
    j = k + k;
    if ( (j < n) && (a[j-1] < a[j]) ) j++;
    if (v >= a[j-1]) break;
    a[k-1] = a[j-1];
    k = j;
  }
  a[k-1] = v;
}

/***********************************************************************
                            HEAPSORT_PAIRS()
************************************************************************/
 
/*
 * The heapsort function from Robert Sedgewick's "Algorithms in C++" p. 156,
 * corrected for the fact that Sedgewick indexes the heap array from 1 to n
 * whereas Java indexes the heap array from 0 to n-1. Note, however, that
 * the heap should be indexed conceptually from 1 to n in order that for
 * any node k the two children are found at nodes 2*k and 2*k+1.  In what
 * follows, the 'for' loop heaporders the array in linear time and the
 * 'while' loop exchanges the largest element with the last element then
 * repairs the heap.
 *
 * Calling parameters are as follows:
 *
 * a - array of indices into the atomic coordinate array x
 * n - the number of items to be sorted
 *
 * Used by the nblist function to sort the pair list arrays.
 */

static
void heapsort_pairs(int *a, int n) {

  int k, v;

  for (k = n/2; k >= 1; k--) downheap_pairs(a, n, k);
  while (n > 1) {
    v = a[0];
    a[0] = a[n-1];
    a[n-1] = v;
    downheap_pairs(a, --n, 1);
  }
}

/***********************************************************************
                            DOWNHEAP_INDEX()
************************************************************************/
 
/*
 * The downheap function from Robert Sedgewick's "Algorithms in C++" p. 152,
 * corrected for the fact that Sedgewick indexes the heap array from 1 to n
 * whereas Java indexes the heap array from 0 to n-1. Note, however, that
 * the heap should be indexed conceptually from 1 to n in order that for
 * any node k the two children are found at nodes 2*k and 2*k+1. Move down
 * the heap, exchanging the node at position k with the larger of its two
 * children if necessary and stopping when the node at k is larger than both
 * children or the bottom of the heap is reached. Note that it is possible
 * for the node at k to have only one child: this case is treated properly.
 * A full exchange is not necessary because the variable 'v' is involved in
 * the exchanges.  The 'while' loop has two exits: one for the case that
 * the bottom of the heap is hit, and another for the case that the heap
 * condition (the parent is greater than or equal to both children) is
 * satisfied somewhere in the interior of the heap.
 *
 * Used by the heapsort_index function which sorts the index arrays indirectly
 * by comparing components of the Cartesian coordinates.
 *
 * Calling parameters are as follows:
 *
 * a - array of indices into the atomic coordinate array x
 * n - the number of items to be sorted
 * k - the exchange node (or element)
 * x - the atomic coordinate array
 * p - the partition (x, y, z or w) on which sorting occurs
 */

static
void downheap_index(int *a, int n, int k, REAL_T *x, int p) {

  int j, v;

  v = a[k-1];
  while (k <= n/2) {
    j = k + k;
    if ( (j < n) && (x[dim * a[j-1] + p] < x[dim * a[j] + p]) ) j++;
    if (x[dim * v + p] >= x[dim * a[j-1] + p]) break;
    a[k-1] = a[j-1];
    k = j;
  }
  a[k-1] = v;
}

/***********************************************************************
                            HEAPSORT_INDEX()
************************************************************************/
 
/*
 * The heapsort function from Robert Sedgewick's "Algorithms in C++" p. 156,
 * corrected for the fact that Sedgewick indexes the heap array from 1 to n
 * whereas Java indexes the heap array from 0 to n-1. Note, however, that
 * the heap should be indexed conceptually from 1 to n in order that for
 * any node k the two children are found at nodes 2*k and 2*k+1.  In what
 * follows, the 'for' loop heaporders the array in linear time and the
 * 'while' loop exchanges the largest element with the last element then
 * repairs the heap.
 *
 * Calling parameters are as follows:
 *
 * a - array of indices into the atomic coordinate array x
 * n - the number of items to be sorted
 * x - the atomic coordinate array
 * p - the partition (x, y, z or w) on which sorting occurs
 *
 * Used by the nblist function to sort the xn, yn, zn, wn and on arrays.
 */

static
void heapsort_index(int *a, int n, REAL_T *x, int p) {

  int k, v;

  for (k = n/2; k >= 1; k--) downheap_index(a, n, k, x, p);
  while (n > 1) {
    v = a[0];
    a[0] = a[n-1];
    a[n-1] = v;
    downheap_index(a, --n, 1, x, p);
  }
}

/***********************************************************************
                            BUILDKDTREE()
************************************************************************/
 
/*
 * Build the kd tree by recursively subdividing the atom number
 * arrays and adding nodes to the tree.  Note that the arrays
 * are permuted cyclically as control moves down the tree in
 * order that sorting occur on x, y, z and (for 4D) w.  Also,
 * if it is desired that the kd tree provide a partial atom
 * order, the sorting will occur on o, x, y, z and (for 4D) w.  The
 * temporary array is provided for the copy and partition operation.
 *
 * Calling parameters are as follows:
 *
 * xn - x sorted array of atom numbers
 * yn - y sorted array of atom numbers
 * zn - z sorted array of atom numbers
 * wn - w sorted array of atom numbers
 * on - ordinal array of atom numbers
 * tn - temporary array for atom numbers
 * start - first element of array 
 * end - last element of array
 * kdpptr - pointer to pointer to kd tree node next available for allocation
 * that - the node currently visited, the equivalent of 'this' in C++
 * x - atomic coordinate array
 * p - the partition (x, y, z, w or o) on which sorting occurs
 */

#define SORT_ATOM_NUMBERS

static
void buildkdtree(int *xn, int *yn, int *zn, int *wn, int *on, int *tn,
                 int start, int end, KDNODE_T **kdpptr, KDNODE_T *that,
                 REAL_T *x, int p)
{
  int i, middle, imedian, lower, upper;
  REAL_T median;

  /*
   * The partition cycles by dim unless SORT_ATOM_NUMBERS is defined,
   * in which case it cycles by dim+1.  Note that if SORT_ATOM_NUMBERS
   * is defined and the partition equals zero, sorting will occur
   * on the ordinal atom number instead of the atom's cartesian
   * coordinate.
   */

#ifndef SORT_ATOM_NUMBERS
  p %= dim;
#else
  p %= (dim + 1);
#endif

  /* If only one element is passed to this function, add it to the tree. */

  if (end == start) {
    that->n = xn[start];
  }

  /*
   * Otherwise, if two elements are passed to this function, determine
   * whether the first element is the low child or the high child.  Or,
   * if neither element is the low child, choose the second element to
   * be the high child.  Allocate a new KDNODE_T and make it one or the
   * other of the children.
   */

  else if (end == start + 1) {

    /* Check whether the first element is the low child. */

#ifdef SORT_ATOM_NUMBERS
    if ( ((p == 0) && (xn[start] < xn[end])) ||
         ((p != 0) && (x[dim * xn[start] + p-1] <
		       x[dim * xn[end] + p-1])) )
#else
    if (x[dim * xn[start] + p] < x[dim * xn[end] + p])
#endif
      {
	that->n = xn[end];
	(*kdpptr)->n = xn[start];
	(*kdpptr)->lo = NULL;
	(*kdpptr)->hi = NULL;
	that->lo = (*kdpptr)++;
      }

    /* Check whether the second element is the low child. */

#ifdef SORT_ATOM_NUMBERS
    else if ( ((p == 0) && (xn[start] > xn[end])) ||
	      ((p != 0) && (x[dim * xn[start] + p-1] >
			    x[dim * xn[end] + p-1])) )
#else
    if (x[dim * xn[start] + p] > x[dim * xn[end] + p])
#endif
      {
	that->n = xn[start];
	(*kdpptr)->n = xn[end];
	(*kdpptr)->lo = NULL;
	(*kdpptr)->hi = NULL;
	that->lo = (*kdpptr)++;
      }

    /* Neither element is the low child so use the second as the high child. */

    else
      {
	that->n = xn[start];
	(*kdpptr)->n = xn[end];
	(*kdpptr)->lo = NULL;
	(*kdpptr)->hi = NULL;
	that->hi = (*kdpptr)++;
      }
  }

  /* Otherwise, more than two elements are passed to this function. */

  else {

    /*
     * The middle element of the xn array is taken as the element about
     * which the yn and zn arrays will be partitioned.  However, search
     * lower elements of the xn array to ensure that the p values of the
     * atomic coordinate array that correspond to these elements are indeed
     * less than the median value because they may be equal to it.  This
     * approach is consistent with partitioning between < and >=.
     *
     * The search described above is not necessary if SORT_ATOM_NUMBERS is
     * defined and p==0 because in this case sorting occurs on the
     * ordinal atom number instead of the atomic coordinate, and the
     * ordinal atom numbers are all unique.
     */

    middle = (start + end)/2;
#ifdef SORT_ATOM_NUMBERS
    if (p == 0) {
      imedian = xn[middle];
    } else {
      median = x[dim * xn[middle] + p-1];
      for (i = middle - 1; i >= start; i--) {
        if (x[dim * xn[i] + p-1] < median) {
          break;
        } else {
          middle = i;
        }
      }
    }
#else
    median = x[dim * xn[middle] + p];
    for (i = middle - 1; i >= start; i--) {
      if (x[dim * xn[i] + p] < median) {
        break;
      } else {
        middle = i;
      }
    }
#endif

    /* Store the middle element at this kd node. */

    that->n = xn[middle];

    /*
     * Scan the yn array in ascending order and compare the p value of
     * each corresponding element of the atomic coordinate array to the
     * median value.  If the p value is less than the median value, copy
     * the element of the yn array into the lower part of the tn array.
     * If the p value is greater than or equal to the median value, copy
     * the element of the yn array into the upper part of the tn array.
     * The lower part of the tn array begins with the start index, and the
     * upper part of the tn array begins one element above the middle index.
     * At the end of this scan and copy operation, the tn array will have
     * been subdivided into three groups: (1) a group of indices beginning
     * with start and continuing up to but not including middle, which indices
     * point to atoms for which the p value is less than the median value;
     * (2) the middle index that has been stored in this node of  the kd tree;
     * and (3) a group of indices beginning one address above middle and
     * continuing up to and including end, which indices point to atoms for
     * which the p value is greater than or equal to the median value.
     *
     * This approach preserves the relative heapsorted order of elements
     * of the atomic coordinate array that correspond to elements of the
     * yn array while those elements are partitioned about the p median.
     *
     * Note: when scanning the yn array, skip the element (i.e., the atom
     * number) that equals the middle element because that atom number has
     * been stored at this node of the kd-tree.
     */

    lower = start - 1;
    upper = middle;
    for (i = start; i <= end; i++) {
      if (yn[i] != xn[middle]) {

#ifdef SORT_ATOM_NUMBERS
        if ( ((p == 0) && (yn[i] < imedian)) ||
             ((p != 0) && (x[dim * yn[i] + p-1] < median)) )
#else
        if (x[dim * yn[i] + p] < median)
#endif
          {
            tn[++lower] = yn[i];
          } else {
            tn[++upper] = yn[i];
          }
      }
    }

    /*
     * All elements of the yn array between start and end have been copied
     * and partitioned into the tn array, so the yn array is available for
     * elements of the zn array to be copied and partitioned into the yn
     * array, in the same manner in which elements of the yn array were
     * copied and partitioned into the tn array.
     *
     * This approach preserves the relative heapsorted order of elements
     * of the atomic coordinate array that correspond to elements of the
     * zn array while those elements are partitioned about the p median.
     *
     * Note: when scanning the zn array, skip the element (i.e., the atom
     * number) that equals the middle element because that atom number has
     * been stored at this node of the kd-tree.
     */

    lower = start - 1;
    upper = middle;
    for (i = start; i <= end; i++) {
      if (zn[i] != xn[middle]) {

#ifdef SORT_ATOM_NUMBERS
        if ( ((p == 0) && (zn[i] < imedian)) ||
             ((p != 0) && (x[dim * zn[i] + p-1] < median)) )
#else
        if (x[dim * zn[i] + p] < median)
#endif
          {
            yn[++lower] = zn[i];
          } else {
            yn[++upper] = zn[i];
          }
      }
    }

    /*
     * Execute the following region of code if SORT_ATOM_NUMBERS is defined,
     * or if SORT_ATOM_NUMBERS is not defined and dim==4.
     */

#ifndef SORT_ATOM_NUMBERS
    if (dim == 4)
#endif

      {

        /*
         * All elements of the zn array between start and end have been copied
         * and partitioned into the yn array, so the zn array is available for
         * elements of the wn array to be copied and partitioned into the zn
         * array, in the same manner in which elements of the zn array were
         * copied and partitioned into the yn array.
         *
         * This approach preserves the relative heapsorted order of elements
         * of the atomic coordinate array that correspond to elements of the
         * wn array while those elements are partitioned about the p median.
         *
         * Note: when scanning the wn array, skip the element (i.e., the atom
         * number) that equals the middle element because that atom number has
         * been stored at this node of the kd-tree.
         */

        lower = start - 1;
        upper = middle;
        for (i = start; i <= end; i++) {
          if (wn[i] != xn[middle]) {

#ifdef SORT_ATOM_NUMBERS
            if ( ((p == 0) && (wn[i] < imedian)) ||
                 ((p != 0) && (x[dim * wn[i] + p-1] < median)) )
#else
            if (x[dim * wn[i] + p] < median)
#endif
              {
                zn[++lower] = wn[i];
              } else {
                zn[++upper] = wn[i];
              }
          }
        }
      }

    /*
     * Execute the following region of code if SORT_ATOM_NUMBERS is defined
     * and dim==4.
     */

#ifdef SORT_ATOM_NUMBERS

    if (dim == 4) {

      /*
       * All elements of the wn array between start and end have been copied
       * and partitioned into the zn array, so the wn array is available for
       * elements of the on array to be copied and partitioned into the wn
       * array, in the same manner in which elements of the wn array were
       * copied and partitioned into the zn array.
       *
       * This approach preserves the relative heapsorted order of elements
       * of the atomic coordinate array that correspond to elements of the
       * wn array while those elements are partitioned about the p median.
       *
       * Note: when scanning the on array, skip the element (i.e., the atom
       * number) that equals the middle element because that atom number has
       * been stored at this node of the kd-tree.
       */

      lower = start - 1;
      upper = middle;
      for (i = start; i <= end; i++) {
        if (on[i] != xn[middle]) {
          if ( ((p == 0) && (on[i] < imedian)) ||
               ((p != 0) && (x[dim * on[i] + p-1] < median)) )
            {
              wn[++lower] = on[i];
            } else {
              wn[++upper] = on[i];
            }
        }
      }
    }

#endif

    /*
     * Recurse down the lo branch of the tree if the lower group of
     * the tn array is non-null.  Note permutation of the xn, yn, zn, wn,
     * on and tn arrays.  In particular, xn was used for partitioning at
     * this level of the tree.  At one level down the tree, yn (which
     * has been copied into tn) will be used for partitioning.  At two
     * levels down the tree, zn (which has been copied into yn) will
     * be used for partitioning.  If SORT_ATOM_NUMBERS is defined, or if
     * SORT_ATOM_NUMBERS is not defined and dim==4, at three levels down the
     * tree, wn (which has been copied into zn) will be used for partitoning.
     * At four levels down the tree, xn will be used for partitioning.
     * In this manner, partitioning cycles through xn, yn, zn and wn
     * at successive levels of the tree.
     *
     * Note that for 3D the wn array isn't allocated so don't permute it
     * cyclically along with the other arrays in the recursive call.
     *
     * Note also that if SORT_ATOM_NUMBERS isn't defined the on array isn't
     * allocated so don't permute it cyclically in the recursive call.
     */

    if (lower >= start) {
      (*kdpptr)->lo = NULL;
      (*kdpptr)->hi = NULL;
      that->lo = (*kdpptr)++;

#ifndef SORT_ATOM_NUMBERS
      if (dim == 4) {
        buildkdtree(tn, yn, zn, xn, on, wn,
                    start, lower, kdpptr, that->lo, x, p+1);
      } else {
        buildkdtree(tn, yn, xn, wn, on, zn,
                    start, lower, kdpptr, that->lo, x, p+1);
      }
#else
      if (dim == 4) {
        buildkdtree(tn, yn, zn, wn, xn, on,
                    start, lower, kdpptr, that->lo, x, p+1);
      } else {
        buildkdtree(tn, yn, zn, xn, on, wn,
                    start, lower, kdpptr, that->lo, x, p+1);
      }
#endif

    }
      
    /*
     * Recurse down the hi branch of the tree if the upper group of
     * the tn array is non-null.  Note permutation of the xn, yn, zn, wn
     * and tn arrays, as explained above for recursion down the lo
     * branch of the tree.
     *
     * Note that for 3D the wn array isn't allocated so don't permute it
     * cyclically along with the other arrays in the recursive call.
     *
     * Note also that if SORT_ATOM_NUMBERS isn't defined the on array isn't
     * allocated so don't permute it cyclically in the recursive call.
     */

    if (upper > middle) {
      (*kdpptr)->lo = NULL;
      (*kdpptr)->hi = NULL;
      that->hi = (*kdpptr)++;

#ifndef SORT_ATOM_NUMBERS
      if (dim == 4) {
        buildkdtree(tn, yn, zn, xn, on, wn,
                    middle+1, end, kdpptr, that->hi, x, p+1);
      } else {
        buildkdtree(tn, yn, xn, wn, on, zn,
                    middle+1, end, kdpptr, that->hi, x, p+1);
      }
#else
      if (dim == 4) {
        buildkdtree(tn, yn, zn, wn, xn, on,
                    middle+1, end, kdpptr, that->hi, x, p+1);
      } else {
        buildkdtree(tn, yn, zn, xn, on, wn,
                    middle+1, end, kdpptr, that->hi, x, p+1);
      }
#endif

    }
  }
}

/***********************************************************************
                            SEARCHKDTREE()
************************************************************************/
 
/*
 * Walk the kd tree and generate the pair lists for the upper and lower
 * triangles.  The pair lists are partially ordered in descending atom
 * number if SORT_ATOM_NUMBERS is defined.  Descending order is preferred
 * by the subsequent heap sort of the pair lists that will occur if
 * HEAP_SORT_PAIRS is defined.
 *
 * Calling parameters are as follows:
 *
 * that - the node currently visited, equivalent to 'this' in C++
 * x - atomic coordinate array
 * p - the partition (x, y, z, w or o) on which sorting occurs
 * q - the query atom number
 * loindexp - pointer to pair count array index for the lower triangle
 * upindexp - pointer to pair count array index for the upper triangle
 * lopairlist - the pair list for the lower triangle
 * uppairlist - the pair list for the upper triangle
 */

static
void searchkdtree( KDNODE_T *that, REAL_T *x, int p, int q,
                   int *loindexp, int *upindexp,
                   int *lopairlist, int *uppairlist )
{
  REAL_T xij, yij, zij, wij, r2;

  /*
   * The partition cycles by dim unless SORT_ATOM_NUMBERS is defined,
   * in which case it cycles by dim+1.  Note that if SORT_ATOM_NUMBERS
   * is defined and the partition equals zero, sorting has occured
   * on the ordinal atom number instead of the atom's cartesian
   * coordinate.
   */

#ifndef SORT_ATOM_NUMBERS
  p %= dim;
#else
  p %= (dim + 1);
#endif

  /*
   * Search the high branch of the tree if the atomic coordinate of the
   * query atom plus the cutoff radius is greater than or equal to the
   * atomic coordinate of the kd node atom.
   *
   * If SORT_ATOM_NUMBERS is defined and p==0, always search the high branch.
   */

#ifdef SORT_ATOM_NUMBERS
  if ( ((p == 0) && (that->hi != NULL)) ||
       ((p != 0) && (that->hi != NULL) &&
        (x[dim * q + p-1] + cut >= x[dim * that->n + p-1])) )
#else
  if ( (that->hi != NULL) &&
       (x[dim * q + p] + cut >= x[dim * that->n + p]) )
#endif
    {
      searchkdtree(that->hi, x, p+1, q, loindexp, upindexp, lopairlist, uppairlist);
    }

  /*
   * If the query atom number does not equal the kd tree node atom number
   * and at least one of the two atoms is not frozen, calculate the interatomic
   * distance and add the kd tree node atom to one of the pair lists if the
   * distance is less than the cutoff distance.  The atom belongs on the lower
   * triangle pair list if the atom number is less than the query node atom
   * number.  Otherwise, it belongs on the upper triangle pair list.
   */

  if ( ( q != that->n ) && ( !frozen[q] || !frozen[that->n] ) ){
    xij = x[dim * q + 0] - x[dim * that->n + 0];
    yij = x[dim * q + 1] - x[dim * that->n + 1];
    zij = x[dim * q + 2] - x[dim * that->n + 2];
    r2 = xij * xij + yij * yij + zij * zij;
    if (dim == 4) {
      wij = x[dim * q + 3] - x[dim * that->n + 3];
      r2 += wij * wij;
    }
    if (r2 < cut2) {
      if (that->n < q) {
        lopairlist[*loindexp] = that->n;
        (*loindexp)++;
      } else {
        uppairlist[*upindexp] = that->n;
        (*upindexp)++;
      }
    }
  }

  /*
   * Search the low branch of the tree if the atomic coordinate of the
   * query atom minus the cutoff radius is less than the atomic coordinate
   * of the kd node atom.
   *
   * If SORT_ATOM_NUMBERS is defined and p==0, always search the low branch.
   */

#ifdef SORT_ATOM_NUMBERS
  if ( ((p == 0) && (that->lo != NULL)) ||
       ((p != 0) && (that->lo != NULL) &&
        (x[dim * q + p-1] - cut < x[dim * that->n + p-1])) )
#else
  if ( (that->lo != NULL) &&
       (x[dim * q + p] - cut < x[dim * that->n + p]) )
#endif
    {
      searchkdtree(that->lo, x, p+1, q, loindexp, upindexp, lopairlist, uppairlist);
    }
}

/***********************************************************************
                            NBLIST()
************************************************************************/
 
/*
 * Create the non-bonded pairlist using a kd-tree.  The kd-tree nodes
 * are allocated from the kdtree array.
 *
 * Calling parameters are as follows:
 *
 * lpears - the number of pairs on the lower triangle pair list
 * upears - the number of pairs on the upper triangle pair list
 * pearlist - the pair list, contiguous for the upper and lower triangles
 * x - atomic coordinate array
 * context_PxQ - the ScaLAPACK context
 * derivs - the derivative flag: -1 for 2nd derivs, 1 for 1st derivs
 *
 * This function returns the total number of pairs.
 */

static
int nblist( int *lpears, int *upears, int **pearlist, REAL_T *x,
            int context_PxQ, int derivs)
{
  int i, j, locnt, upcnt, totpair, numthreads, threadnum;
  int *xn, *yn, *zn, *wn, *on, *tn, *lopairlist, *uppairlist;
  KDNODE_T *kdtree, *kdptr, *root;

#ifdef SCALAPACK
  int myrow, mycol, nprow, npcol, lotot, uptot;
  int *lopearlist, *uppearlist, *divblk, *modrow;
#endif

  /* Square the cutoff distance for use in searchkdtree. */

  cut2 = cut * cut;

  /* Allocate the kdtree array that must hold one node per atom. */

  if ( (kdtree = (KDNODE_T *) malloc(prm->Natom*sizeof(KDNODE_T))) == NULL ) {
    fprintf( nabout,"Error allocate kdnode array in nbtree!\n");
    exit (0);
  }

  /*
   * Allocate, initialize and sort the arrays that hold the results of the
   * heapsort on x,y,z.  These arrays are used as pointers (via array indices)
   * into the atomic coordinate array x.  Allocate an additional temp array
   * so that the buildkdtree function can cycle through x,y,z.  Also allocate
   * and sort an additional array for the w coordinate if dim==4, and
   * allocate an array for the ordinal atom number if SORT_ATOM_NUMBERS is
   * defined.
   *
   * The temp array and the ordinal atom array are not sorted.
   */

  xn = ivector(0, prm->Natom);
  yn = ivector(0, prm->Natom);
  zn = ivector(0, prm->Natom);
  tn = ivector(0, prm->Natom);

  if (dim == 4) {
    wn = ivector(0, prm->Natom);
  }

#ifdef SORT_ATOM_NUMBERS
  on = ivector(0, prm->Natom);
#endif  

  for (i = 0; i < prm->Natom; i++) {
    xn[i] = yn[i] = zn[i] = i;
    if (dim == 4) {
      wn[i] = i;
    }

#ifdef SORT_ATOM_NUMBERS
    on[i] = i;
#endif

  }

  heapsort_index(xn, prm->Natom, x, 0);
  heapsort_index(yn, prm->Natom, x, 1);
  heapsort_index(zn, prm->Natom, x, 2);

  if (dim == 4) {
    heapsort_index(wn, prm->Natom, x, 3);
  }

  /*
   * Build the kd tree.  For 3D the wn array is ignored because it wasn't
   * allocated.  When SORT_ATOM_NUMBERS is not defined the on array is
   * ignored because it wasn't allocated either.  See the recursive calls
   * to the buildkdtree function from within that function to verify that
   * arrays that are ignored do not participate in the cyclic permutation
   * of arrays in the recursive calls.
   *
   * But if SORT_ATOM_NUMBERS is defined the xn, yn, zn, wn and on array order
   * is permuted in the non-recursive call to the buildkdtree function
   * (below) so that the sort at the root node of the tree occurs on the
   * ordinal atom number.
   */

  kdptr = kdtree;
  root = kdptr++;
  root->lo = NULL;
  root->hi = NULL;

#ifndef SORT_ATOM_NUMBERS

  buildkdtree(xn, yn, zn, wn, on, tn, 0, prm->Natom-1, &kdptr, root, x, 0);

#else

  buildkdtree(on, xn, yn, zn, wn, tn, 0, prm->Natom-1, &kdptr, root, x, 0);

#endif

  /*
   * Search the kd tree with each atom and record pairs into temporary
   * arrays for the lower and upper triangle pair lists for one atom.
   * Copy the temporary pair lists into a pair list array that is
   * allocated separately for each atom and that contains the lower
   * and upper triangle pair lists contiguously packed.
   *
   * The pairlist array is an array of pair list arrays.
   */

  totpair = 0;

#pragma omp parallel reduction (+: totpair) \
  private (i, j, locnt, upcnt, lopairlist, uppairlist, threadnum, numthreads)

  {
    /*
     * Get the thread number and the number of threads for multi-threaded
     * execution under OpenMP, ScaLAPACK or MPI.  These variables are not
     * used for single-threaded execution.
     *
     * If MPI is defined, the manner of assignment of the threadnum
     * and numthreads variables depends upon derivs, as follows.
     *
     * If derivs >= 0, the call to nblist is intended to build
     * a pair list for the first derivative calculation.  Therefore,
     * the threadnum and numthreads variables are assigned in
     * a row cyclic manner that is required for the first derivative
     * calculation.  This row cyclic approach makes optimal use of
     * the MPI tasks in that each task has the minimum number of
     * pair lists for the first derivative calculation.
     *
     * If derivs < 0 the call to nblist is intended to build a pair
     * list for the Born second derivative calculation.  However,
     * the Born second derivative calculation is not parallelized
     * for MPI.  Therefore, the pair list will be fully populated
     * for each MPI task.
     *
     * If OPENMP is defined, the threadnum and numthreads variables
     * are assigned in a row cyclic manner for both the first and
     * second derivative calculations.  Thus for both calculations
     * each OpenMP thread has the minimum number of pair lists.
     *
     * If SCALAPACK is defined, the manner of assignment of the threadnum
     * and numthreads variables depends upon derivs variable, as follows.
     *
     * If derivs >= 0, the call to nblist is intended to build
     * a pair list for the first derivative calculation.  Therefore,
     * the threadnum and numthreads variables are assigned in
     * a row cyclic manner that is required for the first derivative
     * calculation.  As in the MPI case, each MPI task has the minimum
     * number of pair list for the first derivative calculation.
     *
     * If derivs < 0, the call to nblist is intended to build a pair
     * list either for the Born second derivative calculation or for
     * the nonbonded second derivative calculation.  For the nonbonded
     * case, the calculation is not parallelized; therefore, the pair
     * list will be fully populated for each ScaLAPACK process.
     *
     * For the Born case, the threadnum and numthreads variables are
     * assigned from the process column and the number of process columns,
     * respectively.  Each row of a particular column receives the same
     * pair list from the search of the kd tree, but thereafter the
     * pair list is culled using the process row and number of process
     * rows so that each row ultimately receives a unique pair list.
     * Thus, although for the Born second derivative calculation each process
     * column receives more pair lists than does each MPI task for the first
     * derivative calculation, the pair lists are shorter, and therefore the
     * total number of pairs that are processed by each process column is the
     * same as in the first derivative calculation, with the following.
     * caveat.  For processes that do not lie of the process grid, the
     * process column is -1 which will result in no pair list for that process
     * due to the fact that the myroc function returns 0.  Hence for the
     * second derivative calculation, fewer processes have pair lists than
     * for the first derivative calculation unless the total number of
     * MPI tasks is a square such as 1, 4, 9, et cetera.
     *
     */

#if defined(OPENMP)

    threadnum = omp_get_thread_num();
    numthreads = omp_get_num_threads();

#elif defined(SCALAPACK)

    if (derivs < 0) {
      if (gb) {
        blacs_gridinfo_(&context_PxQ, &nprow, &npcol, &myrow, &mycol);
        threadnum = mycol;
        numthreads = npcol;
      } else {
        threadnum = 0;
        numthreads = 1;
      }
    } else {
      threadnum = mytaskid;
      numthreads = numtasks;
    }

#elif defined(MPI)

    if (derivs < 0) {
      threadnum = 0;
      numthreads = 1;
    } else {
      threadnum = mytaskid;
      numthreads = numtasks;
    }

#endif

    /*
     * Allocate the temporary arrays for the lower and upper triangles.
     * These arrays must be large enough to hold a maximum size pair
     * list for one atom.  For the ScaLAPACK second derivatives an
     * extra set of temporary arrays is used from which the final pair
     * lists are culled.
     *
     * Also, for ScaLAPACK allocate and initialize lookup tables
     * for division and modulus operations.
     */

    lopairlist = ivector( 0, prm->Natom );
    uppairlist = ivector( 0, prm->Natom );

#ifdef SCALAPACK

    if (derivs < 0 && gb) {
      lopearlist = ivector( 0, prm->Natom );
      uppearlist = ivector( 0, prm->Natom );

      divblk = ivector( 0, prm->Natom );
      modrow = ivector( 0, prm->Natom );

      for (i = 0; i < prm->Natom; i++) {
        divblk[i] = i / BLOCKSIZE;
        modrow[i] = i % nprow;
      }
    }

#endif

    /*
     * Search the kd tree with each atom.  If no pair list array
     * has been allocated and there are pair atoms, allocate a
     * pair list array.  If a pair list array has been allocated
     * but it is too small for the number of pair atoms, deallocate
     * the pair list array and allocate a larger array.  If it is
     * at least 33% larger than is necessary for the pair atoms,
     * deallocate it and allocate a smaller pair list array.  Copy the
     * lower and upper triangle pair lists into the pair list array,
     * packed contiguously with the lower triangle pair list followed
     * by the upper triangle pair list.
     *
     * Explicitly assign threads to loop indices for the following loop,
     * in a manner equivalent to (static, N) scheduling with OpenMP,
     * and identical to the manner in which threads are assigned in
     * nbond, egb and egb2.
     *
     * There is an implied barrier end of this OpenMP parallel region.
     * Because the following loop is the only loop in this parallel
     * region, there is no need for an explicit barrier.  Furthermore,
     * all OpenMP parallelized loops that use the pair list also use
     * loop index to thread mapping that is identical to what is used for
     * this loop.  Hence, no race condition would exist even if OpenMP
     * threads were not synchronized at the end of this parallel region
     * because each thread constructs the pair list that it uses thereafter.
     * This same argument applies to MPI tasks: no synchronization is
     * necessary.
     */

    for (i = 0; i < prm->Natom; i++) {

#if defined(OPENMP) || defined(MPI) || defined(SCALAPACK)

      if ( !myroc(i, BLOCKSIZE, numthreads, threadnum) ) continue;

#endif

#ifdef SCALAPACK

      if (derivs < 0 && gb) {
        lotot = uptot = 0;
        searchkdtree(root, x, 0, i, &lotot, &uptot, lopearlist, uppearlist);
      } else {
        locnt = upcnt = 0;
        searchkdtree(root, x, 0, i, &locnt, &upcnt, lopairlist, uppairlist);
      }

#else

      locnt = upcnt = 0;
      searchkdtree(root, x, 0, i, &locnt, &upcnt, lopairlist, uppairlist);

#endif

      /*
       * If SORT_ATOM_NUMBERS is defined, the upper and lower triangle
       * pair lists are partially sorted by ordinal atom number using
       * the kd tree.
       *
       * If HEAP_SORT_PAIRS is defined, sort the upper and lower triangle
       * pair lists using heap sort.
       *
       * If the pair lists are sorted by ordinal atom number using the
       * kd tree, the subsequent heap sort of the pair lists is quicker,
       * but the kd tree sort is not necessary.
       */

#define HEAP_SORT_PAIRS
#ifdef HEAP_SORT_PAIRS

#ifdef SCALAPACK

      if (derivs < 0 && gb) {
        heapsort_pairs(lopearlist, lotot);
        heapsort_pairs(uppearlist, uptot);
      } else {
        heapsort_pairs(lopairlist, locnt);
        heapsort_pairs(uppairlist, upcnt);
      }

#else

      heapsort_pairs(lopairlist, locnt);
      heapsort_pairs(uppairlist, upcnt);

#endif

#endif

      /*
       * For the ScaLAPACK second derivatives cull the pair lists
       * by copying to the final pair lists only those atoms that
       * are active in a particular process row.  Use the lookup
       * tables for a faster form of calls to myroc of the form:
       *
       *       myroc(lopearlist[j], BLOCKSIZE, npcol, mycol)
       */

#ifdef SCALAPACK

      if (derivs < 0 && gb) {
        locnt = 0;
        for (j = 0; j < lotot; j++) {
          if ( myrow >= 0 && modrow[divblk[lopearlist[j]]] == myrow ) {
            lopairlist[locnt++] = lopearlist[j];
          }
        }
        upcnt = 0;
        for (j = 0; j < uptot; j++) {
          if ( myrow >= 0 && modrow[divblk[uppearlist[j]]] == myrow ) {
            uppairlist[upcnt++] = uppearlist[j];
          }
        }
      }

#endif

      if ( ( pearlist[i] == NULL ) && ( locnt + upcnt > 0 ) ) {
        pearlist[i] = ivector( 0, locnt + upcnt );
      } else if ( ( pearlist[i] != NULL ) &&
                  ( ( locnt + upcnt > lpears[i] + upears[i] ) ||
                    ( 4*(locnt + upcnt) < 3*(lpears[i] + upears[i]) ) ) ) {
        free_ivector( pearlist[i], 0, lpears[i] + upears[i] );
        pearlist[i] = ivector( 0, locnt + upcnt );
      }
      lpears[i] = locnt;
      upears[i] = upcnt;
      for (j = 0; j < locnt; j++) {
        pearlist[i][j] = lopairlist[j];
      }
      for (j = 0; j < upcnt; j++) {
        pearlist[i][locnt + j] = uppairlist[j];
      }
      totpair += locnt + upcnt;
    }

    /*
     * Deallocate the temporary arrays for the lower and upper triangles.
     * For ScaLAPACK deallocate the addtional temporary arrays as well as
     * the lookup tables.
     */

    free_ivector( lopairlist, 0, prm->Natom );
    free_ivector( uppairlist, 0, prm->Natom );

#ifdef SCALAPACK

    if (derivs < 0 && gb) {
      free_ivector( lopearlist, 0, prm->Natom );
      free_ivector( uppearlist, 0, prm->Natom );

      free_ivector( divblk, 0, prm->Natom );
      free_ivector( modrow, 0, prm->Natom );
    }

#endif

  }

  /* Free the temporary arrays. */

  free(kdtree);
  free_ivector(xn, 0, prm->Natom);
  free_ivector(yn, 0, prm->Natom);
  free_ivector(zn, 0, prm->Natom);
  free_ivector(tn, 0, prm->Natom);

  if (dim == 4) {
    free_ivector(wn, 0, prm->Natom);
  }

#ifdef SORT_ATOM_NUMBERS

  free_ivector(on, 0, prm->Natom);

#endif

  return totpair;
}

/***********************************************************************
                            NBOND()
************************************************************************/

/* 
 * Calculate the non-bonded energy and first derivatives.
 * This function is complicated by the fact that it must
 * process two forms of pair lists: the 1-4 pair list and
 * the non-bonded pair list.  The non-bonded pair list
 * must be modified by the excluded atom list whereas the
 * 1-4 pair list is used unmodified.  Also, the non-bonded
 * pair list comprises lower and upper triangles whereas
 * the 1-4 pair list comprises an upper triangle only.
 *
 * Calling parameters are as follows:
 *
 * lpears - the number of pairs on the lower triangle pair list
 * upears - the number of pairs on the upper trianble pair list
 * pearlist - either the 1-4 pair list or the non-bonded pair list
 * N14 - set to 0 for the non-bonded pair list, 1 for the 1-4 pair list
 * x - the atomic coordinate array
 * f - the gradient vector
 * enb - Van der Waals energy return value, passed by reference
 * eel - Coulombic energy return value, passed by reference
 * enbfac - scale factor for Van der Waals energy
 * eelfac - scale factor for Coulombic energy
 */

static
int nbond( int *lpears, int *upears, int **pearlist, int N14,
           REAL_T *x, REAL_T *f, REAL_T *enb, REAL_T *eel,
           REAL_T enbfac, REAL_T eelfac )
{
  int i, j, jn, ic, npr, lpair, iaci, foff, threadnum, numthreads;
  int *iexw;
  REAL_T dumx, dumy, dumz, dumw, cgi, r2inv, df2, r6, r10, f1, f2;
  REAL_T dedx, dedy, dedz, dedw, df, enbfaci, eelfaci, evdw, elec;
  REAL_T xi, yi, zi, wi, xij, yij, zij, wij, r, r2;
  REAL_T dis, kij, d0, diff, rinv, rs, rssq, eps1, epsi, cgijr, pow;
  int nhbpair, ibig, isml;

#define SIG 0.3
#define DIW 78.0
#define C1 38.5

  evdw = 0.;
  elec = 0.;
  enbfaci = 1. / enbfac;
  eelfaci = 1. / eelfac;

#pragma omp parallel reduction (+: evdw, elec) \
  private (i, j, iexw, npr, iaci, \
           xi, yi, zi, wi, xij, yij, zij, wij, dumx, dumy, dumz, dumw, \
           cgi, jn, r2, r2inv, r, rinv, rs, rssq, pow, \
           eps1, epsi, cgijr, df2, ic, r6, f2, f1, df, dis, d0, kij, \
           diff, ibig, isml, dedx, dedy, dedz, dedw, r10, \
           threadnum, numthreads, foff)
  {
    /*
     * Get the thread number and the number of threads for multi-threaded
     * execution under OpenMP.  For all other cases, including ScaLAPACK,
     * MPI and single-threaded execution, use the values that have been
     * stored in mytaskid and numtasks, respectively.
     */

#if defined(OPENMP)
    threadnum = omp_get_thread_num();
    numthreads = omp_get_num_threads();
#else
    threadnum = mytaskid;
    numthreads = numtasks;
#endif

    /*
     * Compute an offset into the gradient array for this thread,
     * but only if OPENMP is defined.
     */

#ifdef OPENMP
    foff = dim * prm->Natom * threadnum;
#else
    foff = 0;
#endif

    /*
     * Allocate and initialize the iexw array used for skipping excluded
     * atoms.  Note that because of the manner in which iexw is used, it
     * is necessary to initialize it before only the first iteration of
     * the following loop.
     */
    
    iexw = ivector(-1, prm->Natom);
    for (i = -1; i < prm->Natom; i++) {
      iexw[i] = -1;
    }

    /*
     * Loop over all atoms i except for the final atom i.
     *
     * Explicitly assign threads to loop indices for the following loop,
     * in a manner equivalent to (static, N) scheduling with OpenMP, and
     * identical to the manner in which threads are assigned in nblist.
     *
     * Synchronization of OpenMP threads will occur following this loop
     * because the parallel region ends after this loop.  Following
     * synchronization, a reduction of the sumdeijda array will be
     * performed.
     *
     * Synchronization of MPI tasks will occur via the MPI_Allreduce
     * function that is called from within mme34.
     */

    for (i = 0; i < prm->Natom - 1; i++) {

#if defined(OPENMP) || defined(MPI) || defined(SCALAPACK)
      if ( !myroc(i, BLOCKSIZE, numthreads, threadnum) ) continue;
#endif

      /* Check whether there are any atoms j on the pair list of atom i. */

      npr = upears[i];
      if (npr <= 0) continue;

      iaci = prm->Ntypes * (prm->Iac[i] - 1);
      dumx = dumy = dumz = 0.0;

      xi = x[dim * i + 0];
      yi = x[dim * i + 1];
      zi = x[dim * i + 2];

      if (dim == 4) {
        dumw = 0.0;
        wi = x[dim * i + 3];
      }

      cgi = eelfaci * prm->Charges[i];

      /*
       * Expand the excluded list into the iexw array by storing i
       * at array address j.
       */

      for (j = 0; j < prm->Iblo[i]; j++) {
        iexw[IexclAt[i][j] - 1] = i;
      }

      /*
       * If the 'N14' calling parameter is clear, use the beginning
       * address of the upper triangle pair list, which happens
       * to be the number of atoms on the lower triangle pair list.
       * If the 'N14' calling parameter is set, the beginning
       * address is zero because no lower triangle pair list is
       * used for the N14 interactions.
       */

      if (N14 == 0) {
        lpair = lpears[i];
      } else {
        lpair = 0;
      }

      /* Select atoms j from the pair list.  Non-graceful error handling. */

      for (jn = 0; jn < npr; jn++) {

        if (pearlist[i] == NULL) {
          fprintf( nabout,"NULL pair list entry in nbond, taskid = %d\n", mytaskid);
          fflush(nabout);
        }
        j = pearlist[i][lpair + jn];

        /*
         * If the 'N14' calling parameter is clear, check whether
         * this i,j pair is exempted by the excluded atom list.
         */

        if (N14 != 0 || iexw[j] != i) {
          xij = xi - x[dim * j + 0];
          yij = yi - x[dim * j + 1];
          zij = zi - x[dim * j + 2];
          r2 = xij * xij + yij * yij + zij * zij;

          if (dim == 4) {
            wij = wi - x[dim * j + 3];
            r2 +=  wij * wij;
          }

          r2inv = 1.0 / r2;
          r = sqrt(r2);
          rinv = r * r2inv;

          /* Calculate the energy and derivatives according to dield. */

          if (dield == -3) {

            /* special code Ramstein & Lavery dielectric, 94 force field */

            rs = SIG * r;
            rssq = rs * rs;
            pow = exp(-rs);
            eps1 = rssq + rs + rs + 2.0;
            epsi = 1.0 / (DIW - C1 * pow * eps1);
            cgijr = cgi * prm->Charges[j] * rinv * epsi;
            elec += cgijr;
            df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
            ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
            if( ic >= 0 ){
              r6 = r2inv * r2inv * r2inv;
              f2 = prm->Cn2[ic] * r6;
              f1 = prm->Cn1[ic] * r6 * r6;
              evdw += (f1 - f2) * enbfaci;
              df = (df2 + (6.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
            } else {
              df = df2 * rinv;
            }

          } else if (dield == -4) {

            /* distance-dependent dielectric code, 94 ff */
            /* epsilon = r  */

            rs = cgi * prm->Charges[j] * r2inv;
            df2 = -2.0 * rs;
            elec += rs;
            ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
            if( ic >= 0 ){
              r6 = r2inv * r2inv * r2inv;
              f2 = prm->Cn2[ic] * r6;
              f1 = prm->Cn1[ic] * r6 * r6;
              evdw += (f1 - f2) * enbfaci;
              df = (df2 + (6.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
            } else {
              df = df2 * rinv;
            }

          } else if (dield == -5) {

            /* non-bonded term from yammp  */

            dis = r;
            ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
            d0 = prm->Cn2[ic];
            if (dis < d0) {
              kij = prm->Cn1[ic];
              diff = dis - d0;
              evdw += kij * diff * diff;
              df = 2.0 * kij * diff;
            } else {
              df = 0.0;
            }
          } else {

            /*
             * Code for various dielectric models.
             * The df2 variable should hold r(dV/dr).
             */

            if (dield == 0) {

              /* epsilon = r  */

              rs = cgi * prm->Charges[j] * r2inv;
              df2 = -2.0 * rs;
              elec += rs;

            } else if (dield == 1) {
              
              /* epsilon = 1  */

              rs = cgi * prm->Charges[j] * rinv;
              df2 = -rs;
              elec += rs;

            } else if (dield == -2) {

              /* Ramstein & Lavery dielectric, PNAS 85, 7231 (1988). */

              rs = SIG * r;
              rssq = rs * rs;
              pow = exp(-rs);
              eps1 = rssq + rs + rs + 2.0;
              epsi = 1.0 / (DIW - C1 * pow * eps1);
              cgijr = cgi * prm->Charges[j] * rinv * epsi;
              elec += cgijr;
              df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
            }

            /* Calculate either Van der Waals or hydrogen bonded term. */

            ic = prm->Cno[iaci + prm->Iac[j] - 1];
            if (ic > 0 || enbfac != 1.0) {
              if (ic > 0) {
                ic--;
              } else {
                ibig = prm->Iac[i] > prm->Iac[j] ?
                  prm->Iac[i] : prm->Iac[j];
                isml = prm->Iac[i] > prm->Iac[j] ?
                  prm->Iac[j] : prm->Iac[i];
                ic = ibig * (ibig - 1) / 2 + isml - 1;
              }
              r6 = r2inv * r2inv * r2inv;
              f2 = prm->Cn2[ic] * r6;
              f1 = prm->Cn1[ic] * r6 * r6;
              evdw += (f1 - f2) * enbfaci;
              df = (df2 + (6.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
#if 0
              if( enbfac != 1.0 ) nb14 += (f1-f2) * enbfaci;
#endif
            } else {
              ic = -ic - 1;
              r10 = r2inv * r2inv * r2inv * r2inv * r2inv;
              f2 = prm->HB10[ic] * r10;
              f1 = prm->HB12[ic] * r10 * r2inv;
              evdw += (f1 - f2) * enbfaci;
              df = (df2 + (10.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
              nhbpair++;
#if 0
              hbener += (f1-f2) * enbfaci;
#endif
            }
          }

          /*
           * Divide by r here instead dividing the derivatives by r,
           * and update the gradient for atom j.
           */

          df *= rinv;

          dedx = df * xij;
          dedy = df * yij;
          dedz = df * zij;

          dumx += dedx;
          dumy += dedy;
          dumz += dedz;

          f[foff + dim * j + 0] -= dedx;
          f[foff + dim * j + 1] -= dedy;
          f[foff + dim * j + 2] -= dedz;

          if (dim == 4) {
            dedw = df * wij;
            dumw += dedw;
            f[foff + dim * j + 3] -= dedw;
          }
        }
      }

      /* For atom i, the gradient is updated in the i-loop only. */

      f[foff + dim * i + 0] += dumx;
      f[foff + dim * i + 1] += dumy;
      f[foff + dim * i + 2] += dumz;

      if (dim == 4) {
        f[foff + dim * i + 3] += dumw;
      }
    }

    /* Deallocate the iexw array within this potentially parallel region. */

    free_ivector(iexw, -1, prm->Natom);
  }

  /* Return evdw and elec through by-reference calling parameters. */

  *enb = evdw;
  *eel = elec;

  return (0);
}

/***********************************************************************
                            EGB()
************************************************************************/

/*
 * Calculate the generalized Born energy and first derivatives.
 *
 * Calling parameters are as follows:
 *
 * lpears - the number of pairs on the lower triangle pair list
 * upears - the number of pairs on the upper trianble pair list
 * pearlist - the pair list, contiguous for the upper and lower triangles
 * x - input: the atomic (x,y,z) coordinates
 * f - updated: the gradient vector
 * fs - input: overlap parameters
 * rborn - input: atomic radii
 * q - input: atomic charges
 * kappa - input: inverse of the Debye-Huckel length
 * diel_ext - input: solvent dielectric constant
 * enb - updated: Lennard-Jones energy
 * eelt - updated: gas-phase electrostatic energy
 * freevectors - if !=0 free the static vectors and return
 */

static
REAL_T egb( int *lpears, int *upears, int **pearlist, REAL_T *x, REAL_T *f,
            REAL_T *fs, REAL_T *rborn, REAL_T *q,
            REAL_T *kappa, REAL_T *diel_ext, REAL_T *enb, REAL_T *eelt,
            INT_T freevectors )

#define BOFFSET 0.09
#define KSCALE 0.73

#define PIx4 12.5663706143591724639918
#define PIx2  6.2831853071795862319959
#define PIx1  3.1415926535897931159979

{
#if defined(MPI) || defined(SCALAPACK)
  int ierror;
  static REAL_T *reductarr=NULL;
#endif

  static REAL_T *reff=NULL, *sumdeijda=NULL, *psi=NULL;

  int i, i34, j, j34, k, threadnum, numthreads, maxthreads, foff, soff;
  int npairs, ic, iaci;
  int *iexw;
  size_t n;
  REAL_T epol, dielfac, qi, qj, qiqj, fgbi, fgbk, rb2, expmkf;
  REAL_T elec, evdw, sumda, daix, daiy, daiz, daiw;
  REAL_T xi, yi, zi, wi, xij, yij, zij, wij;
  REAL_T dedx, dedy, dedz, dedw, de;
  REAL_T dij1i, dij3i, temp1;
  REAL_T qi2h, qid2h, datmp;
  REAL_T theta, ri1i, dij2i;

  REAL_T dij, sumi;
  REAL_T eel, f6, f12, rinv, r2inv, r6inv;
  REAL_T r2, ri, rj, sj, sj2, thi;
  REAL_T uij, efac, temp4, temp5, temp6;

  /* LCPO stuff follows */
  int count, count2, icount;
  REAL_T si, sumAij, sumAjk, sumAijAjk, rij, tmpaij, Aij, dAijddij;
  REAL_T sumdAijddijdxi, sumdAijddijdyi, sumdAijddijdzi, sumdAijddijdwi;
  REAL_T sumdAijddijdxiAjk, sumdAijddijdyiAjk, sumdAijddijdziAjk, sumdAijddijdwiAjk;
  REAL_T dAijddijdxj, dAijddijdyj, dAijddijdzj, dAijddijdwj;
  REAL_T sumdAjkddjkdxj, sumdAjkddjkdyj, sumdAjkddjkdzj, sumdAjkddjkdwj, p3p4Aij;
  REAL_T xk, yk, zk, wk, rjk2, djk1i, rjk, vdw2dif, tmpajk, Ajk, sumAjk2, dAjkddjk;
  REAL_T dAjkddjkdxj, dAjkddjkdyj, dAjkddjkdzj, dAjkddjkdwj;
  REAL_T lastxj, lastyj, lastzj, lastwj;
  REAL_T dAidxj, dAidyj, dAidzj, dAidwj, Ai, dAidxi, dAidyi, dAidzi, dAidwi;
  REAL_T totsasa;
  REAL_T xj, yj, zj, wj;
  REAL_T dumbo,tmpsd;
  REAL_T rgbmax1i, rgbmax2i, rgbmaxpsmax2;

  /*FGB taylor coefficients follow */
  /* from A to H :                 */
  /* 1/3 , 2/5 , 3/7 , 4/9 , 5/11  */
  /* 4/3 , 12/5 , 24/7 , 40/9 , 60/11 */

#define TA 0.33333333333333333333
#define TB 0.4
#define TC 0.42857142857142857143
#define TD 0.44444444444444444444
#define TDD 0.45454545454545454545

#define TE 1.33333333333333333333
#define TF 2.4
#define TG 3.42857142857142857143
#define TH 4.44444444444444444444
#define THH 5.45454545454545454545

  /*
   * Determine the size of the sumdeijda array.  If OPENMP is defined,
   * a copy of this array must be allocated for each thread; otherwise
   * only one copy is allocated.
   */

#ifdef OPENMP
  maxthreads = omp_get_max_threads();
#else
  maxthreads = 1;
#endif

  n = (size_t)prm->Natom;

  /*
   * If freevectors != 0, deallocate the static arrays that have been
   * previously allocated and return.
   */

  if (freevectors != 0) {
    if (reff != NULL) free_vector(reff, 0, n);
    reff = NULL;
    if (sumdeijda != NULL) free_vector(sumdeijda, 0, maxthreads*n);
    sumdeijda = NULL;
    if (psi != NULL) free_vector(psi, 0, n);
    psi = NULL;
#if defined(MPI) || defined(SCALAPACK)
    if (reductarr != NULL) free_vector(reductarr, 0, n);
    reductarr = NULL;
#endif
    return (0.0);
  }

  /*
   * Smooth "cut-off" in calculating GB effective radii.
   * Implementd by Andreas Svrcek-Seiler and Alexey Onufriev.
   * The integration over solute is performed up to rgbmax and includes
   * parts of spheres; that is an atom is not just "in" or "out", as
   * with standard non-bonded cut.  As a result, calclated effective
   * radii are less than rgbmax. This saves time, and there is no
   * discontinuity in dReff/drij.
   *
   * Only the case rgbmax > 5*max(sij) = 5*fsmax ~ 9A is handled; this is
   * enforced in mdread().  Smaller values would not make much physical
   * sense anyway.
   *
   * Note: rgbmax must be less than or equal to cut so that the pairlist
   * generated from cut may be applied to calculation of the effective
   * radius and its derivatives.
   */

  if (rgbmax > cut) {
    fprintf( nabout,"Error in egb: rgbmax = %f is greater than cutoff = %f\n",
           rgbmax, cut);
    exit(1);
  }

  rgbmax1i = 1.0/rgbmax;
  rgbmax2i = rgbmax1i*rgbmax1i;
  rgbmaxpsmax2 = (rgbmax+prm->Fsmax)*(rgbmax+prm->Fsmax);

  /* Allocate some static arrays if they have not been allocated already. */

  if (reff == NULL) reff = vector(0, n);
  if (sumdeijda == NULL) sumdeijda = vector(0, maxthreads*n);
  if ( (psi == NULL) && (gb==2 || gb==5) ) psi = vector(0, n);
#if defined(MPI) || defined(SCALAPACK)
  if (reductarr == NULL) reductarr = vector(0, n);
#endif


  count = 0;
  totsasa = 0.0;

  if (gb_debug)
    fprintf( nabout,"Effective Born radii:\n");

  /* 
   * Get the "effective" Born radii via the approximate pairwise method
   * Use Eqs 9-11 of Hawkins, Cramer, Truhlar, J. Phys. Chem. 100:19824
   * (1996).
   *
   * This computation is not executed in parallel for the case of gbsa==1
   * due to the loop-carried dependence of the 'count' variable.  However,
   * an apparent bug in the Sun implementation of OpenMP may cause this
   * computation to execute more slowly than in serial mode in the case
   * where an 'if (gbsa != 1)' clause is added to the following pragma.
   * To avoid this problem, the code is replicated verbatim and selected
   * with an 'if (gbsa != 1)' statement.  When gbsa != 1 the pragma version
   * is selected; otherwise, the non-pragma version is selected.
   */

  if (gbsa != 1) {

    /*
     * For MPI or ScaLAPACK, initialize all elements of the reff array.
     * Although each task will calculate only a subset of the elements,
     * a reduction is used to combine the results from all tasks.
     * If a gather were used instead of a reduction, no initialization
     * would be necessary.
     */

#if defined(MPI) || defined(SCALAPACK)
    for (i = 0; i < prm->Natom; i++) {
      reff[i] = 0.0;
    }
#endif

#pragma omp parallel \
  private (i, xi, yi, zi, wi, ri, ri1i, sumi, j, k, xij, yij, zij, wij, \
           r2, dij1i, dij, sj, sj2, uij, dij2i, tmpsd, dumbo, theta, \
           threadnum, numthreads)
    {

      /*
       * Get the thread number and the number of threads for multi-threaded
       * execution under OpenMP.  For all other cases, including ScaLAPACK,
       * MPI and single-threaded execution, use the values that have been
       * stored in mytaskid and numtasks, respectively.
       */

#if defined(OPENMP)
      threadnum = omp_get_thread_num();
      numthreads = omp_get_num_threads();
#else
      threadnum = mytaskid;
      numthreads = numtasks;
#endif

      /*
       * Loop over all atoms i.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       *
       * The reff array is written in the following loops.  It is necessary to
       * synchronize the OpenMP threads or MPI tasks that execute these loops
       * following loop execution so that a race condition does not exist for
       * reading the reff array before it is written.  Even if all subsequent
       * loops use loop index to thread or task mapping that is identical to
       * that of the following loop, elements of the reff array are indexed by
       * other loop indices, so synchronization is necessary.
       *
       * OpenMP synchronization is accomplished by the implied barrier
       * at the end of this parallel region.  MPI synchronization is
       * accomplished by MPI_Allreduce.
       */

      for (i = 0; i < prm->Natom; i++) {

#if defined(OPENMP) || defined(MPI) || defined(SCALAPACK)
        if ( !myroc(i, BLOCKSIZE, numthreads, threadnum) ) continue;
#endif

        xi = x[dim * i    ];
        yi = x[dim * i + 1];
        zi = x[dim * i + 2];

        if (dim == 4) {
          wi = x[dim * i + 3];
        }

        ri = rborn[i] - BOFFSET;
        ri1i = 1. / ri;
        sumi = 0.0;

        /* Select atom j from the pair list.  Non-graceful error handling. */

        for (k = 0; k < lpears[i] + upears[i]; k++) {

          if (pearlist[i] == NULL) {
            fprintf( nabout,"NULL pair list entry in egb loop 1, taskid = %d\n", mytaskid);
            fflush(nabout);
          }
          j = pearlist[i][k];

          xij = xi - x[dim * j    ];
          yij = yi - x[dim * j + 1];
          zij = zi - x[dim * j + 2];
          r2 = xij * xij + yij * yij + zij * zij;

          if (dim == 4) {
            wij = wi - x[dim * j + 3];
            r2 += wij * wij;
          }

          if ( r2 > rgbmaxpsmax2 ) continue;
          dij1i = 1.0 / sqrt(r2);
          dij = r2 * dij1i;
          sj = fs[j] * (rborn[j] - BOFFSET);
          sj2 = sj * sj;

          /*
           * ---following are from the Appendix of Schaefer and Froemmel,
           * JMB 216:1045-1066, 1990;  Taylor series expansion for d>>s
           * is by Andreas Svrcek-Seiler; smooth rgbmax idea is from
           * Andreas Svrcek-Seiler and Alexey Onufrie.
           */

          if( dij > rgbmax + sj ) continue;

          if ((dij > rgbmax - sj)) {
            uij = 1./(dij -sj);
            sumi -= 0.125 * dij1i * (1.0 + 2.0 * dij *uij +
                                     rgbmax2i * (r2 - 4.0 * rgbmax * dij - sj2) +
                                     2.0 * log((dij-sj)*rgbmax1i));

          } else if (dij >4.0*sj) {
            dij2i = dij1i*dij1i;
            tmpsd = sj2*dij2i;
            dumbo = TA+tmpsd*(TB+tmpsd*(TC+tmpsd*(TD+tmpsd*TDD)));
            sumi -= sj*tmpsd*dij2i*dumbo;

          } else if (dij > ri + sj) {
            sumi -= 0.5 * (sj / (r2 - sj2) +
                           0.5 * dij1i * log((dij - sj) / (dij + sj)));

          } else if (dij > fabs(ri - sj)) {
            theta = 0.5 * ri1i * dij1i * (r2 + ri * ri - sj2);
            uij = 1. / (dij + sj);
            sumi -= 0.25 * (ri1i * (2. - theta) - uij +
                            dij1i * log(ri * uij));

          } else if (ri < sj) {
            sumi -= 0.5 * (sj / (r2 - sj2) + 2. * ri1i +
                           0.5 * dij1i * log((sj - dij) / (sj + dij)));

          }

        }

        if( gb==1 ) {

          /* "standard" (HCT) effective radii:  */
          reff[i] = 1.0/(ri1i + sumi);
          if( reff[i] < 0.0 ) reff[i] = 30.0;

        } else {

          /* "gbao" formulas:  */

          psi[i] = -ri*sumi;
          reff[i] = 1.0 / (ri1i - tanh( (gbalpha - gbbeta*psi[i] + 
                                         gbgamma*psi[i]*psi[i]) *psi[i])/rborn[i] );
        }

        if (gb_debug)
          fprintf( nabout,"%d\t%15.7f\t%15.7f\n", i + 1, rborn[i], reff[i]);
      }
    }

    /* The MPI synchronization is accomplished via reduction of the reff array. */

#if defined(MPI) || defined(SCALAPACK)

      ierror = MPI_Allreduce(reff, reductarr, prm->Natom,
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (ierror != MPI_SUCCESS) {
        fprintf( nabout,"Error in egb reff reduction, error = %d  mytaskid = %d\n",
               ierror, mytaskid);
      }
      for (i = 0; i < prm->Natom; i++) {
        reff[i] = reductarr[i];
      }

#endif

  } else {

    /*
     * Here is the second replication of the code, always executed serially.
     * Loop over all atoms i.
     */

    for (i = 0; i < prm->Natom; i++) {

      xi = x[dim * i    ];
      yi = x[dim * i + 1];
      zi = x[dim * i + 2];

      if (dim == 4) {
        wi = x[dim * i + 3];
      }

      ri = rborn[i] - BOFFSET;
      ri1i = 1. / ri;
      sumi = 0.0;

      /* Select atom j from the pair list.  Non-graceful error handling. */

      for (k = 0; k < lpears[i] + upears[i]; k++) {

        if (pearlist[i] == NULL) {
          fprintf( nabout,"NULL pair list entry in egb loop 2, taskid = %d\n", mytaskid);
          fflush(nabout);
        }
        j = pearlist[i][k];

        xij = xi - x[dim * j    ];
        yij = yi - x[dim * j + 1];
        zij = zi - x[dim * j + 2];
        r2 = xij * xij + yij * yij + zij * zij;

        if (dim == 4) {
          wij = wi - x[dim * j + 3];
          r2 += wij * wij;
        }

        if ( r2 > rgbmaxpsmax2 ) continue;
        dij1i = 1.0 / sqrt(r2);
        dij = r2 * dij1i;
        sj = fs[j] * (rborn[j] - BOFFSET);
        sj2 = sj * sj;

        /*
         * ---following are from the Appendix of Schaefer and Froemmel,
         * JMB 216:1045-1066, 1990;  Taylor series expansion for d>>s
         * is by Andreas Svrcek-Seiler; smooth rgbmax idea is from
         * Andreas Svrcek-Seiler and Alexey Onufrie.
         */

        if( dij > rgbmax + sj ) continue;

        if ((dij > rgbmax - sj)) {
          uij = 1./(dij -sj);
          sumi -= 0.125 * dij1i * (1.0 + 2.0 * dij *uij +
                                   rgbmax2i * (r2 - 4.0 * rgbmax * dij - sj2) +
                                   2.0 * log((dij-sj)*rgbmax1i));

        } else if (dij >4.0*sj) {
          dij2i = dij1i*dij1i;
          tmpsd = sj2*dij2i;
          dumbo = TA+tmpsd*(TB+tmpsd*(TC+tmpsd*(TD+tmpsd*TDD)));
          sumi -= sj*tmpsd*dij2i*dumbo;

        } else if (dij > ri + sj) {
          sumi -= 0.5 * (sj / (r2 - sj2) +
                         0.5 * dij1i * log((dij - sj) / (dij + sj)));

        } else if (dij > fabs(ri - sj)) {
          theta = 0.5 * ri1i * dij1i * (r2 + ri * ri - sj2);
          uij = 1. / (dij + sj);
          sumi -= 0.25 * (ri1i * (2. - theta) - uij +
                          dij1i * log(ri * uij));

        } else if (ri < sj) {
          sumi -= 0.5 * (sj / (r2 - sj2) + 2. * ri1i +
                         0.5 * dij1i * log((sj - dij) / (sj + dij)));

        }

        if (gbsa == 1) {
          if ((P0[i] + P0[j]) > dij) {
            if (P0[i] > 2.5 && P0[j] > 2.5) {
              ineighbor[count++] = j + 1;
            }
          }
        }

      }

      if( gb==1 ) {

        /* "standard" (HCT) effective radii:  */
        reff[i] = 1.0/(ri1i + sumi);
        if( reff[i] < 0.0 ) reff[i] = 30.0;

      } else {

        /* "gbao" formulas:  */

        psi[i] = -ri*sumi;
        reff[i] = 1.0 / (ri1i - tanh( (gbalpha - gbbeta*psi[i] + 
                                       gbgamma*psi[i]*psi[i]) *psi[i])/rborn[i] );
      }

      if (gb_debug)
        fprintf( nabout,"%d\t%15.7f\t%15.7f\n", i + 1, rborn[i], reff[i]);

      /* gbsa==1 */

      ineighbor[count++] = 0;
    }
    if (gb_debug) {
      fprintf( nabout, "Natom,count\t%d\t%d\t\n", prm->Natom, count );
      for (i = 0; i < prm->Natom; i++) {
        fprintf( nabout, "%d\t%d\t%15.7f\t%15.7f\n", i, ineighbor[i],
                 x[dim * i], x[dim * i +1]);
      }
    }
  }
  assert( count <= INEIGHBOR_SIZE * prm->Natom );

  /*
   * Main LCPO stuff follows.  No parallelization using OpenMP due to
   * the loop-carried dependence of the 'count variable.
   */

  esurf = 0.0;
  if (gbsa == 1) {
    totsasa = 0.0;
    count = 0;
    for (i = 0; i < prm->Natom; i++) {
      if (ineighbor[count] == 0) {
        count = count + 1;
      } else {

        si = PIx4 * P0[i] * P0[i];
        sumAij = 0.0;
        sumAjk = 0.0;
        sumAjk2 = 0.0;
        sumAijAjk = 0.0;

        sumdAijddijdxi = 0.0;
        sumdAijddijdyi = 0.0;
        sumdAijddijdzi = 0.0;
        sumdAijddijdwi = 0.0;

        sumdAijddijdxiAjk = 0.0;
        sumdAijddijdyiAjk = 0.0;
        sumdAijddijdziAjk = 0.0;
        sumdAijddijdwiAjk = 0.0;

        icount = count;

        L70:j = ineighbor[count] - 1;
                                /* mind the -1 , fortran again */
        xi = x[dim * i    ];
        yi = x[dim * i + 1];
        zi = x[dim * i + 2];

        xj = x[dim * j    ];
        yj = x[dim * j + 1];
        zj = x[dim * j + 2];

        r2 = (xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) 
           + (zi - zj) * (zi - zj);

        if (dim == 4) {
          wi = x[dim * i + 3];
          wj = x[dim * j + 3];
          r2 += (wi - wj) * (wi - wj);
        }

        dij1i = 1. / sqrt(r2);
        rij = r2 * dij1i;
        tmpaij = P0[i] - rij * 0.5 - (P0[i] * P0[i] -
                 P0[j] * P0[j]) * 0.5 * dij1i;

        Aij = PIx2 * P0[i] * tmpaij;

        dAijddij = PIx1 * P0[i] * (dij1i * dij1i *
                   (P0[i] * P0[i] - P0[j] * P0[j]) - 1.0);

        dAijddijdxj = dAijddij * (xj - xi) * dij1i;
        dAijddijdyj = dAijddij * (yj - yi) * dij1i;
        dAijddijdzj = dAijddij * (zj - zi) * dij1i;

        if (dim == 4) {
          dAijddijdwj = dAijddij * (wj - wi) * dij1i;
        }

        sumAij = sumAij + Aij;

        count2 = icount;
        sumAjk2 = 0.0;
        sumdAjkddjkdxj =  sumdAjkddjkdyj = sumdAjkddjkdzj = 0.0;

        if (dim == 4) {
          sumdAjkddjkdwj = 0.0;
        }

        p3p4Aij = -surften * (P3[i] + P4[i] * Aij);

        L80:k = ineighbor[count2] - 1;
                                /*same as above, -1 ! */
        if (j == k) goto L85;

        xk = x[dim * k    ];
        yk = x[dim * k + 1];
        zk = x[dim * k + 2];

        rjk2 = (xj - xk) * (xj - xk) + (yj - yk) * (yj - yk)
            + (zj - zk) * (zj - zk);

        if (dim == 4) {
          wk = x[dim * k + 3];
          rjk2 += (wj - wk) * (wj - wk);
        }

        djk1i = 1.0 / sqrt(rjk2);

        rjk = rjk2 * djk1i;

        if (P0[j] + P0[k] > rjk) {
          vdw2dif = P0[j] * P0[j] - P0[k] * P0[k];
          tmpajk = 2.0 * P0[j] - rjk - vdw2dif * djk1i;

          Ajk = PIx1 * P0[j] * tmpajk;

          sumAjk = sumAjk + Ajk;
          sumAjk2 = sumAjk2 + Ajk;

          dAjkddjk = PIx1 * P0[j] * djk1i * (djk1i * djk1i * vdw2dif - 1.0);

          dAjkddjkdxj = dAjkddjk * (xj - xk);
          dAjkddjkdyj = dAjkddjk * (yj - yk);
          dAjkddjkdzj = dAjkddjk * (zj - zk);

          f[dim * k    ] = f[dim * k    ] + dAjkddjkdxj * p3p4Aij;
          f[dim * k + 1] = f[dim * k + 1] + dAjkddjkdyj * p3p4Aij;
          f[dim * k + 2] = f[dim * k + 2] + dAjkddjkdzj * p3p4Aij;

          sumdAjkddjkdxj = sumdAjkddjkdxj + dAjkddjkdxj;
          sumdAjkddjkdyj = sumdAjkddjkdyj + dAjkddjkdyj;
          sumdAjkddjkdzj = sumdAjkddjkdzj + dAjkddjkdzj;

          if (dim == 4) {
            dAjkddjkdwj = dAjkddjk * (wj - wk);
            f[dim * k + 3] = f[dim * k + 3] + dAjkddjkdwj * p3p4Aij;
            sumdAjkddjkdwj = sumdAjkddjkdwj + dAjkddjkdwj;
          }
        }

        L85:count2 = count2 + 1;
        if (ineighbor[count2] != 0) {
          goto L80;
        } else {
          count2 = icount;
        }

        sumAijAjk = sumAijAjk + Aij * sumAjk2;

        sumdAijddijdxi = sumdAijddijdxi - dAijddijdxj;
        sumdAijddijdyi = sumdAijddijdyi - dAijddijdyj;
        sumdAijddijdzi = sumdAijddijdzi - dAijddijdzj;

        sumdAijddijdxiAjk = sumdAijddijdxiAjk - dAijddijdxj * sumAjk2;
        sumdAijddijdyiAjk = sumdAijddijdyiAjk - dAijddijdyj * sumAjk2;
        sumdAijddijdziAjk = sumdAijddijdziAjk - dAijddijdzj * sumAjk2;

        lastxj = dAijddijdxj * sumAjk2 + Aij * sumdAjkddjkdxj;
        lastyj = dAijddijdyj * sumAjk2 + Aij * sumdAjkddjkdyj;
        lastzj = dAijddijdzj * sumAjk2 + Aij * sumdAjkddjkdzj;

        dAidxj = surften * (P2[i] * dAijddijdxj +
                   P3[i] * sumdAjkddjkdxj + P4[i] * lastxj);
        dAidyj = surften * (P2[i] * dAijddijdyj +
                   P3[i] * sumdAjkddjkdyj + P4[i] * lastyj);
        dAidzj = surften * (P2[i] * dAijddijdzj +
                   P3[i] * sumdAjkddjkdzj + P4[i] * lastzj);

        f[dim * j    ] = f[dim * j    ] + dAidxj;
        f[dim * j + 1] = f[dim * j + 1] + dAidyj;
        f[dim * j + 2] = f[dim * j + 2] + dAidzj;

        if (dim == 4) {
          sumdAijddijdwi = sumdAijddijdwi - dAijddijdwj;

          sumdAijddijdwiAjk = sumdAijddijdwiAjk - dAijddijdwj * sumAjk2;

          lastwj = dAijddijdwj * sumAjk2 + Aij * sumdAjkddjkdwj;

          dAidwj = surften * (P2[i] * dAijddijdwj +
                     P3[i] * sumdAjkddjkdwj + P4[i] * lastwj);

          f[dim * j + 3] = f[dim * j + 3] + dAidwj;
        }

        count = count + 1;
        if (ineighbor[count] != 0) {
          goto L70;
        } else {
          count = count + 1;
        }

        Ai = P1[i] * si + P2[i] * sumAij + P3[i] * sumAjk + P4[i] * sumAijAjk;

        dAidxi = surften * (P2[i] * sumdAijddijdxi + P4[i] * sumdAijddijdxiAjk);
        dAidyi = surften * (P2[i] * sumdAijddijdyi + P4[i] * sumdAijddijdyiAjk);
        dAidzi = surften * (P2[i] * sumdAijddijdzi + P4[i] * sumdAijddijdziAjk);

        f[dim * i    ] = f[dim * i    ] + dAidxi;
        f[dim * i + 1] = f[dim * i + 1] + dAidyi;
        f[dim * i + 2] = f[dim * i + 2] + dAidzi;

        if (dim == 4) {
          dAidwi = surften*(P2[i] * sumdAijddijdwi + P4[i] * sumdAijddijdwiAjk);

          f[dim * i + 3] = f[dim * i + 3] + dAidwi;
        }

        totsasa = totsasa + Ai;
      }
    }
    /* fprintf( nabout,"SASA %f , ESURF %f \n",totsasa,totsasa*surften);  */
    esurf = totsasa * surften;
  }

  /* Compute the GB, Coulomb and Lennard-Jones energies and derivatives. */

  epol = elec = evdw = 0.0;

#pragma omp parallel reduction (+: epol, elec, evdw) \
  private (i, i34, ri, qi, qj, expmkf, dielfac, qi2h, qid2h, iaci, \
           xi, yi, zi, wi, k, j, j34, xij, yij, zij, wij, r2, qiqj, \
           rj, rb2, efac, fgbi, fgbk, temp4, temp6, eel, de, temp5, \
           rinv, r2inv, ic, r6inv, f6, f12, dedx, dedy, dedz, dedw, \
           iexw, xj, yj, zj, wj, threadnum, numthreads, foff, soff, \
           sumda, thi, ri1i, dij1i, datmp, daix, daiy, daiz, daiw, \
           dij2i, dij, sj, sj2, temp1, dij3i, tmpsd, dumbo, npairs)
  {

    /*
     * Get the thread number and the number of threads for multi-threaded
     * execution under OpenMP.  For all other cases, including ScaLAPACK,
     * MPI and single-threaded execution, use the values that have been
     * stored in mytaskid and numtasks, respectively.
     */

#if defined(OPENMP)
    threadnum = omp_get_thread_num();
    numthreads = omp_get_num_threads();
#else
    threadnum = mytaskid;
    numthreads = numtasks;
#endif

    /*
     * Compute offsets into the gradient and sumdeijda arrays for this
     * thread, but only if OPENMP is defined.
     */

#ifdef OPENMP
    soff = prm->Natom * threadnum;
    foff = dim * soff;
#else
    soff = 0;
    foff = 0;
#endif

    /*
     * Initialize the sumdeijda array inside of the parallel region.
     *
     * For OpenMP, the "first touch" memory allocation strategy will
     * locate the copy of the sumdeijda array that is initialized by
     * a particular CPU in memory that is local to that CPU.
     *
     * It is not necessary to synchronize OpenMP threads following
     * this loop because each thread initializes the particular copy
     * of the sumdeijda array that it subsequently updates.  A similar
     * argument applies for MPI tasks.
     */

    for (i = 0; i < prm->Natom; i++) {
      sumdeijda[soff + i] = 0.0;
    }

    /*
     * Allocate and initialize the iexw array used for skipping excluded
     * atoms.  Note that because of the manner in which iexw is used, it
     * is necessary to initialize it before only the first iteration of
     * the following loop.
     */

    iexw = ivector(-1, prm->Natom);
    for (i = -1; i < prm->Natom; i++) {
      iexw[i] = -1;
    }

    /*
     * Loop over all atoms i.
     *
     * Explicitly assign threads to loop indices for the following loop,
     * in a manner equivalent to (static, N) scheduling with OpenMP, and
     * identical to the manner in which threads are assigned in nblist.
     *
     * Synchronization of OpenMP threads will occur following this loop
     * because the parallel region ends after this loop.  Following
     * synchronization, a reduction of the sumdeijda array will be
     * performed.
     *
     * Synchronization of MPI tasks will occur via the MPI_Reduce function.
     */

    for (i = 0; i < prm->Natom; i++) {

#if defined(OPENMP) || defined(MPI) || defined(SCALAPACK)
      if ( !myroc(i, BLOCKSIZE, numthreads, threadnum) ) continue;
#endif

      ri = reff[i];
      qi = q[i];

      /*
       * If atom i is not frozen, compute the "diagonal" energy that
       * is a function of only the effective radius Ri but not of the
       * interatomic distance Dij.  Compute also the contribution of
       * the diagonal energy term to the sum by which the derivative
       * of Ri will be multiplied.
       */

      if (!frozen[i]) {
        expmkf = exp(-KSCALE * (*kappa) * ri) / (*diel_ext);
        dielfac = 1.0 - expmkf;
        qi2h = 0.5 * qi * qi;
        qid2h = qi2h * dielfac;
        epol += -qid2h / ri;

        sumdeijda[soff + i] +=
          qid2h - KSCALE * (*kappa) * qi2h * expmkf *ri;
      }

      /*
       * Skip the pair calculations if there are no atoms j on the
       * pair list of atom i.
       */

      npairs = upears[i];
      if ( !npairs ) continue;

      i34 = dim * i;

      xi = x[i34    ];
      yi = x[i34 + 1];
      zi = x[i34 + 2];

      if (dim == 4) {
        wi = x[i34 + 3];
      }

      iaci = prm->Ntypes * (prm->Iac[i] - 1);

      /*
       * Expand the excluded atom list into the iexw array by storing i
       * at array address j.
       */

      for (j = 0; j < prm->Iblo[i]; j++) {
        iexw[IexclAt[i][j] - 1] = i;
      }

      /* Initialize the derivative accumulators. */

      daix = daiy = daiz = daiw = 0.0;

      /* Select atoms j from the pair list.  Non-graceful error handling. */

      for (k = lpears[i]; k < lpears[i] + npairs; k++) {

        if (pearlist[i] == NULL) {
          fprintf( nabout,"NULL pair list entry in egb loop 3, taskid = %d\n", mytaskid);
          fflush(nabout);
        }
        j = pearlist[i][k];

        j34 = dim * j;

        /* Continue computing the non-diagonal energy term. */

        xij = xi - x[j34    ];
        yij = yi - x[j34 + 1];
        zij = zi - x[j34 + 2];
        r2 = xij * xij + yij * yij + zij * zij;

        if (dim == 4) {
          wij = wi - x[j34 + 3];
          r2 += wij * wij;
        }

        /*
         * Because index j is retrieved from the pairlist array it is
         * not constrained to a particular range of values; therefore,
         * the threads that have loaded the reff array must be
         * synchronized prior to the use of reff below.
         */

        qiqj = qi * q[j];
        rj = reff[j];
        rb2 = ri * rj;
        efac = exp(-r2 / (4.0 * rb2));
        fgbi = 1.0 / sqrt(r2 + rb2 * efac);
        fgbk = -(*kappa) * KSCALE / fgbi;

        expmkf = exp(fgbk) / (*diel_ext);
        dielfac = 1.0 - expmkf;

        epol += -qiqj * dielfac * fgbi;

        temp4 = fgbi * fgbi * fgbi;
        temp6 = qiqj * temp4 * (dielfac + fgbk * expmkf);
        de = temp6 * (1.0 - 0.25 * efac);

        temp5 = 0.5 * efac * temp6 * (rb2 + 0.25 * r2);

        /*
         * Compute the contribution of the non-diagonal energy term to the
         * sum by which the derivatives of Ri and Rj will be multiplied.
         */

        sumdeijda[soff + i] += ri * temp5;
        sumdeijda[soff + j] += rj * temp5;

        /*
         * Compute the Van der Waals and Coulombic energies for only
         * those pairs that are not on the excluded atom list.  Any
         * pair on the excluded atom list will have atom i stored at
         * address j of the iexw array.  It is not necessary to reset
         * the elements of the iexw array to -1 between successive
         * iterations in i because an i,j pair is uniquely identified
         * by atom i stored at array address j.  Thus for example, the
         * i+1,j pair would be stored at the same address as the i,j
         * pair but after the i,j pair were used.
         */

        if (iexw[j] != i) {

          rinv = 1. / sqrt(r2);
          r2inv = rinv * rinv;

          /*  gas-phase Coulomb energy:  */

          eel = qiqj * rinv;
          elec += eel;
          de -= eel * r2inv;

          /* Lennard-Jones energy:   */

          ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
          if( ic >= 0 ){
            r6inv = r2inv * r2inv * r2inv;
            f6 = prm->Cn2[ic] * r6inv;
            f12 = prm->Cn1[ic] * r6inv * r6inv;
            evdw += f12 - f6;
            de -= (12. * f12 - 6. * f6) * r2inv;
          }
        }

        /*
         * Sum to the gradient vector the derivatives of Dij that are
         * computed relative to the cartesian coordinates of atoms i and j.
         */

        dedx = de * xij;
        dedy = de * yij;
        dedz = de * zij;

        daix += dedx;
        daiy += dedy;
        daiz += dedz;

        f[foff + j34    ] -= dedx;
        f[foff + j34 + 1] -= dedy;
        f[foff + j34 + 2] -= dedz;

        if (dim == 4) {
          dedw = de * wij;
          daiw += dedw;
          f[foff + j34 + 3] -= dedw;
        }
      }

      /* Update the i elements of the gradient and the sumdeijda array. */

      f[foff + i34    ] += daix;
      f[foff + i34 + 1] += daiy;
      f[foff + i34 + 2] += daiz;

      if (dim == 4) {
        f[foff + i34 + 3] += daiw;
      }
    }

    /* Free the iexw array within this potentially parallel region of code. */

    free_ivector(iexw, -1, prm->Natom);
  }

    /* Perform a reduction of sumdeijda if MPI or SCALAPACK is defined. */

#if defined(MPI) || defined(SCALAPACK)
    ierror = MPI_Allreduce(sumdeijda, reductarr, prm->Natom,
                           MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (ierror != MPI_SUCCESS) {
        fprintf( nabout,"Error in egb sumdeijda reduction, error = %d  mytaskid = %d\n",
               ierror, mytaskid);
      }
    for (i = 0; i < prm->Natom; i++) {
      sumdeijda[i] = reductarr[i];
    }
#endif

    /*
     * Perform a reduction of sumdeijda if OPENMP is defined.
     * In the following, the (j,i) loop nesting is more efficient
     * than (i,j) loop nesting.
     */

#ifdef OPENMP
#pragma omp parallel for private (i) schedule(static, 8)
    for (j = 0; j < prm->Natom; j++) {
      for (i = 1; i < maxthreads; i++) {
        sumdeijda[j] += sumdeijda[prm->Natom * i + j];
      }
    }
#endif

#pragma omp parallel \
  private (i, i34, ri, qi, qj, expmkf, dielfac, qi2h, qid2h, iaci, \
           xi, yi, zi, wi, k, j, j34, xij, yij, zij, wij, r2, qiqj, \
           rj, rb2, efac, fgbi, fgbk, temp4, temp6, eel, de, temp5, \
           rinv, r2inv, ic, r6inv, f6, f12, dedx, dedy, dedz, dedw, \
           iexw, xj, yj, zj, wj, threadnum, numthreads, foff, \
           sumda, thi, ri1i, dij1i, datmp, daix, daiy, daiz, daiw, \
           dij2i, dij, sj, sj2, temp1, dij3i, tmpsd, dumbo, npairs)
  {
    /*
     * Get the thread number and the number of threads for multi-threaded
     * execution under OpenMP.  For all other cases, including ScaLAPACK,
     * MPI and single-threaded execution, use the values that have been
     * stored in mytaskid and numtasks, respectively.
     */

#if defined(OPENMP)
    threadnum = omp_get_thread_num();
    numthreads = omp_get_num_threads();
#else
    threadnum = mytaskid;
    numthreads = numtasks;
#endif

    /*
     * Compute an offset into the gradient array for this thread,
     * but only if OPENMP is defined.  There is no need to
     * compute an offset into the sumdeijda array because all
     * copies of this array have been reduced into copy zero.
     */

#ifdef OPENMP
    foff = prm->Natom * dim * threadnum;
#else
    foff = 0;
#endif

    /*
     * Compute the derivatives of the effective radius Ri of atom i
     * with respect to the cartesian coordinates of each atom j.  Sum
     * all of these derivatives into the gradient vector.
     *
     * Loop over all atoms i.
     *
     * Synchronization of OpenMP threads will occur following this loop
     * because the parallel region ends after this loop.  A reduction
     * of the gradient array will occur in the mme34 function, either
     * for OpenMP or MPI.  This reduction will synchronize the MPI
     * tasks, so an explicit barrier is not necessary at the end of
     * this loop.
     *
     * Explicitly assign threads to loop indices for the following loop,
     * in a manner equivalent to (static, N) scheduling with OpenMP, and
     * identical to the manner in which threads are assigned in nblist.
     */

    for (i = 0; i < prm->Natom; i++) {

#if defined(OPENMP) || defined(MPI) || defined(SCALAPACK)
      if ( !myroc(i, BLOCKSIZE, numthreads, threadnum) ) continue;
#endif

      /*
       * Don't calculate derivatives of the effective radius of atom i
       * if atom i is frozen or if there are no pair atoms j associated
       * with atom i.
       */

      npairs = lpears[i] + upears[i];
      if ( frozen[i] || !npairs ) continue;

      i34 = dim * i;

      xi = x[i34    ];
      yi = x[i34 + 1];
      zi = x[i34 + 2];

      if (dim == 4) {
        wi = x[i34 + 3];
      }

      ri = rborn[i] - BOFFSET;
      ri1i = 1. / ri;

      sumda = sumdeijda[i];

      if ( gb>1 ) {

        ri = rborn[i] - BOFFSET;
        thi = tanh( (gbalpha - gbbeta*psi[i] + gbgamma*psi[i]*psi[i])*psi[i] );
        sumda *= (gbalpha -2.0*gbbeta*psi[i] + 3.0*gbgamma*psi[i]*psi[i])
          *(1.0 - thi*thi)*ri/rborn[i];
      }

      /* Initialize the derivative accumulators. */

      daix = daiy = daiz = daiw = 0.0;

      /* Select atom j from the pair list.  Non-graceful error handling. */

      for (k = 0; k < npairs; k++) {

        if (pearlist[i] == NULL) {
          fprintf( nabout,"NULL pair list entry in egb loop 4, taskid = %d\n", mytaskid);
          fflush(nabout);
        }
        j = pearlist[i][k];

        j34 = dim * j;

        xij = xi - x[j34    ];
        yij = yi - x[j34 + 1];
        zij = zi - x[j34 + 2];
        r2 = xij * xij + yij * yij + zij * zij;

        if (dim == 4) {
          wij = wi - x[j34 + 3];
          r2 += wij * wij;
        }

        if ( r2 > rgbmaxpsmax2 ) continue;

        dij1i = 1.0 / sqrt(r2);
        dij2i = dij1i * dij1i;
        dij = r2 * dij1i;
        sj = fs[j] * (rborn[j] - BOFFSET);
        sj2 = sj * sj;

        /*
         * The following are the numerator of the first derivatives of the
         * effective radius Ri with respect to the interatomic distance Dij.
         * They are derived from the equations from the Appendix of Schaefer
         * and Froemmel as well as from the Taylor series expansion for d>>s
         * by Andreas Svrcek-Seiler.  The smooth rgbmax idea is from Andreas
         * Svrcek-Seiler and Alexey Onufriev.  The complete derivative is
         * formed by multiplying the numerator by -Ri*Ri.  The factor of Ri*Ri
         * has been moved to the terms that are multiplied by the derivative.
         * The negation is deferred until later.  When the chain rule is used
         * to form the first derivatives of the effective radius with respect
         * to the cartesian coordinates, an additional factor of Dij appears
         * in the denominator.  That factor is included in the following
         * expressions.
         */

        if ( dij > rgbmax + sj ) continue;

        if( dij > rgbmax - sj ){
          temp1 = 1./(dij-sj);
          dij3i = dij1i * dij2i;
          datmp = 0.125 * dij3i * ((r2 + sj2) *
                                   (temp1*temp1 - rgbmax2i) -
                                   2.0 * log(rgbmax*temp1));

        } else if(dij >4.0*sj) {
          tmpsd = sj2*dij2i;
          dumbo = TE + tmpsd*(TF+tmpsd*(TG+tmpsd*(TH+tmpsd*THH)));
          datmp = tmpsd*sj*dij2i*dij2i*dumbo;

        } else if (dij > ri + sj) {

          temp1 = 1. / (r2 - sj2);
          datmp = temp1 * sj * (-0.5 * dij2i + temp1)
            + 0.25 * dij1i * dij2i * log((dij - sj) / (dij + sj));

        } else if (dij > fabs(ri - sj)) {

          temp1 = 1. / (dij + sj);
          dij3i = dij1i * dij1i * dij1i;
          datmp =
            -0.25 * (-0.5 * (r2 - ri * ri + sj2) * dij3i * ri1i *
                     ri1i + dij1i * temp1 * (temp1 - dij1i)
                     - dij3i * log(ri * temp1));

        } else if (ri < sj) {

          temp1 = 1. / (r2 - sj2);
          datmp =
            -0.5 * (sj * dij2i * temp1 - 2. * sj * temp1 * temp1 -
                    0.5 * dij2i * dij1i * log((sj - dij) /
                                              (sj + dij)));

        } else {
          datmp = 0.;
        }

        /* Sum the derivatives into daix, daiy and daiz. */

        daix += xij * datmp;
        daiy += yij * datmp;
        daiz += zij * datmp;

        /*
         * Sum the derivatives relative to atom j (weighted by -sumdeijda[i])
         * into the gradient vector.  For example, f[j34 + 2] contains the
         * derivatives of Ri with respect to the z-coordinate of atom j.
         */

        datmp *= sumda;
        f[foff + j34    ] += xij * datmp;
        f[foff + j34 + 1] += yij * datmp;
        f[foff + j34 + 2] += zij * datmp;

        if (dim == 4) {
          daiw += wij * datmp;
          f[foff + j34 + 3] += wij * datmp;
        }
      }

      /*
       * Update the gradient vector with the sums of derivatives of the
       * effective radius Ri with respect to the cartesian coordinates.
       * For example, f[i34 + 1] contains the sum of derivatives of Ri
       * with respect to the y-coordinate of each atom.  Multiply by
       * -sumdeijda[i] here (instead of merely using datmp multiplied by
       * -sumdeijda) in order to distribute the product across the sum of
       * derivatives in an attempt to obtain greater numeric stability.
       */

      f[foff + i34    ] -= sumda * daix;
      f[foff + i34 + 1] -= sumda * daiy;
      f[foff + i34 + 2] -= sumda * daiz;

      if (dim == 4) {
        f[foff + i34 +3] -= sumda * daiw;
      }
    }
  }

  /* Return elec and evdw through the reference parameters eelt and enb. */

  *eelt = elec;
  *enb = evdw;

  /* Free the static arrays if STATIC_ARRAYS is undefined. */

#ifndef STATIC_ARRAYS
  if (reff != NULL) free_vector(reff, 0, n);
  reff = NULL;
  if (sumdeijda != NULL) free_vector(sumdeijda, 0, maxthreads*n);
  sumdeijda = NULL;
  if (psi != NULL) free_vector(psi, 0, n);
  psi = NULL;
#if defined(MPI) || defined(SCALAPACK)
  if (reductarr != NULL) free_vector(reductarr, 0, n);
  reductarr = NULL;
#endif
#endif

  return (epol);
}

/***********************************************************************
                            MME34()
************************************************************************/

/*
 * Here is the mme function for 3D or 4D, depending upon the dim variable.
 *
 * Calling parameters are as follows:
 *
 * x - input: the atomic (x,y,z) coordinates
 * f - updated: the gradient vector
 * iter - the iteration counter, which if negative selects the following:
 *        -1 print some energy values
 *        -3 call egb to deallocate static arrays, then deallocate grad
 *        -(any other value) normal execution
 */

static
REAL_T        mme34( REAL_T *x,  REAL_T *f, int *iter )

{
  extern REAL_T tconjgrad;

  REAL_T ebh, eba, eth, eta, eph, epa, enb, eel, enb14, eel14, ecn, e_gb, frms;
  REAL_T ene[12];
  REAL_T t1, t2;
  static REAL_T *grad=NULL;
  int i, j, k, goff, threadnum, numthreads, maxthreads;
  int dummy = 0;
  size_t n;

#if defined(MPI) || defined(SCALAPACK)
  int ierror;
  REAL_T reductarr[12];
#endif

  t1 = seconds();
  n = (size_t)prm->Natom;

  /*
   * If OPENMP is defined, set maxthreads to the maximum number of
   * OpenMP threads; otherwise, set maxthreads to one.
   */

#ifdef OPENMP
  maxthreads = omp_get_max_threads();
#else
  maxthreads = 1;
#endif

  /*
   * If the iteration count equals -3, call egb to deallocate the
   * static arrays, deallocate the gradient array, then return;
   * otherwise, simply return.
   */

  if (*iter == -3) {
    egb(lpairs, upairs, pairlist, x, grad, prm->Fs, prm->Rborn,
            prm->Charges, &kappa, &epsext, &enb, &eel, 1);
    if (grad != NULL) free_vector(grad, 0, maxthreads * dim * n);
    grad = NULL;
    return (0.0);
  }

  /* If the iteration count equals 0, print the header for task 0 only. */

  if (*iter == 0 && mytaskid == 0) {
    fprintf( nabout, "    iter    Total    bad         vdW     elect.     "
                     "cons.   genBorn      frms\n");
    fflush(nabout);
  }

  /* If the iteration count equals 0, initialize the timing variables. */

  if (*iter == 0) {
    tnonb = tpair = tbond = tangl = tphi = tborn = tcons = tmme = 0.0;
    tconjgrad = tmd = 0.0;
  }

  /*
   * Write the checkpoint file every nchk iteration if the chknm
   * variable is non-NULL, but for task 0 only.
   */

  if (chknm != NULL && (*iter > 0 && *iter % nchk == 0) && mytaskid == 0) {
    checkpoint(chknm, prm->Natom, x, *iter);
  }

  /*
   * Build the pair list if it hasn't already been built; rebuild
   * it every nsnb iteration.
   */

  if (nb_pairs < 0 || (*iter > 0 && *iter % nsnb == 0)) {
    nb_pairs = nblist(lpairs, upairs, pairlist, x, dummy, 1);
    t2 = seconds();
    tpair += t2 - t1;
    t1 = t2;
  }

  /*
   * If OPENMP is defined, allocate a gradient vector for each thread,
   * and let each thread initialize its gradient vector so that the
   * "first touch" strategy will allocate local memory.  If OpenMP
   * is not defined, allocate one gradient vector.
   *
   * Note: the following allocations assume that the dimensionality
   * of the problem does not change during one invocation of NAB.
   * If, for example, mme34 were called with dim==3 and then with dim==4,
   * these allocations would not be repeated for the larger value
   * of n that would be necessitated by dim==4.
   */

  if (grad == NULL) {
    grad = vector(0, maxthreads * dim * n);
  }

#ifdef OPENMP
#pragma omp parallel private (i, goff)
  {
    goff = dim * n * omp_get_thread_num();
    for (i = 0; i < dim * prm->Natom; i++) {
      grad[goff + i] = 0.0;
    }
  }
#else
  for (i = 0; i < dim * prm->Natom; i++) {
    grad[i] = 0.0;
  }
#endif

  t2 = seconds();
  tmme += t2 - t1;
  t1 = t2;

  ebh = ebond(prm->Nbonh, prm->BondHAt1, prm->BondHAt2,
              prm->BondHNum, prm->Rk, prm->Req, x, grad);
  eba = ebond(prm->Mbona, prm->BondAt1, prm->BondAt2,
              prm->BondNum, prm->Rk, prm->Req, x, grad);
  ene[3] = ebh + eba;
  t2 = seconds();
  tbond += t2 - t1;
  t1 = t2;

  eth = eangl(prm->Ntheth, prm->AngleHAt1, prm->AngleHAt2,
              prm->AngleHAt3, prm->AngleHNum, prm->Tk, prm->Teq, x, grad);
  eta = eangl(prm->Ntheta, prm->AngleAt1, prm->AngleAt2,
              prm->AngleAt3, prm->AngleNum, prm->Tk, prm->Teq, x, grad);
  ene[4] = eth + eta;
  t2 = seconds();
  tangl += t2 - t1;
  t1 = t2;

  eph = ephi(prm->Nphih, prm->DihHAt1, prm->DihHAt2,
             prm->DihHAt3, prm->DihHAt4, prm->DihHNum,
             prm->Pk, prm->Pn, prm->Phase, x, grad);
  epa = ephi(prm->Mphia, prm->DihAt1, prm->DihAt2,
             prm->DihAt3, prm->DihAt4, prm->DihNum,
             prm->Pk, prm->Pn, prm->Phase, x, grad);
  ene[5] = eph + epa;
  ene[6] = 0.0;       /*  hbond term not in Amber-94 force field */
  t2 = seconds();
  tphi += t2 - t1;
  t1 = t2;

  /* In the following lpairs is a dummy argument that is not used. */

  nbond(lpairs, prm->N14pairs, N14pearlist, 1, x, grad, &enb14, &eel14,
        scnb, scee);
  ene[7] = enb14;
  ene[8] = eel14;
  t2 = seconds();
  tnonb += t2 - t1;
  t1 = t2;

  if (e_debug) {
    EXPR("%9.3f", enb14);
    EXPR("%9.3f", eel14);
  }

  if (nconstrained) {
    ecn = econs(x, grad);
    t2 = seconds();
    tcons += t2 - t1;
    t1 = t2;
  } else
    ecn = 0.0;
  ene[9] = ecn;

  if (gb) {
    e_gb = egb(lpairs, upairs, pairlist, x, grad, prm->Fs, prm->Rborn,
               prm->Charges, &kappa, &epsext, &enb, &eel, 0);
    t2 = seconds();
    tborn += t2 - t1;
    t1 = t2;
    ene[1] = enb;
    ene[2] = eel;
    ene[10] = e_gb;
    ene[11] = esurf;
    if (e_debug) {
      EXPR("%9.3f", enb);
      EXPR("%9.3f", eel);
      EXPR("%9.3f", e_gb);
      EXPR("%9.3f", esurf);
    }
  } else {
    nbond(lpairs, upairs, pairlist, 0, x, grad, &enb, &eel, 1.0, 1.0);
    t2 = seconds();
    tnonb += t2 - t1;
    t1 = t2;
    ene[1] = enb;
    ene[2] = eel;
    ene[10] = 0.0;
    ene[11] = 0.0;
    if (e_debug) {
      EXPR("%9.3f", enb);
      EXPR("%9.3f", eel);
    }
  }

  /*
   * Perform a reduction over the gradient vector if OPENMP or MPI
   * is defined.
   *
   * If OPENMP is defined, the reduction is performed by the following loop.
   * The (j,i) loop nesting is more efficient than (i,j) nesting.
   *
   * If MPI is defined, the reduction is performed by MPI_Reduce.
   */

#if defined(OPENMP)

    for (i = 0; i < dim * prm->Natom; i++) {
      f[i] = grad[i];
    }
  
    goff = dim * prm->Natom;

#pragma omp parallel for private (i, j)

    for (j = 0; j < dim * prm->Natom; j++) {
      for (i = 1; i < maxthreads; i++) {
        f[j] += grad[goff * i + j];
      }
    }

#elif defined(MPI) || defined(SCALAPACK)

    ierror = MPI_Allreduce(grad, f, dim * prm->Natom,
                           MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (ierror != MPI_SUCCESS) {
      fprintf( nabout,"Error in mme34 grad reduction, error = %d  mytaskid = %d\n",
             ierror, mytaskid);
    }

#else

    for (i = 0; i < dim * prm->Natom; i++) {
      f[i] = grad[i];
    }

#endif

  ene[0] = 0.0;
  for (k = 1; k <= 11; k++)
    ene[0] += ene[k];

  for (k = 0; k < prm->Natom; k++) {        /* zero out frozen forces */
    if (frozen[k]) {
      f[dim * k + 0] = f[dim * k + 1] = f[dim * k + 2] = 0.0;

      if (dim == 4) {
        f[dim * k + 3] = 0.0;
      }
    }
  }

#ifdef PRINT_DERIV
  k = 0;
  for (i = 0; i < 105; i++) {
    k++;
    fprintf( nabout,"%10.5f", f[i]);
    if (k % 8 == 0)
      fprintf( nabout,"\n");
  }
  fprintf( nabout,"\n");
#endif

  frms = 0.0;
  for (i = 0; i < dim * prm->Natom; i++)
    frms += f[i] * f[i];
  frms = sqrt(frms / (dim * prm->Natom));

  /* If MPI is defined perform a reduction of the ene array. */

#if defined(MPI) || defined(SCALAPACK)
    ierror = MPI_Allreduce(ene, reductarr, 12, MPI_DOUBLE,
                           MPI_SUM, MPI_COMM_WORLD);
    if (ierror != MPI_SUCCESS) {
      fprintf( nabout,"Error in mme34 ene reduction, error = %d  mytaskid = %d\n",
             ierror, mytaskid);
    }
    for (i = 0; i < 12; i++) {
      ene[i] = reductarr[i];
    }
#endif

  /*
   * Print the energies and rms gradient but only for task zero,
   * and only for positive values of the iteration counter.
   */

  if (mytaskid == 0) {
    if (*iter > -1 && (*iter == 0 || *iter % ntpr == 0)) {
      fprintf( nabout,"ff:%4d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2e\n",
             *iter, ene[0], ene[3] + ene[4] + ene[5],
             ene[1] + ene[7], ene[2] + ene[8], ene[9]+ene[11], ene[10], frms);
      fflush(nabout);
    }
  }

  /* A value of -1 for the iteration counter is reserved for printing. */

  if( *iter == -1 ){
    fprintf( nabout, "     bond:  %15.9f\n", ene[3] );
    fprintf( nabout, "    angle:  %15.9f\n", ene[4] );
    fprintf( nabout, " dihedral:  %15.9f\n", ene[5] );
    fprintf( nabout, "    enb14:  %15.9f\n", ene[7] );
    fprintf( nabout, "    eel14:  %15.9f\n", ene[8] );
    fprintf( nabout, "      enb:  %15.9f\n", ene[1] );
    fprintf( nabout, "      eel:  %15.9f\n", ene[2] );
    fprintf( nabout, "      egb:  %15.9f\n", ene[10] );
    fprintf( nabout, "    econs:  %15.9f\n", ene[9] );
    fprintf( nabout, "    esurf:  %15.9f\n", ene[11] );
    fprintf( nabout, "    Total:  %15.9f\n", ene[0] );
  }

  /* If STATIC_ARRAYS is not defined, deallocate the gradient array. */

#ifndef STATIC_ARRAYS
  if (grad != NULL) free_vector(grad, 0, maxthreads * dim * n);
  grad = NULL;
#endif

  t2 = seconds();
  tmme += t2 - t1;

  return (ene[0]);
}

/***********************************************************************
                            MME_TIMER()
************************************************************************/
 
/* Print a timing summary but only for task zero. */

int        mme_timer( void )
{
  /* Use the maximum time from all MPI tasks or SCALAPACK processes. */

#if defined(MPI) || defined(SCALAPACK)

  REAL_T timarr[10], reductarr[10];

  timarr[0] = tcons;
  timarr[1] = tbond;
  timarr[2] = tangl;
  timarr[3] = tphi;
  timarr[4] = tpair;
  timarr[5] = tnonb;
  timarr[6] = tborn;
  timarr[7] = tmme;
  timarr[8] = tconjgrad;
  timarr[9] = tmd;

  MPI_Allreduce(timarr, reductarr, 10, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  tcons = reductarr[0];
  tbond = reductarr[1];
  tangl = reductarr[2];
  tphi = reductarr[3];
  tpair = reductarr[4];
  tnonb = reductarr[5];
  tborn = reductarr[6];
  tmme = reductarr[7];
  tconjgrad = reductarr[8];
  tmd = reductarr[9];

#endif

  if (mytaskid == 0) {
    fprintf( nabout,"\nFirst derivative timing summary:\n");
    fprintf( nabout,"   constraints %10.2f\n", tcons);
    fprintf( nabout,"   bonds       %10.2f\n", tbond);
    fprintf( nabout,"   angles      %10.2f\n", tangl);
    fprintf( nabout,"   torsions    %10.2f\n", tphi);
    fprintf( nabout,"   pairlist    %10.2f\n", tpair);
    fprintf( nabout,"   nonbonds    %10.2f\n", tnonb);
    fprintf( nabout,"   gen. Born   %10.2f\n", tborn);
    fprintf( nabout,"   mme         %10.2f\n", tmme);
    fprintf( nabout,"   Total       %10.2f\n\n",
           tcons + tbond + tangl + tphi + tpair + tnonb + tborn + tmme);
    fprintf( nabout,"   conj. grad. %10.2f\n", tconjgrad);
    fprintf( nabout,"   molec. dyn. %10.2f\n\n", tmd);
    fflush(nabout);
  }

  return (0);
}
