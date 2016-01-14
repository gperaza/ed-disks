#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include "structures.h"

double bisection_solve(double(*)(double), double, double);
double sample_vbn(double);
double calc_ts(int);
double tang_vel_sl(double);
double tang_vel_st(double);
double tang_spring_sl(double);
double tang_spring_st(double);
double normal_force(double);
double normal_force_d(double);
double et_sl(double);
double calc_en(double, double);
double calc_et(double, double);
void calc_en_et(double, double, double, double);

extern particle *disk, *wall;
extern gsl_rng *rgen;
extern long colN;
extern double Gamma, wb, mu, kn, kt, dampR, g;
double en, et, gammav, beta, w0, wn, wt, tCol;
double gn0, gt0, t0, tSpr0;

void collision_rule_disk_disk(long i, long j) {
    double Ri = disk[i].R;
    double Rj = disk[j].R;
    double Ii = disk[i].iMoment;
    double Ij = disk[j].iMoment;
    double mi = disk[i].mass;
    double mj = disk[j].mass;

    double rxij = disk[i].x  - disk[j].x;
    double ryij = disk[i].y  - disk[j].y;
    double vxij = disk[i].vx - disk[j].vx;
    double vyij = disk[i].vy - disk[j].vy;
    double Rij =  Ri + Rj;
    double meff = mi*mj/(mi + mj);
    double alpha = 1.0/(1/meff + Ri*Ri/Ii + Rj*Rj/Ij);


    /* Normalized normal. */
    double rxijn = rxij/Rij;
    double ryijn = ryij/Rij;
    /*The tangent versor is (ry12n, -rx12n).*/

    double vn = vxij*rxijn + vyij*ryijn;
    double vt = vxij*ryijn - vyij*rxijn + Ri*disk[i].w + Rj*disk[j].w;
    assert(vn < 0);

    double vnx = vn*rxijn;
    double vny = vn*ryijn;
    double vtx = vt*ryijn;
    double vty = -vt*rxijn;

    calc_en_et(meff, alpha, vn, vt);
    //printf("en: %e et: %e\n", en, et);

    double denom = 1 + meff*(Ri*Ri/Ii + Rj*Rj/Ij);
    double An = meff*(1 + en);
    double At = meff*(1 - et)/denom;

    disk[i].vx -= (An*vnx + At*vtx)/mi;
    disk[i].vy -= (An*vny + At*vty)/mi;
    disk[i].w  -= (Ri/Ii)*At*vt;

    disk[j].vx += (An*vnx + At*vtx)/mj;
    disk[j].vy += (An*vny + At*vty)/mj;
    disk[j].w  -= (Rj/Ij)*At*vt;
}

double bisection_solve(double(*f)(double), double x1, double x2) {
    double xm, f1, f2, fm;

    f1 = f(x1); f2 = f(x2);
    assert(f1*f2 < 0);
    assert(fabs(f1) > 0);
    assert(fabs(f2) > 0);
    xm = (x1 + x2)*0.5;

    while (x1 < xm && x2 > xm) {
        fm = f(xm);
        if (fm == 0) return xm;

        if (f2*fm < 0) {
            x1 = xm;
            f1 = fm;
        }
        else {
            x2 = xm;
            f2 = fm;
        }
        xm = (x1 + x2)*0.5;
    }
    return xm;
}

double sample_vbn(double vn) {
    assert(vn < 0);
    double vbn, t, x, t1, t2, x1, x2;

    double A = Gamma*g/(wb*wb);

    double f_t_lims(double t) {
        return A*sin(wb*t) - vn*t - A*sin(wb*t1) + vn*t1;
    }
    double f_t(double t) {
        return A*sin(wb*t) - vn*(t - x);
    }

    if (fabs(vn) > A*wb) {
        t1 = t2 = M_PI/wb;
    } else {
        t1 = acos(vn/(A*wb))/wb;
        t2 = bisection_solve(&f_t_lims, M_PI/wb, 2.5*M_PI/wb);
    }
    t2 = t2 - 2*M_PI/wb;
    x1 = t1 - A*sin(wb*t1)/vn;
    x2 = t2 - A*sin(wb*t2)/vn;

    do {
        x = x2 + gsl_rng_uniform(rgen)*(x1 - x2);
        t = bisection_solve(&f_t, t2, t1);
        vbn = A*wb*cos(wb*t);
        assert(vn - vbn < 0);
    }  while (vbn <= en/(1 + en)*vn);

    return vbn;
}

void collision_rule_disk_wall(long i, long j) {
    /* i: disk j:wall */
    double nx = wall[j].nx;
    double ny = wall[j].ny;
    double xj = wall[j].x;
    double yj = wall[j].y;
    double xi = disk[i].x;
    double yi = disk[i].y;
    double vxi = disk[i].vx;
    double vyi = disk[i].vy;
    double Ri = disk[i].R;
    double Ii = disk[i].iMoment;

    double rijn = (xi - xj)*nx + (yi - yj)*ny;
    if (rijn < 0) {
        nx = -nx;
        ny = -ny;
        rijn = -rijn;
    }

    double vni = vxi*nx + vyi*ny;
    double vti = vxi*ny - vyi*nx;
    if (vni > 0) printf("vni=%e at collision %ld\n", vni, colN);
    assert(vni < 0);
    /* Sample the wall velocities. Only bottom wall is vibrating. (j
       = 0). Only normal vibration. */
    double vnj = 0;
    double vtj = 0;
    if (j == 0) {
        vnj = sample_vbn(vni);
    } else {
        vnj = 0;
    }

    double vn = vni - vnj;
    double vt = (vti + Ri*disk[i].w) - vtj;

    double vnx = vn*nx;
    double vny = vn*ny;
    double vtx = vt*ny;
    double vty = -vt*nx;

    /* Wall is massive. */
    double mi = disk[i].mass;
    double meff = mi;
    double alpha = 1.0/(1/meff + Ri*Ri/Ii);

    calc_en_et(meff, alpha, vn, vt);
    /* printf("en: %lf et: %lf ratio: %lf\n", en, et, fabs(vt/vn)); */

    double denom = 1 + meff*(Ri*Ri/Ii);
    double An = meff*(1 + en);
    double At = meff*(1 - et)/denom;

    disk[i].vx -= (An*vnx + At*vtx)/mi;
    disk[i].vy -= (An*vny + At*vty)/mi;
    disk[i].w  -= (Ri/Ii)*At*vt;
}

void calc_en_et(double meff, double alpha, double vn, double vt) {
    en = calc_en(meff, alpha);
    et = calc_et(vn, vt);
}

double calc_en(double meff, double alpha){
    gammav = 2*dampR*sqrt(kn*meff);
    beta = gammav/(2.0*meff);
    w0 = sqrt(kn/meff);
    wn = sqrt(fabs(w0*w0 - beta*beta));
    /* wt = sqrt(3*kt/meff); */
    wt = sqrt(kt/alpha);

    /* Calculate collision time. */
     if (beta <= w0/sqrt(2)) {
        tCol = 1/wn*(M_PI - atan(2*beta*wn/(wn*wn - beta*beta)));
    } else if (beta <= w0) {
        tCol = -1/wn*(atan(2*beta*wn/(wn*wn - beta*beta)));
    } else {
        tCol = 1/wn*log((beta + wn)/(beta - wn));
    }

     return exp(-beta*tCol);
}

double calc_et(double gn00, double gt00) {
    double gt = 0, tSpr, ts = 0;

    assert(gn00 < 0);
    double ratio = fabs(gt00/gn00);
    double slLimit = 6*mu - mu*kn/kt;
    /* Head on collision: et = 1. */
    if (ratio == 0) return 1;
    /* Tangential velocity beyond threshold, analytic solutions for
       full sliding. */
    else if (ratio >= slLimit) return 1 - 3*mu*(1 + en)/ratio;

    /* In between we need to solve numerically for et. */
    assert(gt00 != 0);
    gn0 = gn00;
    gt0 = gt00;
    t0 = 0;
    tSpr0 = 0;

    if (gammav == 0 && fabs(gt0/(double)gn0) > mu*kn/kt) {
        /* Start sliding. */
        ts = calc_ts(1);
        gt = tang_vel_sl(ts);
        tSpr = tang_spring_sl(normal_force(ts));
        t0 = ts;
        gt0 = gt;
        tSpr0 = tSpr;
    }
    /* Now we are sticking. */
    int sticking = 1;

    while (ts < tCol) {
        if (sticking) {
            ts = calc_ts(0);
            gt = tang_vel_st(ts);
            tSpr = tang_spring_st(ts);
            sticking = 0;
        } else {
            ts = calc_ts(1);
            gt = tang_vel_sl(ts);
            tSpr = tang_spring_sl(normal_force(ts));
            sticking = 1;
        }
        t0 = ts;
        gt0 = gt;
        tSpr0 = tSpr;
    }

    return gt/(double)gt00;
}

double flt(double t) {
    double sign;
    if (tSpr0 == 0) {
        sign = copysign(1, gt0);
    } else {
        sign = copysign(1, tSpr0);
    }
    return (tang_vel_sl(t) -
            sign*mu/kt*normal_force_d(t));
}

double ftl(double t) {
    return fabs(tang_spring_st(t)) -
        mu/kt*normal_force(t);
}

double calc_ts(int sl) {
    double ti = t0, ts;
    double(*f)(double);
    if (sl) f = &flt;
    else f = &ftl;
    while (fabs(f(ti)) < 1e-15) {
        ti += tCol/1000000.0;
    }
    double tf = ti;
    while (f(ti)*f(tf) > 0) {
        tf += tCol/100.0;
        if (tf > tCol) return tCol;
    }
    ts = bisection_solve(f, ti, tf);
    return ts;
}

double normal_force(double t) {
    double fn = -gn0*exp(-beta*t)/wn*(gammav*wn*cos(wn*t) +
                                      (kn - beta*gammav)*sin(wn*t));
    return fn;
}

double normal_force_d(double t) {
    double dfn = -gn0*exp(-beta*t)/wn*
        ((kn - 2*beta*gammav)*wn*cos(wn*t) -
         (kn*beta + (wn*wn - beta*beta)*gammav)*sin(wn*t));
    return dfn;
}

double tang_vel_st(double t) {
    return -wt*tSpr0*sin(wt*(t - t0)) + gt0*cos(wt*(t - t0));
}

double tang_vel_sl(double t) {
    double sign;

    if (tSpr0 == 0) {
        sign = copysign(1, gt0);
    } else {
        sign = copysign(1, tSpr0);
    }
    return gt0 -
        sign*3*mu*gn0*exp(-beta*(t0))*(exp(-beta*(t - t0))*
                                       (cos(wn*t) - beta/wn*sin(wn*t)) -
                                       (cos(wn*t0) - beta/wn*sin(wn*t0)));
}

double tang_spring_sl(double fn){
    if (tSpr0 == 0) {
        return copysign(1, gt0)*mu*fn/kt;
    } else {
        return copysign(1, tSpr0)*mu*fn/kt;
    }
}

double tang_spring_st(double t) {
    double ts = tSpr0*cos(wt*(t - t0)) + gt0/wt*sin(wt*(t - t0));
    return ts;
}
