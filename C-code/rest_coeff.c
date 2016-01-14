#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/*----------------------------------------------------------------------------*/
/* Functions*/
/*----------------------------------------------------------------------------*/
double calc_t_col();
double calc_en();
double calculate_et(double, double);
double calc_ts(int);
double tang_vel_sl(double);
double tang_vel_st(double);
double tang_spring_sl(double);
double tang_spring_st(double);
double normal_force(double);
double normal_force_d(double);
double bisection_solve(double (*)(double), double, double);
double et_sl(double);
/*----------------------------------------------------------------------------*/

/* Global variables*/
/*----------------------------------------------------------------------------*/
double mu = 0.1;
double kn = 4.5e6;
extern double enTrue;
extern double en;
double kt, gammav, w0, wn, beta, meff, tcol;

void set_constants(double meff_in) {
    meff = meff_in;
    kt = 0.773*kn;
    gammav = 2*0.035*sqrt(kn*meff);
    w0 = sqrt(kn/meff);
    beta = gammav/(2.0*meff);
    wn = sqrt(fabs(w0*w0 - beta*beta));
    tcol = calc_t_col();

    en = calc_en(tcol);
}

double calc_t_col() {
    if (beta <= w0/sqrt(2)) {
        return 1/wn*(M_PI - atan(2*beta*wn/(wn*wn - beta*beta)));
    } else if (beta <= w0) {
        return -1/wn*(atan(2*beta*wn/(wn*wn - beta*beta)));
    } else {
        return 1/wn*log((beta + wn)/(beta - wn));
    }
}

double calc_en() {
    return exp(-beta*tcol);
}

double gn0, gt0, t0, tSpr0, tc;
gsl_interp_accel *acc;
gsl_spline *spline;
double slLimit;

void prep_et() {
    /* Set ups the spline to sample et. */
    double et;
    double gn0 = -1, gt0;

    slLimit = 6*mu - mu*kn/kt;

    double *x = (double*)calloc(100, sizeof(double));
    double *y = (double*)calloc(100, sizeof(double));

    /* Write first point. */
    x[0] = 0; y[0] = cos(sqrt(3*kt/meff)*t_col());
    /* Calculate intermediate points. */
    long i;
    for (i = 1; i < 100; i++) {
        gt0 = slLimit/100.0*i;
        double ratio = fabs(gt0/(double)gn0);
        et = calculate_et(gn0, gt0);
        x[i] = ratio;
        y[i] = et;
    }
    /* Write last point. */
    x[100] = 6*mu - mu*kn/kt;
    y[100] = 1 - 3*mu*(1 + calc_en())/x[100];

    /* Interpolation with cubic spline. */
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, 101);
    gsl_spline_init(spline, x, y, 101);
}

double get_et(double ratio) {
    double et;

    if (ratio == 0) et = 1;
    else if (ratio >= slLimit) et = et_sl(ratio);
    else et = gsl_spline_eval(spline, ratio, acc);

    return et;
}

double calculate_et(double gn00, double gt00) {
    double gt = 0, tSpr, ts = 0;

    assert(gn00 < 0);
    assert(gt00 != 0);
    gn0 = gn00;
    gt0 = gt00;
    t0 = 0;
    tSpr0 = 0;
    tc = t_col();

    if (gammav == 0 && fabs(gt0/(double)gn0) > mu*kn/kt) {
        ts = calc_ts(1);
        gt = tang_vel_sl(ts);
        tSpr = tang_spring_sl(normal_force(ts));
        t0 = ts;
        gt0 = gt;
        tSpr0 = tSpr;
    }
    int sticking = 1;
    while (ts < tc) {
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

double normal_force(double t) {
    double w0 = sqrt(kn/meff);
    double beta = gammav/(2.0*meff);
    double wn = sqrt(fabs(w0*w0 - beta*beta));
    double fn = -gn0*exp(-beta*t)/wn*(gammav*wn*cos(wn*t) +
                                      (kn - beta*gammav)*sin(wn*t));
    return fn;
}

double normal_force_d(double t) {
    double w0 = sqrt(kn/meff);
    double beta = gammav/(2.0*meff);
    double wn = sqrt(fabs(w0*w0 - beta*beta));
    double dfn = -gn0*exp(-beta*t)/wn*
        ((kn - 2*beta*gammav)*wn*cos(wn*t) -
         (kn*beta + (wn*wn - beta*beta)*gammav)*sin(wn*t));
    return dfn;
}

double tang_vel_st(double t) {
    double wt = sqrt(3*kt/meff);
    return -wt*tSpr0*sin(wt*(t - t0)) + gt0*cos(wt*(t - t0));
}

double tang_vel_sl(double t) {
    double wn = sqrt(kn/meff);
    double beta = gammav/(2.0*meff);
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
    double wt = sqrt(3*kt/meff);
    double ts = tSpr0*cos(wt*(t - t0)) + gt0/wt*sin(wt*(t - t0));
    return ts;
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
        ti += tc/1000000.0;
    }
    double tf = ti;
    while (f(ti)*f(tf) > 0) {
        tf += tc/100.0;
        if (tf > tc) return tc;
    }
    ts = bisection_solve(f, ti, tf);
    return ts;
}

double et_sl(double ratio) {
    return 1 - 3*mu*(1 + enTrue)/ratio;
}

double bisection_solve(double(*f)(double), double a, double b) {
    double s;
    double fa = f(a), fb = f(b), fs;
    double ftol = 1e-16;
    assert(fabs(fa) >0);
    assert(fabs(fb) >0);
    assert(fa*fb < 0);
    do {
        s = 0.5*(b + a);
        fa = f(a); fb = f(b); fs = f(s);
        if (fa*fs < 0) b = s;
        else if (fb*fs < 0) a = s;
        else return s;
    } while(fabs(a - b) > ftol);
    return s;
}
