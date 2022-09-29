#include "systemparam.hpp"

// =======================================

void Systemparam::set_params(int N, double xl, double yl, double cutoff) {
    this->N = N;
    this->xl = xl;
    this->yl = yl;
    this->cutoff = cutoff;
}

void Systemparam::calc_params(void) {
    double cutoff = this->cutoff;
    double CL2 = (cutoff * cutoff);
    double RC2 = 1.0 / CL2;
    double RC6 = RC2 * RC2 * RC2;
    double RC12 = RC6 * RC6;
    double C0 = -4.0 * (RC12 - RC6);
    this->CL2 = CL2;
    this->RC2 = RC2;
    this->RC6 = RC6;
    this->RC12 = RC12;
    this->C0 = C0;

    xl = this->xl;
    yl = this->yl;
    this->x_max =  xl/2;
    this->x_min = -xl/2;
    this->y_max =  yl/2;
    this->y_min = -yl/2;
}

void adjust_periodic(double &dx, double &dy, Systemparam sysp) {
    const double xl = sysp.xl;
    const double yl = sysp.yl;
    const double xlh = xl * 0.5;
    const double ylh = yl * 0.5;
    if (dx < -xlh) dx += xl;
    if (dx >  xlh) dx -= xl;
    if (dy < -ylh) dy += yl;
    if (dy >  ylh) dy -= yl;
}

// =======================================