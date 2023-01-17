#pragma once

namespace systemparam {
    extern double xl;
    extern double yl;
    extern double xlh;
    extern double ylh;
    extern unsigned long N;
    extern double cutoff;
    extern double margin;
    extern double co_margin;
    
    extern double CL2;
    extern double RC2;
    extern double RC6;
    extern double RC12;
    extern double C0;

    extern double x_max;
    extern double x_min;
    extern double y_max;
    extern double y_min;

    void calc_params();
    void calc_margin();
} // namespace

// ====================================================

namespace sysp = systemparam;



inline
void periodic_distance(double &dx, double &dy) {
    if (dx < -sysp::xlh) dx += sysp::xl;
    else if (dx >  sysp::xlh) dx -= sysp::xl;
    if (dy < -sysp::ylh) dy += sysp::yl;
    else if (dy >  sysp::ylh) dy -= sysp::yl;
}



inline
void periodic_coordinate(double &x, double &y) {
    if (sysp::x_min > x)      x += sysp::xl;
    else if (sysp::x_max < x) x -= sysp::xl;
    if (sysp::y_min > y)      y += sysp::yl;
    else if (sysp::y_max < y) y -= sysp::yl;
}

// ===================================================