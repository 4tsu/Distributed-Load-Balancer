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


void periodic_distance(double &dx, double &dy);
void periodic_coordinate(double &x, double &y);

// ===================================================