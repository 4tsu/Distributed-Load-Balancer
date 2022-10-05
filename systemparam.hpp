#pragma once
#include "mpiinfo.hpp"
// ====================================================

class Systemparam {
public:
    double xl;
    double yl;
    int N = 0;
    double cutoff;
    double margin;
    int myN;
    
    double CL2;
    double RC2;
    double RC6;
    double RC12;
    double C0;

    double x_max;
    double x_min;
    double y_max;
    double y_min;

    void set_params(int N, double xl, double yl, double cutoff, int procs);
    void calc_params(void);

};

void adjust_periodic(double &dx, double &dy, Systemparam sysp);

// ====================================================