#pragma once
#include "mpiinfo.hpp"

// ====================================================

class Systemparam {
public:
    double xl;
    double yl;
    unsigned long N = 0;
    double cutoff;
    double margin;
    double co_margin;
    
    double CL2;
    double RC2;
    double RC6;
    double RC12;
    double C0;

    double x_max;
    double x_min;
    double y_max;
    double y_min;

    void set_params(unsigned long N, double xl, double yl, double cutoff);
    void calc_params(void);
    void calc_margin(void);

};

void periodic_distance(double &dx, double &dy, Systemparam*);
void periodic_coordinate(double &x, double &y, Systemparam*);

// ====================================================