#pragma once

// ====================================================

class Systemparam {
public:
    double xl;
    double yl;
    double N;
    double cutoff;
    double margin;
    
    double CL2;
    double RC2;
    double RC6;
    double RC12;
    double C0;

    double x_max;
    double x_min;
    double y_max;
    double y_min;

    void set_params(int N, double xl, double yl, double cutoff);
    void calc_params(void);

};

void adjust_periodic(double &dx, double &dy, Systemparam sysp);

// ====================================================