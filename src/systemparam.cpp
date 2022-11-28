#include "systemparam.hpp"

// =======================================


namespace systemparam {
    double xl;
    double yl;
    double xlh;
    double ylh;
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



    void calc_params(void) {
        CL2 = (cutoff * cutoff);
        RC2 = 1.0 / CL2;
        RC6 = RC2 * RC2 * RC2;
        RC12 = RC6 * RC6;
        C0 = -4.0 * (RC12 - RC6);

        x_max =  xl/2;
        x_min = -xl/2;
        y_max =  yl/2;
        y_min = -yl/2;
        xlh = xl/2;
        ylh = yl/2;
    }



    void calc_margin(void) {
        co_margin = margin + cutoff;
    }

}


namespace sysp = systemparam;


void periodic_distance(double &dx, double &dy) {
    if (dx < -sysp::xlh) dx += sysp::xl;
    if (dx >  sysp::xlh) dx -= sysp::xl;
    if (dy < -sysp::ylh) dy += sysp::yl;
    if (dy >  sysp::ylh) dy -= sysp::yl;
}



void periodic_coordinate(double &x, double &y) {
    if (sysp::x_min > x)      x += sysp::xl;
    else if (sysp::x_max < x) x -= sysp::xl;
    if (sysp::y_min > y)      y += sysp::yl;
    else if (sysp::y_max < y) y -= sysp::yl;
}

// =======================================