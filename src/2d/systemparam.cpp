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

// =======================================