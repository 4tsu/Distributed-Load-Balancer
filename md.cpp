#include <iostream>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include <cassert>
#include "md.hpp"
#include "systemparam.hpp"
#include "observer.hpp"
#include "variables.hpp"

// =====================================================

MD::MD(void){
    vars = new Variables();
    sysp = new Systemparam();
    obs = new Observer();
}

MD::~MD(void){
    delete vars;
    delete obs;
}

// -----------------------------------------------------

void MD::set_params(int STEPS, int OB_INTERVAL, double dt) {
    this->STEPS = STEPS;
    this->OB_INTERVAL = OB_INTERVAL;
    this->dt = dt;
}

void MD::set_box(int N, double xl, double yl, double cutoff) {
    sysp->set_params(N, xl, yl, cutoff);
    sysp->calc_params();
}

void MD::set_margin(double margin) {
    sysp->margin = margin;
}

void MD::set_sdd(int sdd_type) {
    this->sdd_type = sdd_type;
}

void MD::makeconf(void) {
    int N = sysp->N;
    double xl = sysp->xl;
    double yl = sysp->yl;
    double x_min = sysp->x_min;
    double y_min = sysp->y_min;
    double x_max = sysp->x_max;
    double y_max = sysp->y_max;
    int xppl = ceil(sqrt(xl*N/yl));
    int yppl = ceil(sqrt(yl*N/xl));
    double pitch = std::min(xl/xppl, yl/yppl);

    // 等間隔分割

    for (int i=0; i<N; i++) {
        int iy = static_cast<int>(i/xppl);
        int ix = i%xppl;
        double x = ix * pitch;
        double y = iy * pitch;

        // どのプロセスに分配するかを判断する

        x += x_min;
        y += y_min;
        vars->add_atoms(i,x,y);
        assert(x_min<=x && x<=x_max);
        assert(y_min<=y && y<=y_max);
    }
}

void MD::run(void) {
    makeconf();
    assert(sysp->N != 0);
    obs->export_cdview(vars->atoms, *sysp);
}
