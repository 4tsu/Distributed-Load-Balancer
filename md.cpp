#include <iostream>
#include <assert.h>
#include <math.h>
#include "md.hpp"
#include "systemparam.hpp"
#include "observer.hpp"
#include "variables.hpp"

// =====================================================

MD::MD(void){
    vars = new Variables();
    obs = new Observer();
}

MD::~MD(void){
    delete vars;
    delete obs;
}

// -----------------------------------------------------
void MD::set_params(int STEPS, int OB_INTERVAL, double, dt) {
    this.STEPS = STEPS
    this.OB_INTERVAL = OB_INTERVAL
    this.dt = dt
}

void MD::makeconf(void) {

}

void MD::run(void) {

    makeconf();
}
