#pragma once
#include "variables.hpp"
#include "observer.hpp"
#include "systemparam.hpp"
#include "mpiinfo.hpp"

// =======================================

class MD {
private:
    Variables *vars;
    Observer *obs;
    Systemparam *sysp;
    MPIinfo mi;
    void makeconf(void);
    void periodic(void);
    void update_position(void);
    void calculate_force(void);
    void calculate(void);
    int steps;
    int ob_interval;
    int dt;
    int sdd_type = 0;

public:
    MD(MPIinfo mi);
    ~MD(void);
    void run(void);
    void set_params(int steps, int ob_interval, double dt);
    void set_box(int N, double xl, double y, double cutoff);
    void set_margin(double margin);
    void set_sdd(int sdd_type);
};

// ---------------------------------------

// =======================================