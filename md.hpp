#pragma once
#include "variables.hpp"
#include "observer.hpp"
#include "systemparam.hpp"

// =======================================

class MD {
private:
    Variables *vars;
    Observer *obs;
    Systemparam *sysp;
    void makeconf(void);
    void periodic(void);
    void update_position(void);
    void calculate_force(void);
    void calculate(void);
    int STEPS;
    int OB_INTERVAL;
    int dt;
    int sdd_type = 0;

public:
    MD(void);
    ~MD(void);
    void run(void);
    void set_params(int STEPS, int OB_INTERVAL, double dt);
    void set_box(int N, double xl, double y, double cutoff);
    void set_margin(double margin);
    void set_sdd(int sdd_type);
};

// =======================================