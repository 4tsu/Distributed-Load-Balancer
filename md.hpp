#pragma once
#include "variables.hpp"
#include "obserber.hpp"

// =======================================

class MD {
private:
    Variables *vars;
    Observer *obs;
    void makeconf(void);
    void periodic(void);
    void update_position(void);
    void calculate_force(void);
    void calculate(void);
    void set_params(int STEPS, int OB_INTERVAL, double dt);
    int STEPS;
    int OB_INTERVAL;
    int dt;

public:
    MD(void);
    ~MD(void);
    void run(void);
};

// =======================================