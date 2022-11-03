#pragma once
#include "variables.hpp"
#include "observer.hpp"
#include "systemparam.hpp"
#include "mpiinfo.hpp"
#include "pairlist.hpp"
#include "subregion.hpp"

// =======================================

class MD {
private:
    Variables *vars;
    Observer *obs;
    Systemparam *sysp;
    SubRegion *sr;
    PairList *pl;
    MPIinfo mi;
    void makeconf(void);
    void periodic(void);
    void update_position(double);
    void calculate_force(void);
    void calculate(void);
    void make_pair(void);
    void check_pairlist(void);
    int steps;
    int begin_step = 0;
    int ob_interval;
    double dt;
    int sdd_type = 0;
    void communicate_atoms(void);
    void communicate_force(void);
    std::string config;
    void read_data(std::string filename, Variables* vars, Systemparam* sysp, const MPIinfo &mi);

public:
    MD(MPIinfo mi);
    ~MD(void);
    void run(void);
    void set_params(int steps, int ob_interval, double dt);
    void set_box(unsigned long N, double xl, double y, double cutoff);
    void set_margin(double margin);
    void set_sdd(int sdd_type);
    void set_config(const std::string);
};

// ---------------------------------------



// =======================================