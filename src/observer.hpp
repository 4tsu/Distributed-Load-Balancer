#pragma once 
#include "variables.hpp"
#include "systemparam.hpp"
#include "mpiinfo.hpp"
#include "pairlist.hpp"

// =================================================

class Observer {
public:
    void export_cdview(Variables*, Systemparam* sysp, MPIinfo mi);
    double kinetic_energy(Variables*, Systemparam*);
    double potential_energy(Variables*, PairList*, Systemparam*);
    void export_cdview_independent(Variables*, Systemparam*, MPIinfo &);
    void concatenate_cdview(int, int);
};

// =================================================

void export_three(const std::string, const int, const double, const double, const double);

// =================================================