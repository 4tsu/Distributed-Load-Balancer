#pragma once 
#include "variables.hpp"
#include "systemparam.hpp"
#include "mpiinfo.hpp"
#include "pairlist.hpp"

// =================================================

class Observer {
public:
    void export_cdview(std::vector<Atom> atoms, Systemparam sysp, MPIinfo mi);
    double kinetic_energy(Variables*, Systemparam*);
    double potential_energy(Variables*, PairList*, Systemparam*);
};

// =================================================

void export_three(const std::string, const int, const double, const double, const double);

// =================================================