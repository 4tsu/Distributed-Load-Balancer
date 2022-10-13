#pragma once 
#include "variables.hpp"
#include "systemparam.hpp"
#include "mpiinfo.hpp"
#include "pairlist.hpp"

// =================================================

class Observer {
public:
    void export_cdview(std::vector<Atom> atoms, Systemparam sysp, MPIinfo mi);
    double kinetic_energy(Variables*, const MPIinfo &mi, Systemparam*);
    double potential_energy(Variables*, PairList*, const MPIinfo &mi, Systemparam*);
};

// =================================================