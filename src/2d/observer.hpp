#pragma once 
#include "variables.hpp"
#include "systemparam.hpp"
#include "mpiinfo.hpp"
#include "pairlist.hpp"

// =================================================

class Observer {
private:
    void export_cdview_independent(Variables*, MPIinfo &);
    void concatenate_cdview(MPIinfo &mi, int);
    void checkpoint_independent(const int, Variables*, MPIinfo &);
    void concatenate_checkpoint(std::string, MPIinfo &);

public:
    void export_cdview(Variables* sysp, MPIinfo mi, int count_begin = 0);
    void export_checkpoint(std::string, const int, Variables* sysp, MPIinfo &);
    double kinetic_energy(Variables*);
    double potential_energy(Variables*, PairList*);
    void export_workload(const int, Variables*, MPIinfo &);
};

// =================================================

void export_three(const std::string, const int, const double, const double, const double);

// =================================================