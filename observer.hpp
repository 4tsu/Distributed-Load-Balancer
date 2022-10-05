#pragma once 
#include "variables.hpp"
#include "systemparam.hpp"
#include "mpiinfo.hpp"

// =================================================

class Observer {
public:
    void export_cdview(std::vector<Atom> atoms, Systemparam sysp, MPIinfo mi);
};

// =================================================