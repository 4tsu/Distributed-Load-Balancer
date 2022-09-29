#pragma once 
#include <vector>
#include "variables.hpp"
#include "systemparam.hpp"


// =================================================

class Observer {
public:
    void export_cdview(std::vector<Atom> atoms, Systemparam systemparam);
};

// =================================================