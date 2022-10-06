#pragma once
#include "mpiinfo.hpp"

// ============================================

struct DomainPair{
    int i,j;
};

// --------------------------------------------

class DomainPairList {
public:
    std::vector<DomainPair> dplist;
    void make_list(MPIinfo);
};

// ============================================