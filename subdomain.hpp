#pragma once
#include "mpiinfo.hpp"

// ============================================

struct DomainPair{
    int i,j;
};

// --------------------------------------------

class SubDomain {
public:
    std::vector<DomainPair> dplist;
    std::vector<DomainPair> dplist_reverse;
    void make_list(MPIinfo);
};

// ============================================