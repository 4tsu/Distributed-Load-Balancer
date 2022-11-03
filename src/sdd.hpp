#pragma once
#include "systemparam.hpp"
#include "variables.hpp"
#include "mpiinfo.hpp"
#include "subregion.hpp"

// ============================================

class Sdd {
    public:
        Sdd(const int sdd_type);
        ~Sdd(void);
        void init(Variables*, Systemparam*, const MPIinfo &, SubRegion*);
        void run(Variables*, Systemparam*, const MPIinfo &, SubRegion*);
        unsigned long ideal(Systemparam*, const MPIinfo &);
    
    private:
        int sdd_type = -1;
        double top, bottom, right, left;
        void migrate_atoms(std::vector<std::vector<Atom>>, Variables*, const MPIinfo &);
        void calc_bounds(Systemparam*, const MPIinfo &);
        void simple(Variables*, Systemparam*, const MPIinfo &);
        void global_sort(Variables*, Systemparam*, const MPIinfo &, SubRegion*);
        void voronoi_init(Variables*, Systemparam*, const MPIinfo &, SubRegion*);
        void voronoi(Variables*, Systemparam*, const MPIinfo &, SubRegion*);
};

// ============================================