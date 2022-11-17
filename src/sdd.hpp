#pragma once
#include "systemparam.hpp"
#include "variables.hpp"
#include "mpiinfo.hpp"
#include "subregion.hpp"

// ============================================

struct Workload {
    int rank;
    unsigned long counts;
};

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
        unsigned long ideal_count;
        void set_limits(double, double, double, double);
        void migrate_atoms(std::vector<std::vector<Atom>>, Variables*, const MPIinfo &);
        void calc_bounds(Systemparam*, const MPIinfo &);
        void simple(Variables*, Systemparam*, const MPIinfo &);
        void global_sort(Variables*, Systemparam*, const MPIinfo &, SubRegion*);
        void voronoi_init(Variables*, Systemparam*, const MPIinfo &, SubRegion*);
        void voronoi(Variables*, Systemparam*, const MPIinfo &, SubRegion*, int, double, double);
        void voronoi_allocate(Variables*, Systemparam*, const MPIinfo &, SubRegion*);
        void center_atom_distance(int, double &, int &, const Atom, SubRegion*, Systemparam*);
};

// ============================================