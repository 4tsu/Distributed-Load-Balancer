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
        void init(Variables*, const MPIinfo &, SubRegion*);
        void run(Variables*, const MPIinfo &, SubRegion*);
        unsigned long ideal(const MPIinfo &);
    
    private:
        int sdd_type = -1;
        double top, bottom, right, left;
        unsigned long ideal_count;
        std::vector<double> all_biases;
        void set_limits(double, double, double, double);
        void migrate_atoms(std::vector<std::vector<Atom>>, Variables*, const MPIinfo &);
        void calc_bounds(const MPIinfo &);
        void simple(Variables*, const MPIinfo &);
        void global_sort(Variables*, const MPIinfo &);
        void voronoi_init(Variables*, const MPIinfo &, SubRegion*);
        void voronoi(Variables*, const MPIinfo &, SubRegion*, int, double, double);
        void voronoi_allocate(Variables*, const MPIinfo &, SubRegion*);
        void center_atom_distance(int, double &, int &, const Atom, SubRegion*);
        void voronoi_figure(Variables*, const MPIinfo &);
};

// ============================================