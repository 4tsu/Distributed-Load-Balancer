#pragma once
#include "mpiinfo.hpp"

// ========================================================================

struct Atom {
    int id;
    double  x,  y,  z;
    double vx, vy, vz;
};

class Variables {
public:
    std::vector<Atom> atoms;
    double time;
    Variables(void) {time = 0.0;}
    void add_atoms(int id, double x, double y);
    void export_cdview(void);
    int number_of_atoms(void) {return static_cast<int>(atoms.size());}
    void set_initial_velocity(const double, MPIinfo);

    double xp_max;
    double xp_min;
    double yp_max;
    double yp_min;

};

// ========================================================================