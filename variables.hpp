#pragma once
#include "mpiinfo.hpp"

// ========================================================================

struct Atom {
    int id;
    double  x,  y;
    double vx, vy;
};




struct Force {
    int id;
    double vx, vy;
};



class Variables {
public:
    std::vector<Atom> atoms;
    std::vector<std::vector<Atom>> other_atoms;
    std::vector<std::vector<int>> com_recv_list;
    std::vector<std::vector<int>> com_send_list;
    std::vector<std::vector<Force>> sending_force;
    
    double time;
    Variables(void) {time = 0.0;}
    void add_atoms(int id, double x, double y);
    void export_cdview(void);
    int number_of_atoms(void) {return static_cast<int>(atoms.size());}
    void set_initial_velocity(const double, MPIinfo);
    double margin_life;
    double max_velocity(void);

    double xp_max;
    double xp_min;
    double yp_max;
    double yp_min;

};

// ========================================================================