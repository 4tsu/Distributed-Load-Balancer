#include <iostream>
#include <fstream>
#include <random>
#include "systemparam.hpp"
#include "variables.hpp"

// =================================

void Variables::add_atoms(int id, double x, double y) {
    Atom a;
    a.id = id;
    a.x = x;
    a.y = y;
    a.vx = 0.0;
    a.vy = 0.0;
    atoms.push_back(a);
}


