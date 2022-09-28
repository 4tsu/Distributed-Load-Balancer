#include <iostream>
#include <fstream>
#include <random>
#include "systemparam.hpp"
#include "variables.hpp"

// =================================

void Variables::add_atoms(double x, double y, double z) {
    Atom a;
    a.x = x;
    a.y = y;
    a.z = z;
    a.vx = 0.0;
    a.vy = 0.0;
    a.vz = 0.0;
    atoms.push_back(a);
}

void Variables::export_cdview(void) {
    static int count = 0;
    char filename[256];
    sprintf(filename, "conf%03d.cdv", count);
    ++count;
    std::ofstream ofs(filename);
    int i = 0;
    for (auto &a : atoms) {
        ofs << i << " ";

        ofs << std::endl;
        ++i;
    }
}
