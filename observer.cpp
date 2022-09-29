#include <cstdio>
#include <fstream>
#include <iostream>
#include "systemparam.hpp"
#include "observer.hpp"
#include "variables.hpp"

// ======================================================

void Observer::export_cdview(std::vector<Atom> atoms) {
    static int count = 0;
    char filename[256];
    sprintf(filename, "conf%03d.cdv", count);
    ++count;
    std::ofstream ofs(filename);
    int i = 0;
    for (auto &a : atoms) {
        ofs << i   << " ";
        ofs << "0" << " ";
        ofs << a.x << " ";
        ofs << a.y << " ";
        ofs << std::endl;
        ++i;
    }
}

// ======================================================