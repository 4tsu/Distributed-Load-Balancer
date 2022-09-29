#include <cstdio>
#include <fstream>
#include <iostream>
#include "observer.hpp"

// ======================================================

void Observer::export_cdview(std::vector<Atom> atoms, Systemparam sysp) {
    static int count = 0;
    char filename[256];
    sprintf(filename, "conf%03d.cdv", count);
    ++count;
    std::ofstream ofs(filename);
    ofs << "#box_sx=" << sysp.x_min << std::endl;
    ofs << "#box_sy=" << sysp.y_min << std::endl;
    ofs << "#box_ex=" << sysp.x_max << std::endl;
    ofs << "#box_ey=" << sysp.y_max << std::endl;
    ofs << "#box_sz=0" << std::endl;
    ofs << "#box_ez=0" << std::endl;
    int i = 0;
    for (auto &a : atoms) {
        ofs << i   << " ";
        ofs << "0" << " ";
        ofs << a.x << " ";
        ofs << a.y << " ";
        ofs << "0" << " ";
        ofs << std::endl;
        ++i;
    }
}

// ======================================================