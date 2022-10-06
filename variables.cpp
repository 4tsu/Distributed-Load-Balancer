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



void Variables::set_initial_velocity(const double V0, MPIinfo mi) {
    std::mt19937 mt(mi.rank);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    double local_avx = 0.0;
    double local_avy = 0.0;
    for (auto &a : atoms) {
        double phi = 2.0 * ud(mt) * M_PI;
        double vx = V0 * cos(phi);
        double vy = V0 * sin(phi);
        a.vx = vx;
        a.vy = vy;
        local_avx += vx;
        local_avy += vy;
    }
    const int pn = atoms.size();
    local_avx /= static_cast<double>(pn);
    local_avy /= static_cast<double>(pn);
    
    // 全粒子平均速度を得るための通信
    double avx_sum, avy_sum;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&local_avx, &avx_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_avy, &avy_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double avx = avx_sum / static_cast<double>(mi.procs);
    double avy = avy_sum / static_cast<double>(mi.procs);

    // 全粒子平均速度を引いて、平均速度をゼロにする
    for (auto &a: atoms) {
        a.vx -= avx;
        a.vy -= avy;
    }
}

// =================================