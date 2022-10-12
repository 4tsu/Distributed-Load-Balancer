#include "systemparam.hpp"
#include "variables.hpp"

// =================================

/* やってみた残骸
MPI_Datatype MPI_Atom;
const int blocklength[] = {1,3,3};
const int displacements[] = {8,48,0};
MPI_Datatype types[] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
MPI_Type_create_struct(3, blocklength, displacements, types, &MPI_Atom);
MPI_Type_commit(&MPI_Atom);
*/

// ---------------------------------

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



double Variables::max_velocity(void) {
    double max_v = 0;
    for (auto &a : atoms) {
        double v = sqrt(a.vx*a.vx + a.vy*a.vy);
        if (max_v < v)
            max_v = v;
    }
    return max_v;
}



void Variables::pack_send_atoms(void) {
    this->send_atoms.clear();
    for (auto& one_send_list : this->send_list){
        std::vector<Atom*> one_send_atom;
        int atom_index = 0;
        for (int i=0; i<this->atoms.size(); i++) {
            if (atoms.at(i).id == one_send_list.at(atom_index)) {
                one_send_atom.push_back(&atoms.at(i));
                atom_index++;
            }
        }
        this->send_atoms.push_back(one_send_atom);
    }
}

// =================================