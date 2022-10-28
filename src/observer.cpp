#include "observer.hpp"

// ======================================================

void Observer::export_cdview(std::vector<Atom> atoms, Systemparam sysp, MPIinfo mi) {
    static int count = 0;
    char filename[256];
#ifdef FS
    std::filesystem::create_directory("./cdv");
    sprintf(filename, "cdv/conf%03d.cdv", count);
#else
    sprintf(filename, "conf%03d.cdv", count);
#endif
    ++count;
    std::ofstream ofs(filename, std::ios::app);
    if (mi.rank==0) {
        ofs << "#box_sx=" << sysp.x_min << std::endl;
        ofs << "#box_sy=" << sysp.y_min << std::endl;
        ofs << "#box_ex=" << sysp.x_max << std::endl;
        ofs << "#box_ey=" << sysp.y_max << std::endl;
        ofs << "#box_sz=0" << std::endl;
        ofs << "#box_ez=0" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (auto &a : atoms) {
        ofs << a.id       << " ";
        ofs << mi.rank%9  << " ";   // cdviewの描画色が9色なので
        ofs << a.x        << " ";
        ofs << a.y        << " ";
        ofs << "0"        << " ";
        ofs << std::endl;
    }
}



double Observer::kinetic_energy(Variables *vars, Systemparam *sysp) {
    double k = 0;
    for (auto& a : vars->atoms) {
        k += a.vx * a.vx;
        k += a.vy * a.vy;
    }

    double k_global;
    MPI_Allreduce(&k, &k_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    k_global /= static_cast<double>(sysp->N)*2;
    return k_global;
}



// ポテンシャルエネルギーの算出
/// ペアリストを使用するので、事前に構築しておくこと
double Observer::potential_energy(Variables *vars, PairList *pl, Systemparam *sysp) {
    double v = 0;

    for (auto& l : pl->list) {
        Atom ia = vars->atoms.at(l.i);
        Atom ja = vars->atoms.at(l.j);
        assert(ia.id == l.idi);
        assert(ja.id == l.idj);
        double dx = ja.x - ia.x;
        double dy = ja.y - ia.y;
        periodic_distance(dx, dy, sysp);
        double r = sqrt(dx*dx + dy*dy);
        if (r > sysp->cutoff)
            continue;
        double r6 = pow(r, 6);
        double r12 = r6*r6;
        v += 4.0 * (1.0/r12 - 1.0/r6) + sysp->C0;
    }
    for (std::size_t i=0; i<pl->other_list.size(); i++) {
        for (auto& l : pl->other_list.at(i)) {
            Atom ia = vars->atoms.at(l.i);
            Atom ja = vars->other_atoms.at(i).at(l.j);
            assert(ia.id == l.idi);
            assert(ja.id == l.idj);
            double dx = ja.x - ia.x;
            double dy = ja.y - ia.y;
            periodic_distance(dx, dy, sysp);
            double r = sqrt(dx*dx + dy*dy);
            if (r > sysp->cutoff)
                continue;
            double r6 = pow(r, 6);
            double r12 = r6*r6;
            v += 4.0 * (1.0/r12 - 1.0/r6) + sysp->C0;
        }
    }
    double v_global;
    MPI_Allreduce(&v, &v_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    v_global /= static_cast<double>(sysp->N);
    return v_global;
}

// ======================================================

void export_three(const std::string filename, const int s, const double a, const double b, const double c) {
    std::ofstream ofs(filename, std::ios::app);
    ofs << s << " ";
    ofs << a << " ";
    ofs << b << " ";
    ofs << c << " ";
    ofs << "\n";
}

// =================================================