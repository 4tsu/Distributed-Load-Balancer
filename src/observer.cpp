#include "observer.hpp"

// ======================================================

void Observer::export_cdview(Variables* vars, Systemparam* sysp, MPIinfo mi, int count_begin) {
    this->export_cdview_independent(vars, sysp, mi);
    MPI_Barrier(MPI_COMM_WORLD);
    this->concatenate_cdview(mi, count_begin);
}



void Observer::export_cdview_independent(Variables* vars, Systemparam* sysp, MPIinfo &mi) {
    static int count = 0;
    char filename[256];
#ifdef FS
    std::filesystem::create_directory("./cdv");
    sprintf(filename, "cdv/tconf%03d_%d.temp", count, mi.rank);
#else
    sprintf(filename, "tconf%03d_%d.temp", count, mi.rank);
#endif
    ++count;
    std::ofstream ofs(filename, std::ios::out);
    if (mi.rank==0) {
        ofs << "#box_sx=" << sysp->x_min << std::endl;
        ofs << "#box_sy=" << sysp->y_min << std::endl;
        ofs << "#box_ex=" << sysp->x_max << std::endl;
        ofs << "#box_ey=" << sysp->y_max << std::endl;
        ofs << "#box_sz=0" << std::endl;
        ofs << "#box_ez=0" << std::endl;
    }
    for (auto &a : vars->atoms) {
        ofs << a.id       << " ";
        ofs << mi.rank%9  << " ";   // cdviewの描画色が9色なので
        ofs << a.x        << " ";
        ofs << a.y        << " ";
        ofs << "0"        << " ";
        ofs << std::endl;
    }
}



void Observer::concatenate_cdview(MPIinfo &mi, int count_begin) {
    if(mi.rank==0) {
		static int cdv_count = count_begin;
		static int count = 0;
		char output[256];
#ifdef FS
		sprintf(output, "cdv/conf%03d.cdv", cdv_count);
#else
		sprintf(output, "conf%03d.cdv", cdv_count);
#endif
        std::ofstream ofs(output, std::ios::out);
		for (int i=0; i<mi.procs; i++) {
            char filename[256];
#ifdef FS
            sprintf(filename, "cdv/tconf%03d_%d.temp", count, i);
#else
            sprintf(filename, "tconf%03d_%d.temp", count, i);
#endif
            std::ifstream reading_file;
            reading_file.open(filename, std::ios::in);
            std::string line;
            while(std::getline(reading_file, line)) {
                ofs << line << std::endl;
            }
		}
#ifdef FS
        for (const auto & file : std::filesystem::directory_iterator("./cdv/")) {
            std::string path = file.path();
            size_t word_pos = path.find(".temp");
            if (word_pos != std::string::npos) {
                std::filesystem::remove(path);
                continue;
            }
        }
#endif
        count++;
        cdv_count++;
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



void Observer::export_checkpoint(std::string filename, const int step, Variables* vars, Systemparam* sysp, MPIinfo & mi) {
    this->checkpoint_independent(step, vars, sysp, mi);
    MPI_Barrier(MPI_COMM_WORLD);
    this->concatenate_checkpoint(filename, mi);
}



void Observer::checkpoint_independent(const int step, Variables* vars, Systemparam* sysp, MPIinfo &mi) {
    char filename[256];
#ifdef FS
    std::filesystem::create_directory("./ckpt");
    sprintf(filename, "ckpt/ckpt%d.temp", mi.rank);
#else
    sprintf(filename, "ckpt%d.temp", mi.rank);
#endif
    std::ofstream ofs(filename, std::ios::out);
    if (mi.rank==0) {
        ofs << "ITEM: TIMESTEP" << std::endl;
        ofs << step             << std::endl;
        ofs << "ITEM: NUMBER OF ATOMS" << std::endl;
        ofs << sysp->N                 << std::endl;
        ofs << "ITEM: BOX BOUNDS pp pp pp" << std::endl;
        ofs << sysp->x_min << " " << sysp->x_max << std::endl;
        ofs << sysp->y_min << " " << sysp->y_max << std::endl;
        ofs << "0 0"                             << std::endl;
        ofs << "ITEM: ATOMS x y z vx vy vz" << std::endl;
    }
    for (auto &a : vars->atoms) {
        ofs << a.x        << " ";
        ofs << a.y        << " ";
        ofs << "0"        << " ";
        ofs << a.vx        << " ";
        ofs << a.vy        << " ";
        ofs << "0"        << " ";
        ofs << std::endl;
    }
}

void Observer::concatenate_checkpoint(std::string filename, MPIinfo & mi) {
    if(mi.rank==0) {
#ifdef FS
        filename = "ckpt/" + filename;
#endif
        std::ofstream ofs(filename, std::ios::out);
		for (int i=0; i<mi.procs; i++) {
            char filename[256];
#ifdef FS
            sprintf(filename, "ckpt/ckpt%d.temp", i);
#else
            sprintf(filename, "ckpt%d.temp", i);
#endif
            std::ifstream reading_file;
            reading_file.open(filename, std::ios::in);
            std::string line;
            while(std::getline(reading_file, line)) {
                ofs << line << std::endl;
            }
		}
#ifdef FS
        for (const auto & file : std::filesystem::directory_iterator("./ckpt/")) {
            std::string path = file.path();
            size_t word_pos = path.find(".temp");
            if (word_pos != std::string::npos) {
                std::filesystem::remove(path);
                continue;
            }
        }
#endif
    }
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