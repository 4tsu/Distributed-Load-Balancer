#include "sdd.hpp"



// ============================================

Sdd::Sdd(const int sdd_type) {
    this->sdd_type = sdd_type;
}

Sdd::~Sdd(void) {
}

// -----------------------------------------------------

void Sdd::init(Variables* vars, Systemparam* sysp
               , const MPIinfo &mi, SubRegion* sr) {
    if        (sdd_type==0) {
        calc_bounds(sysp, mi);
        return;

    } else if (sdd_type==1) {
        global_sort(vars, sysp, mi, sr);

    } else if (sdd_type==2) {
        calc_bounds(sysp, mi);
        simple(vars, sysp, mi);
        voronoi_init(vars, sysp, mi, sr);
        return;
    }
}



void Sdd::run(Variables* vars, Systemparam* sysp, const MPIinfo &mi, SubRegion* sr) {

    if        (sdd_type==0) {
        simple(vars, sysp, mi);

    } else if (sdd_type==1) {
        global_sort(vars, sysp, mi, sr);

    } else if (sdd_type==2) {
        voronoi(vars, sysp, mi, sr);
    }
}



unsigned long Sdd::ideal(Systemparam* sysp, const MPIinfo &mi) {
    return static_cast<unsigned long>(sysp->N/mi.procs);
}



void Sdd::migrate_atoms(std::vector<std::vector<Atom>> migration_atoms, Variables* vars, const MPIinfo &mi) {
    // 移動する原子の通信
    // 最初にAllgatherで移動粒子数のテーブルを共有（1次元で送って戻す）
    // 次に必要なペアでだけ通信をする

    std::vector<unsigned long> num_migration(mi.procs);
    for (int i=0; i<mi.procs; i++) {
        num_migration.at(i) = migration_atoms.at(i).size();
    }
    // 移動する原子の通信
    std::vector<std::vector<unsigned long>> migration_table(mi.procs, std::vector<unsigned long>(mi.procs));
    std::vector<unsigned long> recvbuf(mi.procs*mi.procs);
    MPI_Allgather(num_migration.data(), mi.procs, MPI_UNSIGNED_LONG, recvbuf.data(), mi.procs, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    int bias = 0;
    unsigned long sum_recv = 0;
    for (auto& one_mt : migration_table) {
        std::copy(recvbuf.begin()+bias, recvbuf.begin()+bias+mi.procs, one_mt.begin());
        bias += mi.procs;
        sum_recv += one_mt.at(mi.rank);
    }

    MPI_Request ireq;
    std::vector<MPI_Request> mpi_send_requests;
    std::vector<MPI_Request> mpi_recv_requests;
    unsigned long nm_size;
    std::vector<Atom> recvatoms(sum_recv);
    unsigned long recv_head = 0;
    for (int p=0; p<mi.procs; p++) {
        if (num_migration.at(p)!=0) {
            nm_size = num_migration.at(p)*sizeof(Atom);
            MPI_Isend(migration_atoms.at(p).data(), nm_size, MPI_CHAR, p, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);
        }
    }

    for (int p=0; p<mi.procs; p++) {
        unsigned long m = migration_table.at(p).at(mi.rank);
        if (m!=0) {
            nm_size = m*sizeof(Atom);
            MPI_Irecv(&recvatoms.at(recv_head), nm_size, MPI_CHAR, p, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
            recv_head += m;
        }
    }
    MPI_Status st;
    for (auto& ireq : mpi_recv_requests) {
        MPI_Wait(&ireq, &st);
    }
    for (auto& ireq : mpi_send_requests) {
        MPI_Wait(&ireq, &st);
    }

    for (Atom atom : recvatoms) {
        vars->atoms.push_back(atom);
    }
}



void Sdd::calc_bounds(Systemparam* sysp, const MPIinfo &mi) {
    int ix = mi.rank%mi.npx;
    int iy = mi.rank/mi.npx;
    const double lpx = sysp->xl/static_cast<double>(mi.npx);
    const double lpy = sysp->yl/static_cast<double>(mi.npy);
    this->bottom = lpy*static_cast<double>(iy) + sysp->y_min;
    this->top    = lpy*static_cast<double>(iy+1) + sysp->y_min;
    this->left   = lpx*static_cast<double>(ix) + sysp->x_min;
    this->right  = lpx*static_cast<double>(ix+1) + sysp->x_min;
}



void Sdd::simple(Variables* vars, Systemparam* sysp, const MPIinfo &mi) {
    std::vector<std::vector<Atom>> migration_atoms(mi.procs);
    std::vector<Atom> new_atoms;
    const double lpx = sysp->xl/static_cast<double>(mi.npx);
    const double lpy = sysp->yl/static_cast<double>(mi.npy);
    for (const Atom atom : vars->atoms) {
        if (atom.x>right || atom.x<left || atom.y>top || atom.y<bottom) {
            int new_rank = static_cast<int>(floor((atom.y-sysp->y_min)/lpy)*mi.npx + floor((atom.x-sysp->x_min)/lpx));
            migration_atoms.at(new_rank).push_back(atom);
        } else {
            new_atoms.push_back(atom);
        }
    }
    vars->atoms = new_atoms;
    migrate_atoms(migration_atoms, vars, mi);
}
        
        

// 未実装
#pragma GCC diagnostic ignored "-Wunused-parameter"
void Sdd::global_sort(Variables* vars, Systemparam* sysp, const MPIinfo &mi, SubRegion* sr) {

}



void Sdd::voronoi_init(Variables* vars, Systemparam* sysp, const MPIinfo &mi, SubRegion* sr) {

}



void Sdd::voronoi(Variables* vars, Systemparam* sysp, const MPIinfo &mi, SubRegion* sr) {

}
#pragma GCC diagnostic warning "-Wunused-parameter"

// ============================================