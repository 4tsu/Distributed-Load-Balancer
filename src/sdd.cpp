#include "sdd.hpp"

namespace sysp = systemparam;

// ============================================

Sdd::Sdd(const int sdd_type) {
    this->sdd_type = sdd_type;
}

Sdd::~Sdd(void) {
}

// -----------------------------------------------------

void Sdd::init(Variables* vars
               , const MPIinfo &mi, SubRegion* sr) {
    if        (sdd_type==0) {
        calc_bounds(mi);
        return;

    } else if (sdd_type==1) {
        global_sort(vars, mi);

    } else if (sdd_type==2) {
        std::vector<double> v(6);
        v = calc_limit(vars);
        this->left   = v.at(0);
        this->right  = v.at(1);
        this->front  = v.at(2);
        this->back   = v.at(3);
        this->bottom = v.at(4);
        this->top    = v.at(5);
        voronoi_init(vars, mi, sr);
        return;

    } else if (sdd_type==3) {
        return;

    } else if (sdd_type==4) {
        odp_init(vars, mi);
        return;

    } else if (sdd_type==5) {
        sb_init(vars, mi);
        return;
    }
}



void Sdd::run(Variables* vars, const MPIinfo &mi, SubRegion* sr) {

    if        (sdd_type==0) {
        simple(vars, mi);

    } else if (sdd_type==1) {
        global_sort(vars, mi);

    } else if (sdd_type==2) {
        voronoi(vars, mi, sr, 300, 0.050, 0.02);
    
    } else if (sdd_type==3) {
        rcb(vars, mi);

    } else if (sdd_type==4) {
        one_d_parallel(vars, mi, 300, 0.0001, 0.02);

    } else if (sdd_type==5) {
        skew_boundary(vars, mi, 300, 0.0001, 0.02);
    }
}



unsigned long Sdd::ideal(const MPIinfo &mi) {
    return static_cast<unsigned long>(sysp::N/mi.procs);
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
    int head = 0;
    unsigned long sum_recv = 0;
    for (auto& one_mt : migration_table) {
        std::copy(recvbuf.begin()+head, recvbuf.begin()+head+mi.procs, one_mt.begin());
        head += mi.procs;
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



void Sdd::calc_bounds(const MPIinfo &mi) {
    int ix = mi.rank%mi.npx;
    int iy = (mi.rank/mi.npx)%mi.npy;
    int iz = mi.rank/(mi.npx*mi.npy);
    const double lpx = sysp::xl/static_cast<double>(mi.npx);
    const double lpy = sysp::yl/static_cast<double>(mi.npy);
    const double lpz = sysp::zl/static_cast<double>(mi.npz);
    this->left   = lpx*static_cast<double>(ix) + sysp::x_min;
    this->right  = lpx*static_cast<double>(ix+1) + sysp::x_min;
    this->back   = lpy*static_cast<double>(iy+1) + sysp::y_min;
    this->front  = lpy*static_cast<double>(iy) + sysp::y_min;
    this->bottom = lpz*static_cast<double>(iz) + sysp::z_min;
    this->top    = lpz*static_cast<double>(iz+1) + sysp::z_min;
}



void Sdd::simple(Variables* vars, const MPIinfo &mi) {
    std::vector<std::vector<Atom>> migration_atoms(mi.procs);
    std::vector<Atom> new_atoms;
    const double lpx = sysp::xl/static_cast<double>(mi.npx);
    const double lpy = sysp::yl/static_cast<double>(mi.npy);
    const double lpz = sysp::zl/static_cast<double>(mi.npz);
    for (const Atom atom : vars->atoms) {
        if (atom.x>right || atom.x<left || atom.y>back || atom.y<front || atom.z>top || atom.z<bottom) {
            int new_rank = static_cast<int>(floor((atom.z-sysp::z_min)/lpz)*mi.npx*mi.npy + floor((atom.y-sysp::y_min)/lpy)*mi.npx + floor((atom.x-sysp::x_min)/lpx));
            migration_atoms.at(new_rank).push_back(atom);
        } else {
            new_atoms.push_back(atom);
        }
    }
    vars->atoms = new_atoms;
    migrate_atoms(migration_atoms, vars, mi);
}
        
        

// X,Y軸方向のそれぞれについて、各プロセスの粒子数が均等になるように壁を平行に動かす
void Sdd::global_sort(Variables* vars, const MPIinfo &mi) {
    std::vector<std::vector<Atom>> migration_atoms(mi.procs);
    unsigned long l = vars->atoms.size();
    unsigned long max_l;
    MPI_Allreduce(&l, &max_l, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    vars->atoms.resize(max_l);
    unsigned long diff_n = max_l - l;
    std::vector<unsigned long> diffs(mi.procs);
    MPI_Gather(&diff_n, 1, MPI_UNSIGNED_LONG, diffs.data(), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    std::vector<Atom> all_atoms;
    if (mi.rank==0) {
        all_atoms.resize(max_l*mi.procs);
    }
    unsigned long s = max_l*sizeof(Atom);
    MPI_Gather(vars->atoms.data(), s, MPI_CHAR, all_atoms.data(), s, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (mi.rank==0) {
        auto compare_x  = [](const Atom & a1, const Atom & a2) {return a1.x < a2.x;};
        auto compare_y  = [](const Atom & a1, const Atom & a2) {return a1.y < a2.y;};
        auto compare_z  = [](const Atom & a1, const Atom & a2) {return a1.z < a2.z;};
        unsigned long sum = 0;
        for (int i=0; i<mi.procs; i++) {
            sum += max_l;
            all_atoms.erase(all_atoms.begin()+sum-diffs.at(i), all_atoms.begin()+sum);
            sum -= diffs.at(i);
        }
        // まずはZ軸についてソート
        std::sort(all_atoms.begin(), all_atoms.end(), compare_z);
        unsigned long naz = std::floor(sysp::N/mi.npz) + 1;
        unsigned long head = 0;
        unsigned long z_count = naz;
        for (int i=0; i<mi.npz; i++) {
            if (head+z_count>sysp::N)
                z_count = sysp::N - head;

            // Y軸についてソート
            std::sort(all_atoms.begin()+head, all_atoms.begin()+head+z_count, compare_y);
            unsigned long nay = std::floor(z_count/mi.npy) + 1;
            unsigned long y_count = nay;
            unsigned long y_sum = 0;
            for (int j=0; j<mi.npy; j++) {
                if (y_sum+y_count>z_count)
                    y_count = z_count-y_sum;

                // Z軸についてソート
                std::sort(all_atoms.begin()+head, all_atoms.begin()+head+y_count, compare_x);
                unsigned long nax = std::floor(y_count/mi.npx) + 1;
                unsigned long x_count = nax;
                unsigned long x_sum = 0;
                for (int k=0; k<mi.npx; k++) {
                    int p = k + j*mi.npx + i*(mi.npx*mi.npy);
                    if (x_sum+x_count>y_count)
                        x_count = y_count-x_sum;
                    migration_atoms.at(p).resize(x_count);
                    std::copy(all_atoms.begin()+head+y_sum+x_sum, all_atoms.begin()+head+y_sum+x_sum+x_count, migration_atoms.at(p).data());
                    x_sum += x_count;
                }
                y_sum += y_count;
            }
            head += z_count;
        }
    vars->atoms.clear();
    vars->atoms.resize(migration_atoms.at(0).size());
    std::copy(migration_atoms.at(0).begin(), migration_atoms.at(0).end(), vars->atoms.data());
    migration_atoms.at(0).clear();
    } else {
        for (auto & oma : migration_atoms)
                oma.clear();
        vars->atoms.clear();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    migrate_atoms(migration_atoms, vars, mi);
}



// 領域が空だとvoronoiに参加できないので、最も重い領域たちから粒子を分けてもらう
void Sdd::voronoi_init(Variables* vars, const MPIinfo &mi, SubRegion* sr) {
    unsigned long na = vars->number_of_atoms();
    Workload wl;
    wl.rank = mi.rank;
    wl.counts = na;
    std::vector<Workload> loads(mi.procs);
    MPI_Allgather(&wl, sizeof(wl), MPI_CHAR, loads.data(), sizeof(wl), MPI_CHAR, MPI_COMM_WORLD);
    
    // 昇順比較
    auto compare_wl = [](const Workload & wl1, const Workload & wl2) {return wl1.counts > wl2.counts;};
    std::sort(loads.begin(), loads.end(), compare_wl);

    while (loads.at(mi.procs-1).counts==0) {
        int zero_regions = 0;
        int non_zero_regions = 0;
        int my_position;
        for (int i=0; i<mi.procs; i++) {
            Workload wlt = loads.at(i);
            if (wlt.counts==0) {
                zero_regions++;
            } else {
                non_zero_regions++;
            }
            if(wlt.rank==mi.rank)
                my_position = i;
        }
        int zero_head = mi.procs - zero_regions;
        if (non_zero_regions < zero_regions) {
            zero_regions = non_zero_regions;
        }
        if (zero_head <= my_position && my_position < zero_head+zero_regions) {
            // 相手から半分もらう作業
            int source_proc = loads.at(my_position - zero_head).rank;
            unsigned long num_recv;
            MPI_Status st;
            MPI_Recv(&num_recv, 1, MPI_UNSIGNED_LONG, source_proc, 0, MPI_COMM_WORLD, &st);
            vars->atoms.resize(num_recv);
            MPI_Recv(vars->atoms.data(), num_recv*sizeof(Atom), MPI_CHAR, source_proc, 0, MPI_COMM_WORLD, &st);
        } else if (my_position < zero_regions) {
            // 相手に半分渡す作業
            int target_proc = loads.at(zero_head + my_position).rank;
            std::vector<double> ls = calc_limit(vars);
            double center_line = (ls.at(0)+ls.at(1))/2;
            std::vector<Atom> send_atoms;
            std::vector<Atom> new_atoms;
            for (const auto &a : vars->atoms) {
                if (a.x<center_line) {
                    send_atoms.push_back(a);
                } else {
                    new_atoms.push_back(a);
                }
            }
            unsigned long num_send = send_atoms.size();
            MPI_Send(&num_send, 1, MPI_UNSIGNED_LONG, target_proc, 0, MPI_COMM_WORLD);
            MPI_Send(send_atoms.data(), num_send*sizeof(Atom), MPI_CHAR, target_proc, 0, MPI_COMM_WORLD);
            vars->atoms.clear();
            vars->atoms.resize(new_atoms.size());
            vars->atoms = new_atoms;
        }
        wl.counts = vars->number_of_atoms();
        loads.clear(); 
        loads.resize(mi.procs);
        MPI_Allgather(&wl, sizeof(wl), MPI_CHAR, loads.data(), sizeof(wl), MPI_CHAR, MPI_COMM_WORLD);
        std::sort(loads.begin(), loads.end(), compare_wl);
    } 
    sr->bias = 0;
    this->ideal_count = ideal(mi);
}



void Sdd::voronoi(Variables* vars, const MPIinfo &mi, SubRegion* sr,
                  int iteration, double alpha, double early_stop_range) {
    std::vector<unsigned long> counts(mi.procs);
    unsigned long ideal_count_max = std::ceil(ideal_count*(1+early_stop_range));

    all_biases.clear();
    all_biases.resize(mi.procs);
    MPI_Allgather(&sr->bias, 1, MPI_DOUBLE, all_biases.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    sr->calc_center(vars);
    sr->calc_radius(vars);
    Sdd::voronoi_allocate(vars, mi, sr);
    sr->calc_center(vars);
    sr->calc_radius(vars);

    // ボロノイの分割最適化の様子出力
    // voronoi_figure(vars, mi);

    for (int s=1; s<=iteration; s++) {
        unsigned long pn = vars->number_of_atoms();
        MPI_Allgather(&pn, 1, MPI_UNSIGNED_LONG, counts.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        unsigned long max_count = *std::max_element(counts.begin(), counts.end());

        /*
        // ボロノイ最適化中のロードバランスの変化を出力
        if (mi.rank==0) {
            std::printf("%d ", s);
            for (auto c : counts) {
                std::printf("%ld ", c);
            }
            std::printf("\n");
        }
        */
        
        // early stop
       if (max_count <= ideal_count_max || max_count-ideal_count < 10) {
            // std::fprintf(stderr, "***early stop (iter #%d)***\n", s);
            break;
        }

        sr->make_dplist(mi, vars);
        std::vector<DomainPair> dp_all;
        for (auto dp:sr->dplist)
            dp_all.push_back(dp);
        for (auto dp:sr->dplist_reverse)
            dp_all.push_back(dp);
        
        for (auto dp : dp_all) {
            int i = dp.i;
            int j = dp.j;
            assert(i==mi.rank);
            double c_i = static_cast<double>(counts.at(i));
            double c_j = static_cast<double>(counts.at(j));
            double db = pow((alpha*(c_i-c_j)), 3);
            sr->bias -= db;
        }
        std::vector<double> limits(6);
        limits = calc_limit(vars);
        double min_xyzl = std::min({limits.at(1)-limits.at(0), limits.at(3)-limits.at(2), limits.at(5)-limits.at(4)});
        if (sr->bias < -sr->radius) {
            sr->bias = -sr->radius;
        } else if (sr->bias > min_xyzl) {
            sr->bias = min_xyzl;
        }
        all_biases.clear();
        all_biases.resize(mi.procs);
        MPI_Allgather(&sr->bias, 1, MPI_DOUBLE, all_biases.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

        /*
        // ボロノイ最適化中のバイアスの変化を出力
        if (mi.rank==0) {
            std::printf("%d ", s);
            for (auto b : all_biases) {
                std::printf("%lf ", b);
            }
            std::printf("\n");
        }
        */
        
        voronoi_allocate(vars, mi, sr);
        sr->calc_center(vars);
        sr->calc_radius(vars);

        // voronoi_figure(vars, mi);

    }
}



void Sdd::voronoi_allocate(Variables* vars, const MPIinfo &mi, SubRegion* sr) {
    std::vector<std::vector<Atom>> migration_atoms(mi.procs);
    std::vector<Atom> new_atoms;
    sr->communicate_centradi(mi);

    std::vector<DomainPair> dp_all;
    for (auto dp:sr->dplist) {
        assert(dp.i == mi.rank);
        dp_all.push_back(dp);
    }
    for (auto dp:sr->dplist_reverse) {
        assert(dp.i == mi.rank);
        dp_all.push_back(dp);
    }
    int closest_proc;
    double min_distance;
    for (const Atom atom : vars->atoms) {
        min_distance = sysp::xl*sysp::xl+sysp::yl*sysp::yl+sysp::zl*sysp::zl;
        for (auto dp : dp_all) {
            center_atom_distance(dp.j, min_distance, closest_proc, atom, sr);
        }
        center_atom_distance(mi.rank, min_distance, closest_proc, atom, sr);
        migration_atoms.at(closest_proc).push_back(atom);
    }
    vars->atoms.clear();
    vars->atoms.resize(migration_atoms.at(mi.rank).size());
    std::copy(migration_atoms.at(mi.rank).begin(), migration_atoms.at(mi.rank).end(), vars->atoms.data());
    migration_atoms.at(mi.rank).clear();
    MPI_Barrier(MPI_COMM_WORLD);
    migrate_atoms(migration_atoms, vars, mi);
}



void Sdd::center_atom_distance(int rank, double & min_distance, int & closest_proc,
                               const Atom atom, SubRegion* sr) {
    double dx = atom.x - sr->centers.at(rank).at(0);
    double dy = atom.y - sr->centers.at(rank).at(1);
    double dz = atom.z - sr->centers.at(rank).at(2);
    periodic_distance(dx, dy, dz);
    double r2 = dx*dx + dy*dy + dz*dz - all_biases.at(rank);
    if (r2<min_distance) {
        min_distance = r2;
        closest_proc = rank;
    }
}



void Sdd::voronoi_figure(Variables* vars, const MPIinfo &mi) {
    static int s = 0;
    char filename[256];
    std::sprintf(filename, "voronoi_%03d.cdv", s);
    std::ofstream ofs(filename, std::ios::app);
    if (mi.rank==0) {
        ofs << "#box_sx=" << sysp::x_min << std::endl;
        ofs << "#box_sy=" << sysp::y_min << std::endl;
        ofs << "#box_sz=" << sysp::z_min << std::endl;
        ofs << "#box_ex=" << sysp::x_max << std::endl;
        ofs << "#box_ey=" << sysp::y_max << std::endl;
        ofs << "#box_ez=" << sysp::z_max << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (auto &a : vars->atoms) {
        ofs << a.id       << " ";
        ofs << mi.rank%9  << " ";   // cdviewの描画色が9色なので
        ofs << a.x        << " ";
        ofs << a.y        << " ";
        ofs << a.z        << " ";
        ofs << std::endl;
    }
    s++;
}



void Sdd::rcb(Variables* vars, const MPIinfo &mi) {
    // 初期分割
    std::vector<std::vector<Atom>> migration_atoms(mi.procs);
    if (mi.rank!=0) {
        migration_atoms.at(0).resize(vars->number_of_atoms());
        std::copy(vars->atoms.begin(), vars->atoms.end(), migration_atoms.at(0).begin());
        vars->atoms.clear();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    migrate_atoms(migration_atoms, vars, mi);

    std::vector<int> directions(mi.procs);
    std::fill(directions.begin(), directions.end(), 0);
            
    // 準備
    auto compare_x  = [](const Atom & a1, const Atom & a2) {return a1.x < a2.x;};
    auto compare_y  = [](const Atom & a1, const Atom & a2) {return a1.y < a2.y;};
    auto compare_z  = [](const Atom & a1, const Atom & a2) {return a1.z < a2.z;};

    for (int i=1; i<mi.procs; i++) {
        std::vector<unsigned long> counts(mi.procs);
        unsigned long pn = vars->number_of_atoms();
        MPI_Allgather(&pn, 1, MPI_UNSIGNED_LONG, counts.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        unsigned long s = std::accumulate(counts.begin(), counts.end(), 0);
        assert(s==sysp::N);
        std::vector<unsigned long>::iterator max_it = std::max_element(counts.begin(), counts.end());
        int target = static_cast<int>(std::distance(counts.begin(), max_it));

        migration_atoms.clear();
        migration_atoms.resize(mi.procs);
        
        if (target==mi.rank) {
            // x軸分割
            if (directions.at(target) == 0) {
                directions.at(target) = 1;
                directions.at(i)      = 1;
                std::sort(vars->atoms.begin(), vars->atoms.end(), compare_x);
            // y軸分割
            } else if (directions.at(target) == 1) {
                directions.at(target) = 2;
                directions.at(i)      = 2;
                std::sort(vars->atoms.begin(), vars->atoms.end(), compare_y);
            } else {
                directions.at(target) = 0;
                directions.at(i)      = 0;
                std::sort(vars->atoms.begin(), vars->atoms.end(), compare_z);
            }
            unsigned long hi = vars->number_of_atoms()/2;
            migration_atoms.at(i).resize(hi);
            std::copy(vars->atoms.begin(), vars->atoms.begin()+hi, migration_atoms.at(i).begin());
            vars->atoms.erase(vars->atoms.begin(), vars->atoms.begin()+hi);
            vars->atoms.shrink_to_fit();
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        migrate_atoms(migration_atoms, vars, mi);
    }
}



void Sdd::set_np(Neighbor_Process &proc, int rank, unsigned long counts, double left, double right) {
    proc.rank = rank;
    proc.counts = static_cast<double>(counts);
    proc.left  = left;
    proc.right = right;
}



// 初期分割
void Sdd::odp_init(Variables* vars, const MPIinfo &mi) {
    double lxp = sysp::xl/static_cast<double>(mi.procs);
    this->back   = sysp::y_max;
    this->front  = sysp::y_min;
    this->top    = sysp::z_max;
    this->bottom = sysp::z_min;
    double rank_double = static_cast<double>(mi.rank);
    this->right  = (rank_double+1)*lxp + sysp::x_min;
    this->left   = rank_double*lxp + sysp::x_min;

    std::vector<std::vector<Atom>> migration_atoms(mi.procs);
    for (auto& a : vars->atoms) {
        migration_atoms.at(static_cast<int>((a.x-sysp::x_min)/lxp)).push_back(a);
    }
    vars->atoms.clear();
    vars->atoms.resize(migration_atoms.at(mi.rank).size());
    std::copy(migration_atoms.at(mi.rank).begin(), migration_atoms.at(mi.rank).end(), vars->atoms.begin());
    migration_atoms.at(mi.rank).clear();
    MPI_Barrier(MPI_COMM_WORLD);
    migrate_atoms(migration_atoms, vars, mi);

    std::vector<unsigned long> counts(mi.procs);
    unsigned long pn = vars->number_of_atoms();
    MPI_Allgather(&pn, 1, MPI_UNSIGNED_LONG, counts.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    unsigned long s = std::accumulate(counts.begin(), counts.end(), 0);
    assert(s==sysp::N);
}



void Sdd::one_d_parallel(Variables* vars, const MPIinfo &mi, int iteration, double alpha, double early_stop_range) {
    // いったん領域境界をはみ出た粒子を整理する
    unsigned long ideal_count_max = std::ceil(static_cast<double>(sysp::N/mi.procs)*(1+early_stop_range));
    std::vector<std::vector<Atom>> migration_atoms(mi.procs);
    Neighbor_Process my_proc, next_proc, previous_proc;
    unsigned long pn = vars->number_of_atoms();
    set_np(my_proc, mi.rank, pn, this->left, this->right);
    std::vector<Neighbor_Process> all_procs(mi.procs);
    MPI_Allgather(&my_proc, sizeof(Neighbor_Process), MPI_CHAR, all_procs.data(), sizeof(Neighbor_Process), MPI_CHAR, MPI_COMM_WORLD);
    for (auto& a : vars->atoms) {
        if (my_proc.left <= a.x && a.x <= my_proc.right) {
            migration_atoms.at(mi.rank).push_back(a);
        } else {
            for (int i=1; i<=mi.procs/2; i++) {
                int next = ((mi.rank-i)%mi.procs+mi.procs)%mi.procs;
                int prev = (mi.rank+i)%mi.procs;
                if (all_procs.at(next).left <= a.x && a.x <= all_procs.at(next).right) {
                    migration_atoms.at(next).push_back(a);
                    break;
                } else if (all_procs.at(prev).left <= a.x && a.x <= all_procs.at(prev).right) {
                    migration_atoms.at(prev).push_back(a);
                    break;
                }
            }
        }
    }
    vars->atoms.clear();
    vars->atoms.resize(migration_atoms.at(mi.rank).size());
    std::copy(migration_atoms.at(mi.rank).begin(), migration_atoms.at(mi.rank).end(), vars->atoms.begin());
    migration_atoms.at(mi.rank).clear();
    MPI_Barrier(MPI_COMM_WORLD);
    migrate_atoms(migration_atoms, vars, mi);

    // iteration
    for (int s=0; s<iteration; s++) {
        pn = vars->number_of_atoms();
        set_np(my_proc, mi.rank, pn, this->left, this->right);
        std::vector<MPI_Request> mpi_sendlist, mpi_recvlist;
        MPI_Request ireq;
        MPI_Status st;

        if (0<mi.rank) {
            MPI_Isend(&my_proc, sizeof(Neighbor_Process), MPI_CHAR, mi.rank-1, 0, MPI_COMM_WORLD, &ireq);
            mpi_sendlist.push_back(ireq);
            MPI_Irecv(&previous_proc, sizeof(Neighbor_Process), MPI_CHAR, mi.rank-1, 0, MPI_COMM_WORLD, &ireq);
            mpi_recvlist.push_back(ireq);
        }
        if (mi.rank<mi.procs-1) {
            MPI_Isend(&my_proc, sizeof(Neighbor_Process), MPI_CHAR, mi.rank+1, 0, MPI_COMM_WORLD, &ireq);
            mpi_sendlist.push_back(ireq);
            MPI_Irecv(&next_proc, sizeof(Neighbor_Process), MPI_CHAR, mi.rank+1, 0, MPI_COMM_WORLD, &ireq);
            mpi_recvlist.push_back(ireq);
        }
        for (auto& ireq : mpi_sendlist) {
            MPI_Wait(&ireq, &st);
        }
        for (auto& ireq : mpi_recvlist) {
            MPI_Wait(&ireq, &st);
        }
        
        double dx;
        double left_limit;
        double right_limit;
        if (0<mi.rank) {
            dx = -1.0*std::pow(alpha*(previous_proc.counts - my_proc.counts), 1);
            left_limit = -1.0*(previous_proc.right - previous_proc.left)/2.0;
            right_limit = (my_proc.right - my_proc.left)/2.0;
            if (dx < left_limit) {
                dx = left_limit;
            } else if (dx > right_limit) {
                dx = right_limit;
            }
            this->left += dx;
        }
        if (mi.rank<mi.procs-1) {
            dx = -1.0*std::pow(alpha*(my_proc.counts - next_proc.counts), 1);
            left_limit = -1.0*(my_proc.right - my_proc.left)/2.0;
            right_limit = (next_proc.right - next_proc.left)/2.0;
            if (dx < left_limit) {
                dx = left_limit;
            } else if (dx > right_limit) {
                dx = right_limit;
            }
            this->right += dx;
        }

        migration_atoms.clear();
        migration_atoms.resize(mi.procs);
        for (auto& a : vars->atoms) {
            if (a.x < this->left) {
                migration_atoms.at(mi.rank-1).push_back(a);
            } else if(a.x > this->right) {
                migration_atoms.at(mi.rank+1).push_back(a);
            } else {
                migration_atoms.at(mi.rank).push_back(a);
            }
        }

        vars->atoms.clear();
        vars->atoms.resize(migration_atoms.at(mi.rank).size());
        std::copy(migration_atoms.at(mi.rank).begin(), migration_atoms.at(mi.rank).end(), vars->atoms.begin());
        migration_atoms.at(mi.rank).clear();
        MPI_Barrier(MPI_COMM_WORLD);
        migrate_atoms(migration_atoms, vars, mi);

        std::vector<unsigned long> counts(mi.procs);
        pn = vars->number_of_atoms();
        MPI_Allgather(&pn, 1, MPI_UNSIGNED_LONG, counts.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        unsigned long sum_counts = std::accumulate(counts.begin(), counts.end(), 0);
        assert(sum_counts==sysp::N);
        
        // early stop
        unsigned long max_count = *std::max_element(counts.begin(), counts.end());
        if (max_count <= ideal_count_max || max_count-ideal_count < 10) {
            // std::fprintf(stderr, "***early stop (iter #%d)***\n", s);
            break;
        }
    }
}



void Sdd::sb_init(Variables* vars, const MPIinfo& mi) {
    // 初期分割は等間隔
    calc_bounds(mi);
    simple(vars, mi);
    // ただし、境界の定義だけ変更する
    double py = static_cast<double>(mi.rank/mi.npx);
    this->right += py*sysp::xl - sysp::x_min;
    this->left  += py*sysp::xl - sysp::x_min;

    std::vector<unsigned long> counts(mi.procs);
    unsigned long pn = vars->number_of_atoms();
    MPI_Allgather(&pn, 1, MPI_UNSIGNED_LONG, counts.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    unsigned long s = std::accumulate(counts.begin(), counts.end(), 0);
    assert(s==sysp::N);
}



void Sdd::skew_boundary(Variables* vars, const MPIinfo& mi, int iteration, double alpha, double early_stop_range) {
    unsigned long ideal_count_max = std::ceil(static_cast<double>(sysp::N/mi.procs)*(1+early_stop_range));
    const double lpy = sysp::yl/static_cast<double>(mi.npy);
    const double lpz = sysp::zl/static_cast<double>(mi.npz);
    
    // いったん領域境界をはみ出た粒子を整理する
    std::vector<std::vector<Atom>> migration_atoms(mi.procs);
    Neighbor_Process my_proc, next_proc, previous_proc;
    unsigned long pn = vars->number_of_atoms();
    set_np(my_proc, mi.rank, pn, this->left, this->right);
    std::vector<Neighbor_Process> all_procs(mi.procs);
    MPI_Allgather(&my_proc, sizeof(Neighbor_Process), MPI_CHAR, all_procs.data(), sizeof(Neighbor_Process), MPI_CHAR, MPI_COMM_WORLD);
    for (auto& a : vars->atoms) {
        double sbx = std::floor((a.z-sysp::z_min)/lpz)*mi.npy*sysp::xl + std::floor((a.y-sysp::y_min)/lpy)*sysp::xl + a.x-sysp::x_min;
        if (my_proc.left <= sbx && sbx <= my_proc.right) {
            migration_atoms.at(mi.rank).push_back(a);
        } else {
            for (int i=1; i<=mi.procs/2; i++) {
                int next = ((mi.rank-i)%mi.procs+mi.procs)%mi.procs;
                int prev = (mi.rank+i)%mi.procs;
                if (all_procs.at(next).left <= sbx && sbx <= all_procs.at(next).right) {
                    migration_atoms.at(next).push_back(a);
                    break;
                } else if (all_procs.at(prev).left <= sbx && sbx <= all_procs.at(prev).right) {
                    migration_atoms.at(prev).push_back(a);
                    break;
                }
            }
        }
    }
    vars->atoms.clear();
    vars->atoms.resize(migration_atoms.at(mi.rank).size());
    std::copy(migration_atoms.at(mi.rank).begin(), migration_atoms.at(mi.rank).end(), vars->atoms.begin());
    migration_atoms.at(mi.rank).clear();
    MPI_Barrier(MPI_COMM_WORLD);
    migrate_atoms(migration_atoms, vars, mi);

    // iteration
    for (int s=0; s<iteration; s++) {
        pn = vars->number_of_atoms();
        set_np(my_proc, mi.rank, pn, this->left, this->right);
        std::vector<MPI_Request> mpi_sendlist, mpi_recvlist;
        MPI_Request ireq;
        MPI_Status st;

        if (0<mi.rank) {
            MPI_Isend(&my_proc, sizeof(Neighbor_Process), MPI_CHAR, mi.rank-1, 0, MPI_COMM_WORLD, &ireq);
            mpi_sendlist.push_back(ireq);
            MPI_Irecv(&previous_proc, sizeof(Neighbor_Process), MPI_CHAR, mi.rank-1, 0, MPI_COMM_WORLD, &ireq);
            mpi_recvlist.push_back(ireq);
        }
        if (mi.rank<mi.procs-1) {
            MPI_Isend(&my_proc, sizeof(Neighbor_Process), MPI_CHAR, mi.rank+1, 0, MPI_COMM_WORLD, &ireq);
            mpi_sendlist.push_back(ireq);
            MPI_Irecv(&next_proc, sizeof(Neighbor_Process), MPI_CHAR, mi.rank+1, 0, MPI_COMM_WORLD, &ireq);
            mpi_recvlist.push_back(ireq);
        }
        for (auto& ireq : mpi_sendlist) {
            MPI_Wait(&ireq, &st);
        }
        for (auto& ireq : mpi_recvlist) {
            MPI_Wait(&ireq, &st);
        }
        
        double dx;
        double left_limit;
        double right_limit;
        if (0<mi.rank) {
            dx = -1.0*std::pow(alpha*(previous_proc.counts - my_proc.counts), 3);
            left_limit = -1.0*(previous_proc.right - previous_proc.left)/2.0;
            right_limit = (my_proc.right - my_proc.left)/2.0;
            if (dx < left_limit) {
                dx = left_limit;
            } else if (dx > right_limit) {
                dx = right_limit;
            }
            this->left += dx;
        }
        if (mi.rank<mi.procs-1) {
            dx = -1.0*std::pow(alpha*(my_proc.counts - next_proc.counts), 3);
            left_limit = -1.0*(my_proc.right - my_proc.left)/2.0;
            right_limit = (next_proc.right - next_proc.left)/2.0;
            if (dx < left_limit) {
                dx = left_limit;
            } else if (dx > right_limit) {
                dx = right_limit;
            }
            this->right += dx;
        }

        migration_atoms.clear();
        migration_atoms.resize(mi.procs);
        for (auto& a : vars->atoms) {
            double sbx = std::floor((a.z-sysp::z_min)/lpz)*mi.npy*sysp::xl + std::floor((a.y-sysp::y_min)/lpy)*sysp::xl + a.x-sysp::x_min;
            if (sbx < this->left) {
                migration_atoms.at(mi.rank-1).push_back(a);
            } else if(sbx > this->right) {
                migration_atoms.at(mi.rank+1).push_back(a);
            } else {
                migration_atoms.at(mi.rank).push_back(a);
            }
        }

        vars->atoms.clear();
        vars->atoms.resize(migration_atoms.at(mi.rank).size());
        std::copy(migration_atoms.at(mi.rank).begin(), migration_atoms.at(mi.rank).end(), vars->atoms.begin());
        migration_atoms.at(mi.rank).clear();
        MPI_Barrier(MPI_COMM_WORLD);
        migrate_atoms(migration_atoms, vars, mi);

        std::vector<unsigned long> counts(mi.procs);
        pn = vars->number_of_atoms();
        MPI_Allgather(&pn, 1, MPI_UNSIGNED_LONG, counts.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        unsigned long sum_counts = std::accumulate(counts.begin(), counts.end(), 0);
        assert(sum_counts==sysp::N);
        
        // early stop
        unsigned long max_count = *std::max_element(counts.begin(), counts.end());
        if (max_count <= ideal_count_max || max_count-ideal_count < 10) {
            // std::fprintf(stderr, "***early stop (iter #%d)***\n", s);
            break;
        }
    }
}

// ============================================
