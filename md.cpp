#include "md.hpp"

MD::MD(MPIinfo mi){
    vars = new Variables();
    sysp = new Systemparam();
    obs = new Observer();
    dpl = new DomainPairList();
    pl = new PairList();
    this->mi = mi;
}

MD::~MD(void){
    delete vars;
    delete obs;
    delete sysp;
    delete dpl;
    delete pl;
}

// -----------------------------------------------------

// MDクラスメンバ関数
void MD::set_params(int steps, int ob_interval, double dt) {
    this->steps = steps;
    this->ob_interval = ob_interval;
    this->dt = dt;
}



void MD::set_box(int N, double xl, double yl, double cutoff) {
    sysp->set_params(N, xl, yl, cutoff, mi.procs);
    sysp->calc_params();
}



void MD::set_margin(double margin) {
    sysp->margin = margin;
    sysp->calc_margin();
}



void MD::set_sdd(int sdd_type) {
    this->sdd_type = sdd_type;
}



void MD::makeconf(void) {
    int N = sysp->N;
    int myN = sysp->myN;
    double xl = sysp->xl;
    double yl = sysp->yl;
    double x_min = sysp->x_min;
    double y_min = sysp->y_min;
    double x_max = sysp->x_max;
    double y_max = sysp->y_max;
    int xppl = ceil(sqrt(xl*N/yl));
    int yppl = ceil(sqrt(yl*N/xl));
    double pitch = std::min(xl/xppl, yl/yppl);

    // 等間隔配置・分割
    for (int i=0; i<N; i++) {
        int iy = static_cast<int>(i/xppl);
        int ix = i%xppl;
        double x = ix * pitch;
        double y = iy * pitch;

        // どのプロセスに分配するかを判断する
        int lpx = sysp->xl/mi.npx;
        int lpy = sysp->yl/mi.npy;
        int ip = static_cast<int>(floor(y/lpy)*mi.npx + floor(x/lpx));
        if (ip==mi.rank) {
            x += x_min;
            y += y_min;
            vars->add_atoms(i,x,y);
            assert(x_min<=x && x<=x_max);
            assert(y_min<=y && y<=y_max);
        }
    }
}



// ペアリスト作成・更新
void MD::make_pair(void) {
    // centerとradiusを計算済みである事を仮定
    // domainpairlistの作成
    dpl->make_list(mi);
    // 他領域粒子情報をすべて持ってくる
    if (mi.procs > 1) {
        // あらかじめ、送受信のデータ容量だけやりとりしておく
        MPI_Request ireq;
        MPI_Status st;
        std::vector<MPI_Request> mpi_send_requests;
        for (auto &l : dpl->dplist_reverse) {
            assert(l.i == mi.rank);
            int my_n = vars->number_of_atoms();
            MPI_Isend(&my_n, 1, MPI_INT, l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);
        }

        std::vector<int> other_n_vec(dpl->dplist.size());
        std::vector<MPI_Request> mpi_recv_requests;
        int other_count = 0;
        for (auto &l : dpl->dplist) {
            assert(l.i == mi.rank);
            int other_n;
            MPI_Irecv(&other_n_vec[other_count], 1, MPI_INT, l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
            other_count++;
        }
        for (auto &req : mpi_recv_requests) {
            MPI_Wait(&req, &st);
        }
        for (auto &req : mpi_send_requests) {
            MPI_Wait(&req, &st);
        }

        // まずは自分が送る分
        mpi_send_requests.clear();
        // 構造体をバラしてから送信する
        std::vector<int> sendbuf_id(vars->number_of_atoms());
        std::vector<double> sendbuf(vars->number_of_atoms()*4);
        Atom *atoms = vars->atoms.data();
        for (int i=0; i<vars->number_of_atoms(); i++) {
            sendbuf_id[i] = atoms[i].id;
            sendbuf[i*4+0] = atoms[i].x;
            sendbuf[i*4+1] = atoms[i].y;
            sendbuf[i*4+2] = atoms[i].vx;
            sendbuf[i*4+3] = atoms[i].vy;
        }
        for (auto &l : dpl->dplist_reverse) {
            MPI_Isend(sendbuf_id.data(), sendbuf_id.size(), MPI_INT, 
                    l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);
            MPI_Isend(sendbuf.data(), sendbuf.size(), MPI_DOUBLE, 
                    l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);
        }
        
        // そして受け取る分
        mpi_recv_requests.clear();
        int sum_recv = std::accumulate(other_n_vec.begin(), other_n_vec.end(), 0);
        other_count = 0;
        int start_position = 0;
        std::vector<int> recvbuf_id(sum_recv);
        std::vector<double> recvbuf(sum_recv*4);
        for (auto &l : dpl->dplist) {
            MPI_Irecv(&recvbuf_id[start_position], other_n_vec[other_count], MPI_INT, 
                    l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
            MPI_Irecv(&recvbuf[start_position*4], other_n_vec[other_count]*4, MPI_DOUBLE, 
                    l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
            start_position += other_n_vec[other_count];
            other_count++;
        }
        
        for (auto &req : mpi_recv_requests) {
            MPI_Wait(&req, &st);
        }
        for (auto &req : mpi_send_requests) {
            MPI_Wait(&req, &st);
        }
        // 受け取ったデータを構造体に戻す
        vars->other_atoms.clear();
        other_count = 0;
        int bias = 0;
        for (auto &dp : dpl->dplist) {
            std::vector<Atom> one_other_atom;
            for (int i=0; i<other_n_vec[other_count]; i++) {
                Atom a;
                a.id = recvbuf_id[bias + i];
                a.x = recvbuf[(bias+i)*4+0];
                a.y = recvbuf[(bias+i)*4+1];
                a.vx = recvbuf[(bias+i)*4+2];
                a.vy = recvbuf[(bias+i)*4+3];
                one_other_atom.push_back(a);
            }
            bias += other_n_vec[other_count];
            vars->other_atoms.push_back(one_other_atom);
        }
    }
    // ローカルにペアリスト作成
    pl->make_pair(vars, sysp, dpl);
    // std::cerr << vars->other_atoms.size() << " == " << dpl->dplist.size() << std::endl;

    // 今後通信すべき粒子のリストを近くの領域間で共有する
    if (mi.procs > 1) {
    /// まず、リストのサイズを送受信する
    /// 通信が必要なくなった場合はサイズ0を送り、これをもって互いのDomainPairListからこのペアは削除する
        MPI_Request ireq;
        MPI_Status st;
        std::vector<MPI_Request> mpi_send_requests;
        /// 自領域がほしい粒子情報のリクエストであることに注意。
        for (auto &l : dpl->dplist) {
            assert(l.i == mi.rank);
            int want_n = vars->other_atoms.size();
        }
    }
}



// ペアリスト作成・更新
void MD::check_pairlist(void) {
    double local_max = vars->max_velocity();
    double global_max;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    vars->margin_life -= global_max*2.0*dt;
    if (vars->margin_life < 0) {
        this->make_pair();
        vars->margin_life = sysp->margin;
    }
    double ml = vars->margin_life;
    MPI_Bcast(&ml, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    vars->margin_life = ml;
}



void MD::update_position(double coefficient) {
    for (auto& atom : vars->atoms) {
        double x = atom.x + atom.vx * dt * coefficient;
        double y = atom.y + atom.vy * dt * coefficient;
        periodic_coordinate(x, y, sysp);
        atom.x = x;
        atom.y = y;
        assert(sysp->x_min < x < sysp->x_max);
        assert(sysp->y_min < y < sysp->y_max);
    }
}



void MD::calculate_force(void) {
    // 自領域内粒子同士
    Atom *atoms = vars->atoms.data();
    for (auto &pl : pl->list) {
        Atom ia = atoms[pl.i];
        Atom ja = atoms[pl.j];
        assert(pl.idi == ia.id);
        assert(pl.idj == ja.id);
        double dx = ja.x - ia.x;
        double dy = ja.y - ia.y;
        periodic_distance(dx, dy, sysp);
        double r = sqrt(dx*dx + dy*dy);
        if (r > sysp->cutoff)
            continue;
        
        double df = (24.0 * pow(r, 6) - 48.0) / pow(r, 14) * dt;
        atoms[pl.i].x += df * dx;
        atoms[pl.i].y += df * dy;
        atoms[pl.j].x -= df * dx;
        atoms[pl.j].y -= df * dy;
    }
    vars->sending_force.clear();
    // 自領域-他領域粒子ペア
    for (int i=0; i<pl->other_list.size(); i++) {
        Atom *one_other_atoms = vars->other_atoms[i].data();
        std::vector<Force> one_sending_force;
        for (auto &pl : pl->other_list[i]) {
            Atom ia = atoms[pl.i];
            Atom ja = one_other_atoms[pl.j];
            assert(pl.idi == ia.id);
            assert(pl.idj == ja.id);
            double dx = ja.x - ia.x;
            double dy = ja.y - ia.y;
            periodic_distance(dx, dy, sysp);
            double r = sqrt(dx*dx + dy*dy);
            if (r > sysp->cutoff)
                continue;
            
            double df = (24.0 * pow(r, 6) - 48.0) / pow(r, 14) * dt;
            atoms[pl.i].x += df * dx;
            atoms[pl.i].y += df * dy;
            Force sf;
            sf.id = ja.id;
            sf.vx = ja.vx - df * dx;
            sf.vy = ja.vy - df * dy;
            one_sending_force.push_back(sf);
        }
        vars->sending_force.push_back(one_sending_force);
    }
}



void MD::communicate_atoms(void) {

}



void MD::communicate_force(void) {

}



void MD::run(void) {
    // 結果出力が追記なので、同名ファイルは事前に削除しておく
    if (mi.rank == 0) {
        for (const auto & file : std::filesystem::directory_iterator(".")) {
            std::string path = file.path();
            int word_pos = path.find(".cdv");
            if (word_pos != std::string::npos) {
                std::filesystem::remove(path);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    
    
    /// MD
    // 初期配置orデータ読み込み
    makeconf();
    
     // ロードバランサー選択
    vars->set_initial_velocity(1.0, mi); // 初速決定
     //最初のペアリスト作成
    assert(sysp->N != 0);
    this->make_pair();
    // std::cout << mi.rank << " members " << vars->atoms.size() << std::endl;
    /*
    for (auto& one_other_list : pl->other_list) {
        std::cout << mi.rank << " pairlist other" << one_other_list.size() << std::endl;
    }
    */
    /*
    for (auto &one_other_atoms : vars->other_atoms) {
        std::cout << mi.rank << " other_atoms " << one_other_atoms.size() << std::endl;
    }
    */
    // std::cout << mi.rank << " pairlist " << pl->list.size() << std::endl;
    /*
    for (auto &pl : pl->list) {
        std::cout << pl.idi << " " << pl.idj << std::endl;
    }
    for (auto &opl : pl->other_list) {
        for (auto& pl : opl) {
            std::cout << std::min(pl.idi, pl.idj) << " " << std::max(pl.idi, pl.idj) << std::endl;
        }
    }
    */
    // step 0 情報の出力
    obs->export_cdview(vars->atoms, *sysp, mi);

    // 計算ループ
    for (int step=1; step<=steps; step++) {
        if (mi.rank==0) printf("step %d\n", step);
        vars->time += dt;
        // シンプレクティック積分
        this->update_position(0.5);
        // this->communicate_atoms();
        // this->calculate_force();
        // this->communicate_force();
        this->update_position(0.5);
        // this->communicate_atoms();
        // if (mi.rank==0) std::cerr << "(" << vars->atoms[0].x << vars->atoms[0].y << ")" << std::endl;
        // 情報の出力
        if (step % ob_interval == 0) {
            obs->export_cdview(vars->atoms, *sysp, mi);
        }
        this->check_pairlist();
    }
}

// =====================================================