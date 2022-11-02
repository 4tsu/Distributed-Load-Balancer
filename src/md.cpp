#include "md.hpp"

MD::MD(MPIinfo mi){
    vars = new Variables();
    sysp = new Systemparam();
    obs = new Observer();
    sr = new SubRegion();
    pl = new PairList();
    this->mi = mi;
}

MD::~MD(void){
    delete vars;
    delete obs;
    delete sysp;
    delete sr;
    delete pl;
}

// -----------------------------------------------------

// MDクラスメンバ関数
void MD::set_params(int steps, int ob_interval, double dt) {
    this->steps = steps;
    this->ob_interval = ob_interval;
    this->dt = dt;
}



void MD::set_box(unsigned long N, double xl, double yl, double cutoff) {
    sysp->set_params(N, xl, yl, cutoff);
    if (N<std::numeric_limits<unsigned long>::min() || std::numeric_limits<unsigned long>::max()<N) {
        fprintf(stderr, "=== input 'N' is too large! ===\n");
        exit(EXIT_FAILURE);
    }
    sysp->calc_params();
}



void MD::set_margin(double margin) {
    sysp->margin = margin;
    sysp->calc_margin();
    vars->set_margin_life(margin);
}



void MD::set_sdd(int sdd_type) {
    this->sdd_type = sdd_type;
}



void MD::makeconf(void) {
    const unsigned long N = sysp->N;
    const double xl = sysp->xl;
    const double yl = sysp->yl;
    const double x_min = sysp->x_min;
    const double y_min = sysp->y_min;
    const double x_max = sysp->x_max;
    const double y_max = sysp->y_max;
    const double lpx = sysp->xl/mi.npx;
    const double lpy = sysp->yl/mi.npy;
    const unsigned long xppl = ceil(sqrt(xl*N/yl));
    const unsigned long yppl = ceil(sqrt(yl*N/xl));
    const double pitch = std::min(xl/xppl, yl/yppl);

   // 等間隔配置・分割
    for (unsigned long i=0; i<N; i++) {
        unsigned long iy = static_cast<unsigned long>(i/xppl);
        unsigned long ix = i%xppl;
        double x = ix * pitch;
        double y = iy * pitch;

        // どのプロセスに分配するかを判断する
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
    sr->make_dplist(mi, vars, sysp);

    // 他領域粒子情報をすべて持ってくる
    if (mi.procs > 1) {
        // あらかじめ、送受信のデータ容量だけやりとりしておく
        MPI_Request ireq;
        MPI_Status st;
        std::vector<MPI_Request> mpi_send_requests;
        unsigned long my_size = vars->number_of_atoms()*sizeof(Atom);
        for (auto &l : sr->dplist_reverse) {
            assert(l.i == mi.rank);
            MPI_Isend(&my_size, 1, MPI_UNSIGNED_LONG, l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);
        }

        std::vector<unsigned long> other_size_vec(sr->dplist.size());
        std::vector<MPI_Request> mpi_recv_requests;
        int proc_count = 0;
        for (auto &l : sr->dplist) {
            assert(l.i == mi.rank);
            MPI_Irecv(&other_size_vec[proc_count], 1, MPI_UNSIGNED_LONG, l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
            proc_count++;
        }
        for (auto &req : mpi_recv_requests)
            MPI_Wait(&req, &st);
        for (auto &req : mpi_send_requests)
            MPI_Wait(&req, &st);

        // まずは自分が送る分
        mpi_send_requests.clear();
        // 構造体をバラしてから送信する
        std::vector<Atom> sendbuf = vars->atoms;
        unsigned long send_size = sendbuf.size()*sizeof(Atom);
        for (auto &l : sr->dplist_reverse) {
            MPI_Isend(sendbuf.data(), send_size, MPI_CHAR, 
                    l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);
        }
        
        // そして受け取る分
        mpi_recv_requests.clear();
        unsigned long sum_recv = std::accumulate(other_size_vec.begin(), other_size_vec.end(), 0);
        std::vector<Atom> recvbuf(sum_recv/sizeof(Atom));
        unsigned long recv_position = 0;
        proc_count = 0;
        for (auto &l : sr->dplist) {
            MPI_Irecv(&recvbuf.at(recv_position), other_size_vec.at(proc_count), MPI_CHAR, 
                    l.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
            recv_position += other_size_vec[proc_count]/sizeof(Atom);
            proc_count++;
        }
        for (auto &req : mpi_recv_requests)
            MPI_Wait(&req, &st);
        for (auto &req : mpi_send_requests)
            MPI_Wait(&req, &st);

// for (auto a : recvbuf)
// fprintf(stderr, "%d ", a.id);
        vars->other_atoms.clear();
        unsigned long bias = 0;
        unsigned long recv_range;
        for (std::size_t i=0; i<sr->dplist.size(); i++) {
            recv_range = other_size_vec[i]/sizeof(Atom);
            std::vector<Atom> one_other_atoms(recv_range);
            std::copy(recvbuf.begin()+bias, recvbuf.begin()+bias+recv_range, one_other_atoms.begin());
            vars->other_atoms.push_back(one_other_atoms);
            bias += recv_range;
        }
    }



    // ローカルにペアリスト作成
    pl->make_pair(vars, sysp);
    // std::cerr << vars->other_atoms.size() << " == " << sr->dplist.size() << std::endl;



    // 今後通信すべき粒子リストを近くの領域間で共有しておく
    if (mi.procs > 1) {

        /// まず、リストのサイズrecv_sizeを送受信する。
        /// 同時にDomainPairListの整理を行う。
        /// 通信が必要なくなった場合はサイズ0を送り、これをもって互いのDomainPairListからこのペアは削除する
        MPI_Request ireq;
        MPI_Status st;
        std::vector<MPI_Request> mpi_send_requests;
        DomainPair *dplist = sr->dplist.data();
        std::vector<DomainPair> new_dplist;
        for (std::size_t i=0; i<sr->dplist.size(); i++) {
            assert(dplist[i].i == mi.rank);
            MPI_Isend(&vars->recv_size.at(i), 1, MPI_UNSIGNED_LONG, dplist[i].j, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);
            if (vars->recv_size.at(i) != 0)
                new_dplist.push_back(dplist[i]);
        }
        sr->dplist = new_dplist;

        std::vector<MPI_Request> mpi_recv_requests;
        vars->send_size.resize(sr->dplist_reverse.size());
        for (std::size_t i=0; i<sr->dplist_reverse.size(); i++) {
            DomainPair dp = sr->dplist_reverse.at(i);
            assert(dp.i == mi.rank);
            MPI_Irecv(&vars->send_size.at(i), 1, MPI_UNSIGNED_LONG, dp.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
        }
        
        for (auto& req : mpi_recv_requests)
            MPI_Wait(&req, &st);
        for (auto& req : mpi_send_requests)
            MPI_Wait(&req, &st);

        // 逆DomainPairの整理
        std::vector<DomainPair> new_dplist_r;
        for (std::size_t i=0; i<vars->send_size.size(); i++) {
            if (vars->send_size.at(i) != 0)
                new_dplist_r.push_back(sr->dplist_reverse.at(i));
        }
        sr->dplist_reverse = new_dplist_r;

        // send_sizeとrecv_size, recv_list, other_atomsの整理
        std::vector<unsigned long> new_recv_size;
        std::vector<std::vector<unsigned long>> new_recv_list;
        std::vector<std::vector<Atom>> new_other_atoms;
        std::vector<unsigned long> rs = vars->recv_size;
        std::vector<std::vector<unsigned long>> rl = vars->recv_list;
        std::vector<std::vector<Atom>> oa = vars->other_atoms;
        std::vector<std::vector<Pair>> ol = pl->other_list;
        for (std::size_t i=0; i<rs.size(); i++) {
            if (rs.at(i)!=0) {
                new_recv_size.push_back(rs.at(i));
                new_recv_list.push_back(rl.at(i));
                new_other_atoms.push_back(oa.at(i));
            }
        } 
        std::vector<unsigned long> new_send_size;
        for (auto& s : vars->send_size){
            if (s!=0)
                new_send_size.push_back(s);
        }
        vars->recv_size = new_recv_size;
        vars->recv_list = new_recv_list;
        vars->send_size = new_send_size;
        vars->other_atoms = new_other_atoms;
        assert(sr->dplist.size() == vars->recv_size.size());
        assert(sr->dplist.size() == vars->recv_list.size());
        assert(sr->dplist.size() == vars->other_atoms.size());
        assert(sr->dplist.size() == pl->other_list.size());
        assert(sr->dplist_reverse.size() == vars->send_size.size());
        

        // リストrecv_listの送信
        mpi_send_requests.clear();
        std::vector<std::vector<unsigned long>> sendbuf;
        for (std::size_t i=0; i<sr->dplist.size(); i++) {
            sendbuf.push_back(vars->recv_list.at(i));
            MPI_Isend(sendbuf.at(i).data(), sendbuf.at(i).size(), MPI_UNSIGNED_LONG, sr->dplist[i].j, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);
        }

        // リストrecv_listを受信して、send_listに格納
        mpi_recv_requests.clear();
        unsigned long send_list_total = std::accumulate(vars->send_size.begin(), vars->send_size.end(), 0) / sizeof(Atom);
        std::vector<unsigned long> recvbuf(send_list_total);
        unsigned long recv_index = 0;
        for (std::size_t i=0; i<sr->dplist_reverse.size(); i++) {
            unsigned long list_size = vars->send_size.at(i) / sizeof(Atom);
            MPI_Irecv(&recvbuf.at(recv_index), list_size, MPI_UNSIGNED_LONG, sr->dplist_reverse[i].j, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
            recv_index += list_size;
        }
        assert(recv_index == send_list_total);

        for (auto& req : mpi_recv_requests)
            MPI_Wait(&req, &st);
        for (auto& req : mpi_send_requests)
            MPI_Wait(&req, &st);

        /// recvbufのsend_listへの展開
        vars->send_list.clear();
        unsigned long bias = 0;
        for (std::size_t i=0; i<sr->dplist_reverse.size(); i++) {
            unsigned long recv_range = vars->send_size.at(i) / sizeof(Atom);
            std::vector<unsigned long> one_send_list(recv_range);
            std::copy(recvbuf.begin()+bias, recvbuf.begin()+bias+recv_range, one_send_list.begin());
            vars->send_list.push_back(one_send_list);
            bias += recv_range;
        }

        // send_listを受けて、send_atomsを詰めておく
        vars->pack_send_atoms();
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
        vars->set_margin_life(sysp->margin);
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
        assert(sysp->x_min <= x && x <= sysp->x_max);
        assert(sysp->y_min <= y && x <= sysp->y_max);
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
        atoms[pl.i].vx += df * dx;
        atoms[pl.i].vy += df * dy;
        atoms[pl.j].vx -= df * dx;
        atoms[pl.j].vy -= df * dy;
    }

    vars->sending_force.clear();
    // 自領域-他領域粒子ペア
    for (std::size_t i=0; i<pl->other_list.size(); i++) {

        // Atom *one_other_atoms = vars->other_atoms[i].data();
        std::vector<Atom> one_other_atoms = vars->other_atoms.at(i);
        std::vector<Force> one_sending_force;
        for (auto &pl : pl->other_list[i]) {
            Atom ia = atoms[pl.i];
            Atom ja = one_other_atoms.at(pl.j);
            assert(pl.idi == ia.id);
            assert(pl.idj == ja.id);
            double dx = ja.x - ia.x;
            double dy = ja.y - ia.y;
            periodic_distance(dx, dy, sysp);
            double r = sqrt(dx*dx + dy*dy);
            if (r > sysp->cutoff){
                continue;
            }
            
            double df = (24.0 * pow(r, 6) - 48.0) / pow(r, 14) * dt;
            atoms[pl.i].vx += df * dx;
            atoms[pl.i].vy += df * dy;
            Force sf;
            sf.id = ja.id;
            sf.vx =  - df * dx;
            sf.vy =  - df * dy;
            one_sending_force.push_back(sf);
        }
        vars->sending_force.push_back(one_sending_force);
    }
}



// 変化した座標を他領域と共有する
/// ペアリストにある粒子はすべてやりとりする
void MD::communicate_atoms(void) {
    // 必要な通信サイズを予め共有 <- ペアリスト構築時点で共有できるはず
    MPI_Request ireq;
    MPI_Status st;
    
    // 自領域粒子の情報を他領域に送る
    std::vector<MPI_Request> mpi_send_requests;
    std::vector<std::vector<Atom>> sendbuf(sr->dplist_reverse.size());
    for (std::size_t i=0; i<sr->dplist_reverse.size(); i++) {
        DomainPair dp = sr->dplist_reverse.at(i);
        assert(dp.i == mi.rank);
        sendbuf.at(i).resize(vars->send_atoms.at(i).size());
        for (std::size_t j=0; j<sendbuf.at(i).size(); j++) {
            sendbuf.at(i).at(j) = *vars->send_atoms.at(i).at(j);
        }
        MPI_Isend(sendbuf.at(i).data(), vars->send_size.at(i), MPI_CHAR, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_send_requests.push_back(ireq);
    }
    
    // 自領域の計算で使う他領域粒子の情報をもらう
    std::vector<MPI_Request> mpi_recv_requests;
    unsigned long total_recv_size = std::accumulate(vars->recv_size.begin(), vars->recv_size.end(), 0);
    std::vector<Atom> recvbuf(total_recv_size/sizeof(Atom));
    unsigned long recv_index = 0;
    for (std::size_t i=0; i<sr->dplist.size(); i++) {
        DomainPair dp = sr->dplist[i];
        assert(dp.i == mi.rank);
        MPI_Irecv(&recvbuf.at(recv_index), vars->recv_size.at(i), MPI_CHAR, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_recv_requests.push_back(ireq);
        recv_index += vars->recv_size[i]/sizeof(Atom);
    }

    for (auto& req : mpi_recv_requests)
        MPI_Wait(&req, &st);
    for (auto& req : mpi_send_requests)
        MPI_Wait(&req, &st);

    // 展開
    recv_index = 0;
    unsigned long range;
    for (std::size_t i=0; i<sr->dplist.size(); i++) {
        vars->other_atoms.at(i).clear();
        range = static_cast<unsigned long>(vars->recv_size.at(i) / sizeof(Atom));
        vars->other_atoms.at(i).resize(range);
        std::copy(recvbuf.begin()+recv_index, recvbuf.begin()+recv_index+range, vars->other_atoms.at(i).begin());
        recv_index += range;
    }
}



// 他領域粒子に関する力積を書き戻すための通信
/// ペアリストにある粒子全部 or 実際に相互作用した粒子のみ？
/// 事前に、粒子リストの行先とそのサイズ、送信元とそのサイズがわかっていればよい。
void MD::communicate_force(void) {
    MPI_Request ireq;
    MPI_Status st;
    // まず、sending_forceのサイズを共有する
    std::vector<MPI_Request> mpi_send_requests;

    for (std::size_t i=0; i<sr->dplist.size(); i++) {
        DomainPair dp = sr->dplist[i];
        assert(dp.i == mi.rank);
        unsigned long sending_force_size = vars->sending_force.at(i).size()*sizeof(Force);
        MPI_Isend(&sending_force_size, 1, MPI_UNSIGNED_LONG, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_send_requests.push_back(ireq);
    }
   
    std::vector<MPI_Request> mpi_recv_requests;
    std::vector<unsigned long> recv_force_size(sr->dplist_reverse.size());
    for (std::size_t i=0; i<sr->dplist_reverse.size(); i++) {
        DomainPair dp = sr->dplist_reverse[i];
        assert(dp.i == mi.rank);
        MPI_Irecv(&recv_force_size[i], 1, MPI_UNSIGNED_LONG, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_recv_requests.push_back(ireq);
    }

    for (auto& req : mpi_recv_requests)
        MPI_Wait(&req, &st);
    for (auto& req : mpi_send_requests)
        MPI_Wait(&req, &st);

    // sending_force自体の通信
    /// 送信
    mpi_send_requests.clear();
    for (std::size_t i=0; i<sr->dplist.size(); i++) {
        DomainPair dp = sr->dplist[i];
        unsigned long force_size = vars->sending_force.at(i).size()*sizeof(Force);
        if (force_size == 0)
            continue;
        MPI_Isend(vars->sending_force[i].data(), force_size, MPI_CHAR, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_send_requests.push_back(ireq);
    }

    /// 受信
    mpi_recv_requests.clear();
    unsigned long total_recv_size = static_cast<unsigned long>(std::accumulate(recv_force_size.begin(), recv_force_size.end(), 0));
    std::vector<Force> recv_force(total_recv_size/sizeof(Force));
    unsigned long recv_index = 0;
    for (std::size_t i=0; i<sr->dplist_reverse.size(); i++) {
        DomainPair dp = sr->dplist_reverse[i];
        if (recv_force_size[i] == 0)
            continue;
        MPI_Irecv(&recv_force[recv_index], recv_force_size[i], MPI_CHAR, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_recv_requests.push_back(ireq);
        recv_index += recv_force_size[i]/sizeof(Force);
    }

    for (auto& req : mpi_recv_requests)
        MPI_Wait(&req, &st);
    for (auto& req : mpi_send_requests)
        MPI_Wait(&req, &st);
    
   // 力の書き戻し
    if (recv_force.size() != 0) {
        unsigned long f_index = 0;
        Force f = recv_force.at(0);
        bool flag = false;
        for (std::size_t i=0; i<sr->dplist_reverse.size(); i++){
            for (Atom& a : vars->atoms){
                while (a.id==f.id) {
                    a.vx += f.vx;
                    a.vy += f.vy;
                    f_index++;
                    if (f_index >= recv_force.size()){
                        flag = true;
                        break;
                    }
                    f = recv_force.at(f_index);
                }
                if (flag) break;
            }
            if (flag) break;
        }
        assert(f_index == recv_force.size());
    }
}



// ----------------------------------------------------------------------



void MD::run(void) {
#ifdef FS
    // std::filesystemは使用できない環境もある．コンパイル時に選択可
    // 結果出力が追記なので、同名ファイルは事前に削除しておく
    if (mi.rank == 0) {
        for (const auto & file : std::filesystem::directory_iterator("./")) {
            std::string path = file.path();
            size_t word_pos = path.find("/cdv");
            if (word_pos != std::string::npos) {
                std::filesystem::remove_all(path);
                continue;
            } else if (path == "./energy.dat") {
                std::filesystem::remove(path);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    
    if (mi.procs<std::numeric_limits<int>::min() || std::numeric_limits<int>::max()<mi.procs) {
        fprintf(stderr, "=== Too many processes! ===\n");
        exit(EXIT_FAILURE);
    }
    
    /// MD
    // 初期配置orデータ読み込み
    makeconf();
   
     // ロードバランサー選択
    vars->set_initial_velocity(1.0, mi, sysp); // 初速決定
    obs->export_cdview(vars->atoms, *sysp, mi);

    //最初のペアリスト作成
    assert(sysp->N != 0);
    this->make_pair();
    /*
    for (auto dp : sr->dplist) {
        fprintf(stderr, "# %d %d-%d\n", mi.rank, dp.i, dp.j);
    }
    */
    // fprintf(stderr, "# %d %ld\n", mi.rank, sr->dplist.size());
    // fprintf(stderr, "# %d members=%ld\n", mi.rank, vars->atoms.size());
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


    // 計算ループ
    for (int step=1; step<=steps; step++) {
        if (mi.rank==0 && step%ob_interval==0) fprintf(stderr, "step %d\n", step);
        vars->time += dt;
        // シンプレクティック積分
        this->update_position(0.5);
        this->communicate_atoms();
/*MPI_Barrier(MPI_COMM_WORLD);
for (int p=0; p<mi.procs; p++) {
for (int i=0; i<vars->other_atoms.size(); i++){
for (auto atom : vars->other_atoms.at(i)) {
    fprintf(stderr, "%d\n", atom.id);
}fprintf(stderr, "\n\n");
}sleep(1);}*/
        this->calculate_force();
        this->communicate_force();
        this->update_position(0.5);
        this->communicate_atoms();
        // if (mi.rank==0) std::cerr << "(" << vars->atoms[0].x << vars->atoms[0].y << ")" << std::endl;
        // 情報の出力
        double k = obs->kinetic_energy(vars, sysp);
        double v = obs->potential_energy(vars, pl, sysp);
        if (mi.rank==0) {
            export_three("energy.dat", step, k, v, k+v);
        }
        if (step % ob_interval == 0) {
            obs->export_cdview(vars->atoms, *sysp, mi);
        }
        this->check_pairlist();
    }
}

// =====================================================