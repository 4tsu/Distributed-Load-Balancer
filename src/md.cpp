#include "md.hpp"

namespace sysp = systemparam;

MD::MD(MPIinfo mi){
    vars = new Variables();
    obs = new Observer();
    sr = new SubRegion();
    pl = new PairList();
    calctimer = new CalcTimer();
    grosstimer = new CalcTimer();
    commtimer = new CalcTimer();
    sddtimer = new CalcTimer();
    wholetimer = new CalcTimer();
    this->mi = mi;
}

MD::~MD(void){
    delete vars;
    delete obs;
    delete sr;
    delete pl;
    delete sdd;
    delete calctimer;
    delete grosstimer;
    delete commtimer;
    delete sddtimer;
    delete wholetimer;
}

// -----------------------------------------------------

// MDクラスメンバ関数
void MD::set_params(int steps, int ob_interval, double dt) {
    this->steps = steps;
    this->ob_interval = ob_interval;
    this->dt = dt;
}



void MD::set_box(unsigned long N, double xl, double yl) {
    sysp::N = N;
    sysp::xl = xl;
    sysp::yl = yl;
    if (N<std::numeric_limits<unsigned long>::min() || std::numeric_limits<unsigned long>::max()<N) {
        std::fprintf(stderr, "=== input 'N' is too large! ===\n");
        exit(EXIT_FAILURE);
    }
}



void MD::set_cutoff(double cutoff) {
    sysp::cutoff = cutoff;
}



void MD::set_margin(double margin) {
    sysp::margin = margin;
    sysp::calc_margin();
    vars->set_margin_life(margin);
}



void MD::set_sdd(int sdd_type) {
    sdd = new Sdd(sdd_type);
}



void MD::set_config(const std::string conf) {
    this->config = conf;
}



void MD::makeconf(void) {
    const unsigned long N = sysp::N;
    const double xl = sysp::xl;
    const double yl = sysp::yl;
    const double x_min = sysp::x_min;
    const double y_min = sysp::y_min;
    const double x_max = sysp::x_max;
    const double y_max = sysp::y_max;
    const double lpx = sysp::xl/mi.npx;
    const double lpy = sysp::yl/mi.npy;
    const unsigned long xppl = ceil(std::sqrt(xl*N/yl));
    const unsigned long yppl = ceil(std::sqrt(yl*N/xl));
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
    sr->make_dplist(mi, vars);

    // 領域内でペアリスト作成
    // 粒子の空間ソートも同時に行うので、通信の前に行うのが良い
    pl->make_pair(vars);

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



    // ローカルに領域間ペアリスト作成
    pl->make_pair_ext(vars);



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
        sddtimer->start();
        sdd->run(vars, mi, sr);   // 領域分割の切り直しも同時に行う
        sddtimer->stop();
        this->make_pair();
        vars->set_margin_life(sysp::margin);
    }
    double ml = vars->margin_life;
    MPI_Bcast(&ml, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    vars->margin_life = ml;
}



void MD::update_position(double coefficient) {
    calctimer->start();
    for (auto& atom : vars->atoms) {
        if (atom.vx>9 || atom.vy>9) {
            std::fprintf(stderr, "Abnormal velocity!(rank#%d atom#%ld:[%lf, %lf])",mi.rank, atom.id, atom.vx, atom.vy);
            abort();
        }
        double x = atom.x + atom.vx * dt * coefficient;
        double y = atom.y + atom.vy * dt * coefficient;
        periodic_coordinate(x, y);
        atom.x = x;
        atom.y = y;
        assert(sysp::x_min <= x && x <= sysp::x_max);
        assert(sysp::y_min <= y && y <= sysp::y_max);
    }
    calctimer->stop();
}



void MD::calculate_force(void) {
    calctimer->start();
    // 自領域内粒子同士
    Atom *atoms = vars->atoms.data();
    for (auto &pair : pl->list) {
        Atom ia = atoms[pair.i];
        Atom ja = atoms[pair.j];
        assert(pair.idi == ia.id);
        assert(pair.idj == ja.id);
        double dx = ja.x - ia.x;
        double dy = ja.y - ia.y;
        periodic_distance(dx, dy);
        double r = std::sqrt(dx*dx + dy*dy);
        double df = 0.0;
        if (r <= sysp::cutoff) {
            df = (24.0 * pow(r, 6) - 48.0) / pow(r, 14) * dt;
        }
        // if ((df*df)>1.5) {
        //    std::fprintf(stderr, "Abnormal Force! (rank#%d pair[%ld-%ld])\n", mi.rank, ia.id, ja.id);
        //    abort();
        // }
        atoms[pair.i].vx += df * dx;
        atoms[pair.i].vy += df * dy;
        atoms[pair.j].vx -= df * dx;
        atoms[pair.j].vy -= df * dy;
    }

    vars->sending_force.clear();
    // 自領域-他領域粒子ペア
    for (std::size_t i=0; i<pl->other_list.size(); i++) {

        // Atom *one_other_atoms = vars->other_atoms[i].data();
        std::vector<Atom> one_other_atoms = vars->other_atoms.at(i);
        std::vector<Force> one_sending_force(pl->other_list[i].size());
        for (std::size_t j=0; j<pl->other_list[i].size(); j++) {
            Pair pair = pl->other_list[i][j];
            Atom ia = atoms[pair.i];
            Atom ja = one_other_atoms[pair.j];
            assert(pair.idi == ia.id);
            assert(pair.idj == ja.id);
            double dx = ja.x - ia.x;
            double dy = ja.y - ia.y;
            periodic_distance(dx, dy);
            double r = std::sqrt(dx*dx + dy*dy);
            Force sf;
            sf.id = ja.id;
            sf.vx = 0.0;
            sf.vy = 0.0;
            if (r <= sysp::cutoff){
                double df = (24.0 * pow(r, 6) - 48.0) / pow(r, 14) * dt;
                atoms[pair.i].vx += df * dx;
                atoms[pair.i].vy += df * dy;
                sf.vx =  - df * dx;
                sf.vy =  - df * dy;
            }
           one_sending_force[j] = sf;
        }
        vars->sending_force.push_back(one_sending_force);
    }
    calctimer->stop();
}



// 変化した座標を他領域と共有する
/// ペアリストにある粒子はすべてやりとりする
void MD::communicate_atoms(void) {
    commtimer->start();
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
    commtimer->stop();
}



// 他領域粒子に関する力積を書き戻すための通信
/// ペアリストにある粒子全部 or 実際に相互作用した粒子のみ？
/// 事前に、粒子リストの行先とそのサイズ、送信元とそのサイズがわかっていればよい。
void MD::communicate_force(void) {
    commtimer->start();
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
    commtimer->stop();
}



void MD::read_data(const std::string filename, Variables* vars, const MPIinfo &mi) {
    std::ifstream reading_file;
    reading_file.open(filename, std::ios::in);
    std::string line;
    bool is_step = false;
    bool is_init = false;
    bool is_num = false;
    int is_bounds = 0;
    double lpx, lpy;
    unsigned long id = 0;
    while(std::getline(reading_file, line)) {
        // LAMMPS出力ではじめに出てくる初期配置のdumpであれば無視する。
        if (line=="ITEM: TIMESTEP") {
            is_step = true;
            is_init = false;
            continue;
        } else if (is_step) {
            is_step = false;
            if (line=="0") {
                is_init = true;
            } else {
                begin_step = std::stoi(line);
            }
            continue;
        }
        if (is_init) {
            continue;
        }
        
        // ユーザーが定義したdump、もしくはLAMMPSの最終結果のdumpを読み込む
        if (line=="ITEM: NUMBER OF ATOMS") {
            is_num = true;
            continue;
        } else if (is_num) {
            is_num = false;
            sysp::N = std::stoul(line);
            continue;
        }

        if (std::equal(line.begin(), line.begin()+16, "ITEM: BOX BOUNDS")) {
            is_bounds = true;
            continue;
        } else if (std::equal(line.begin(), line.begin()+10, "ITEM: ATOM")) {
            is_bounds = 0;
            continue;
        } else if (is_bounds) {
            std::string var1, var2;
            auto space_pos = line.find(" ", 0);
            auto length = line.size();
            for (std::size_t i=0; i<space_pos; i++)
                var1 += line[i];
            for (std::size_t i=space_pos; i<length; i++)
                var2 += line[i];
            if (is_bounds==1) {
                sysp::x_min = std::stod(var1);
                sysp::x_max = std::stod(var2);
                sysp::xl = sysp::x_max - sysp::x_min;
                lpx = sysp::xl/mi.npx;
            } else if (is_bounds==2) {
                sysp::y_min = std::stod(var1);
                sysp::y_max = std::stod(var2);
                sysp::yl = sysp::y_max - sysp::y_min;
                lpy = sysp::yl/mi.npy;
            }
            is_bounds++;
            continue;
        }

        double x, y, vx, vy;
        std::vector<std::string> var;
        auto offset = std::string::size_type(0);
        while (true) {
            auto pos = line.find(" ", offset);
            if (pos==std::string::npos) {
                var.push_back(line.substr(offset));
                break;
            }
            var.push_back(line.substr(offset, pos-offset));
            offset = pos + 1;
        }
        x  = std::stod(var.at(0));
        y  = std::stod(var.at(1));
        vx = std::stod(var.at(3));
        vy = std::stod(var.at(4));
        periodic_coordinate(x, y);
        int ip = static_cast<int>(floor((y-sysp::y_min)/lpy)*mi.npx + floor((x-sysp::x_min)/lpx));
        if (ip==mi.rank) {
            Atom a;
            a.id = id;
            a.x  = x;
            a.y  = y;
            a.vx = vx;
            a.vy = vy;
            vars->atoms.push_back(a);
        }
        id++;
    }
}



void MD::get_exec_time(int step, CalcTimer* ct, const std::string filename) {
    std::vector<double> results(mi.procs);
    ct->share(results.data());

    if (mi.rank==0) {
        double max_time = *std::max_element(results.begin(), results.end());
        double min_time = *std::min_element(results.begin(), results.end());
        double avg_time = std::accumulate(results.begin(), results.end(), 0.0);
        avg_time /= static_cast<double>(mi.procs);
        export_three(filename, step, min_time, max_time, avg_time);
    }
}

// ----------------------------------------------------------------------



void MD::run(int trial) {
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
            } else if (path == "./time_gross.dat") {
                std::filesystem::remove(path);
            } else if (path == "./time_net.dat") {
                std::filesystem::remove(path);
            } else if (path == "./time_sdd.dat") {
                std::filesystem::remove(path);
            } else if (path == "./time_whole.dat") {
                std::filesystem::remove(path);
            } else if (path == "./time_comm.dat") {
                std::filesystem::remove(path);
            } else if (path == "./load_balance.dat") {
                std::filesystem::remove(path);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    
    if (mi.procs<std::numeric_limits<int>::min() || std::numeric_limits<int>::max()<mi.procs) {
        std::fprintf(stderr, "=== Too many processes! ===\n");
        exit(EXIT_FAILURE);
    }
    
    // 出力ファイル準備
    std::string net_time_out   = "time_net";
    std::string gross_time_out = "time_gross";
    std::string sdd_time_out   = "time_sdd";
    std::string whole_time_out = "time_whole";
    std::string comm_time_out  = "time_comm";
    if (trial > 0) {
        net_time_out   += "_" + std::to_string(trial);
        gross_time_out += "_" + std::to_string(trial);
        sdd_time_out   += "_" + std::to_string(trial);
        whole_time_out += "_" + std::to_string(trial);
        comm_time_out  += "_" + std::to_string(trial);
    }
    net_time_out   += ".dat";
    gross_time_out += ".dat";
    sdd_time_out   += ".dat";
    whole_time_out += ".dat";
    comm_time_out  += ".dat";

    /// MD
    // 初期配置orデータ読み込み
    if (this->config=="make") {
        sysp::calc_params();
        makeconf();
        vars->set_initial_velocity(1.0, mi); // 初速決定
    } else {
        this->read_data(config, vars, mi);
        sysp::calc_params();
    }

     // ロードバランサー選択
    sddtimer->start();
    sdd->init(vars, mi, sr);
    sdd->run(vars, mi, sr);
    sddtimer->stop();
    this->get_exec_time(0, sddtimer, sdd_time_out);
    sddtimer->reset();
    
    obs->export_cdview(vars, mi, std::ceil(begin_step/ob_interval));

    //最初のペアリスト作成
    assert(sysp::N != 0);
    this->make_pair();

    // 計算ループ
    for (int step=1+begin_step; step<=steps+begin_step; step++) {
        if (mi.rank==0 && step%ob_interval==0) std::fprintf(stderr, "step %d\n", step);
        wholetimer->start();
        vars->time += dt;

        grosstimer->start();
        // シンプレクティック積分
        this->update_position(0.5);
        this->communicate_atoms();
        this->calculate_force();
        this->communicate_force();
        this->update_position(0.5);
        this->communicate_atoms();
        grosstimer->stop();

        // 情報の出力
        double k = obs->kinetic_energy(vars);   // rank0が値を受け取る
        double v = obs->potential_energy(vars, pl);   // rank0が値を受け取る
        if (mi.rank==0) {
            export_three("energy.dat", step, k, v, k+v);
        }
        if (step % ob_interval == 0) {
            obs->export_cdview(vars, mi);
            if (step % (ob_interval*10) == 0) {
                obs->export_checkpoint("1.ckpt", step, vars, mi);
            } else if (step % (ob_interval*5) == 0) {
                obs->export_checkpoint("2.ckpt", step, vars, mi);
            }
        }

        this->get_exec_time(step, grosstimer, gross_time_out);
        this->get_exec_time(step, calctimer,  net_time_out);
        this->get_exec_time(step, commtimer,  comm_time_out);
        grosstimer->reset();
        calctimer->reset();
        commtimer->reset();
        obs->export_workload(step, vars, mi);

        this->check_pairlist();
        this->get_exec_time(step, sddtimer, sdd_time_out);
        sddtimer->reset();

        wholetimer->stop();
        this->get_exec_time(step, wholetimer, whole_time_out);
        wholetimer->reset();
    }
}

// =====================================================
