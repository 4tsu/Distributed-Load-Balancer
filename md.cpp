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



    // 今後通信すべき粒子リストを近くの領域間で共有しておく
    if (mi.procs > 1) {

        /// まず、リストのサイズrecv_sizeを送受信する。
        /// 同時にDomainPairListの整理を行う。
        /// 通信が必要なくなった場合はサイズ0を送り、これをもって互いのDomainPairListからこのペアは削除する
        MPI_Request ireq;
        MPI_Status st;
        std::vector<MPI_Request> mpi_send_requests;
        DomainPair *dplist = dpl->dplist.data();
        std::vector<DomainPair> new_dplist;
        for (int i=0; i<dpl->dplist.size(); i++) {
            assert(dplist[i].i == mi.rank);
            int one_recv_size = vars->recv_size.at(i);
            MPI_Isend(&one_recv_size, 1, MPI_INT, dplist[i].j, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);
            if (one_recv_size != 0)
                new_dplist.push_back(dplist[i]);

        }
        dpl->dplist = new_dplist;



        std::vector<MPI_Request> mpi_recv_requests;
        vars->send_size.resize(dpl->dplist_reverse.size());
        for (int i=0; i<dpl->dplist_reverse.size(); i++) {
            DomainPair dp = dpl->dplist_reverse.at(i);
            assert(dp.i == mi.rank);
            MPI_Irecv(&vars->send_size.at(i), 1, MPI_INT, dp.j, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
        }
        
        for (auto& req : mpi_recv_requests) {
            MPI_Wait(&req, &st);
        }
        for (auto& req : mpi_send_requests) {
            MPI_Wait(&req, &st);
        }

        // 逆DomainPairの整理
        std::vector<DomainPair> new_dplist_r;
        for (int i=0; i<vars->send_size.size(); i++) {
            if (vars->send_size.at(i) != 0)
                new_dplist_r.push_back(dpl->dplist_reverse.at(i));
        }
        dpl->dplist_reverse = new_dplist_r;

        // send_sizeとrecv_sizeの整理
        std::vector<int> new_recv_size;
        for (auto& s : vars->recv_size){
            if (s!=0)
                new_recv_size.push_back(s);
        } 
        std::vector<int> new_send_size;
        for (auto& s : vars->send_size){
            if (s!=0)
                new_recv_size.push_back(s);
        } 
        
        // リストrecv_listの送信
        mpi_send_requests.clear();
        for (int i=0; i<dpl->dplist.size(); i++) {
            std::vector<int> one_recv_list = vars->recv_list.at(i);
            int list_size = vars->recv_size.at(i) / sizeof(Atom);
            MPI_Isend(one_recv_list.data(), one_recv_list.size(), MPI_INT, dpl->dplist[i].j, 0, MPI_COMM_WORLD, &ireq);
            mpi_send_requests.push_back(ireq);

sleep(0.5);
if (mi.rank==0 && dplist[i].j==1) {fprintf(stderr, "\nsend buf check ");
for(auto& e : vars->recv_list.at(i)) {
fprintf(stderr, "%d ", e);
}fprintf(stderr, "\n");}
        }
sleep(0.5);
std::cerr << mi.rank << " c.1" << std::endl;

        // リストrecv_listの受信
        mpi_recv_requests.clear();
        int send_list_total = std::accumulate(vars->send_size.begin(), vars->send_size.end(), 0) / sizeof(Atom);
        std::vector<int> recvbuf(send_list_total);
        int recv_index = 0;
        for (int i=0; i<dpl->dplist_reverse.size(); i++) {
            int list_size = vars->send_size.at(i) / sizeof(Atom);
            MPI_Irecv(&recvbuf.at(recv_index), list_size, MPI_INT, dpl->dplist_reverse[i].j, 0, MPI_COMM_WORLD, &ireq);
            mpi_recv_requests.push_back(ireq);
            recv_index += list_size;
        }

std::cerr << mi.rank << " c.2" << std::endl;

       
        for (auto& req : mpi_recv_requests) {
            MPI_Wait(&req, &st);
        }
        for (auto& req : mpi_send_requests) {
            MPI_Wait(&req, &st);
        }

std::cerr << mi.rank << " c.3" << std::endl;
sleep(0.5);
if (mi.rank==1) {fprintf(stderr, "\n recv buf check ");
for(int e : recvbuf) {
fprintf(stderr, "%d ", e);
}fprintf(stderr, "\n");}
sleep(0.5);

        vars->send_list.clear();
        int bias = 0;
        for (int i=0; i<dpl->dplist_reverse.size(); i++) {
            int recv_range = vars->send_size.at(i) / sizeof(Atom);
            std::vector<int> one_send_list(recv_range);
            std::copy(recvbuf.begin()+bias, recvbuf.begin()+bias+recv_range-1, one_send_list.begin());
            vars->send_list.push_back(one_send_list);
            bias += recv_range;
        }

fprintf(stderr, "send_list.at(0).at(0) %d\n", vars->send_list.at(0).at(0));
std::cerr << mi.rank << " c.4" << std::endl;

        // send_listを受けて、send_atomsを詰めておく
        vars->pack_send_atoms();

std::cerr << mi.rank << " c.5" << std::endl;
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
            Force sf;
            sf.id = ja.id;
            assert(pl.idi == ia.id);
            assert(pl.idj == ja.id);
            double dx = ja.x - ia.x;
            double dy = ja.y - ia.y;
            periodic_distance(dx, dy, sysp);
            double r = sqrt(dx*dx + dy*dy);
            if (r > sysp->cutoff){
                sf.vx = 0;
                sf.vy = 0;
                continue;
            }
            
            double df = (24.0 * pow(r, 6) - 48.0) / pow(r, 14) * dt;
            atoms[pl.i].x += df * dx;
            atoms[pl.i].y += df * dy;
            sf.vx = ja.vx - df * dx;
            sf.vy = ja.vy - df * dy;
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
    std::vector<MPI_Request> mpi_send_requests;
    // 自領域粒子の情報を他領域に送る
    for (int i=0; i<dpl->dplist_reverse.size(); i++) {
        DomainPair dp = dpl->dplist_reverse.at(i);
        assert(dp.i == mi.rank);
        std::vector<Atom> sendbuf(vars->send_atoms.size());
        for (int j=0; j<sendbuf.size(); j++)
            sendbuf.at(j) = *vars->send_atoms.at(i).at(j);
        MPI_Isend(&sendbuf, vars->send_size.at(i), MPI_CHAR, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_send_requests.push_back(ireq);
    }
    
    // 自領域の計算で使う他領域粒子の情報をもらう
    std::vector<MPI_Request> mpi_recv_requests;
    int total_recv_size = std::accumulate(vars->recv_size.begin(), vars->recv_size.end(), 0);
    std::vector<Atom> recv_atoms(static_cast<int>(total_recv_size)/sizeof(Atom));
    int recv_index = 0;

    for (int i=0; i<dpl->dplist.size(); i++) {
        DomainPair dp = dpl->dplist[i];
        assert(dp.i == mi.rank);
        std::vector<Atom> one_recv_atoms;
        MPI_Irecv(&recv_atoms[recv_index], vars->recv_size[i], MPI_CHAR, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_recv_requests.push_back(ireq);
        recv_index += static_cast<int>(vars->recv_size[i] / sizeof(Atom));
    }
    for (auto& req : mpi_recv_requests) {
        MPI_Wait(&req, &st);
    }
    for (auto& req : mpi_send_requests) {
        MPI_Wait(&req, &st);
    }
    recv_index = 0;
    int range;
    for (int i=0; i<dpl->dplist.size(); i++) {
        range = static_cast<int>(vars->recv_size[i] / sizeof(Atom));
        std::copy(recv_atoms.begin()+recv_index, recv_atoms.begin()+recv_index+range-1, vars->other_atoms[i].begin());
        assert(range == vars->other_atoms[i].size());
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
    for (int i=0; i<dpl->dplist.size(); i++) {
        DomainPair dp = dpl->dplist[i];
        assert(dp.i == mi.rank);
        int sending_force_size = sizeof(vars->sending_force[i]);
        MPI_Isend(&sending_force_size, 1, MPI_INT, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_send_requests.push_back(ireq);
    }
    std::vector<MPI_Request> mpi_recv_requests;
    std::vector<int> recv_force_size(dpl->dplist_reverse.size());
    for (int i=0; i<dpl->dplist_reverse.size(); i++) {
        DomainPair dp = dpl->dplist[i];
        assert(dp.i == mi.rank);
        MPI_Irecv(&recv_force_size[i], 1, MPI_INT, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_recv_requests.push_back(ireq);
    }
    for (auto& req : mpi_recv_requests) {
        MPI_Wait(&req, &st);
    }
    for (auto& req : mpi_send_requests) {
        MPI_Wait(&req, &st);
    }

    // sending_force自体の通信
    /// 送信
    mpi_send_requests.clear();
    for (int i=0; i<dpl->dplist.size(); i++) {
        DomainPair dp = dpl->dplist[i];
        Force *sending_force = vars->sending_force[i].data();
        int force_size = sizeof(vars->sending_force[i]);
        MPI_Isend(sending_force, force_size, MPI_CHAR, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_send_requests.push_back(ireq);
    }
    /// 受信
    mpi_recv_requests.clear();
    int total_recv_size = static_cast<int>(std::accumulate(recv_force_size.begin(), recv_force_size.end(), 0));
    std::vector<Force> recv_force(total_recv_size/sizeof(Force));
    int recv_index = 0;
    for (int i=0; i<dpl->dplist_reverse.size(); i++) {
        DomainPair dp = dpl->dplist[i];
        MPI_Irecv(&recv_force[recv_index], recv_force_size[i], MPI_CHAR, dp.j, 0, MPI_COMM_WORLD, &ireq);
        mpi_recv_requests.push_back(ireq);
        recv_index += recv_force_size[i]/sizeof(Force);
    }
    for (auto& req : mpi_recv_requests) {
        MPI_Wait(&req, &st);
    }
    for (auto& req : mpi_send_requests) {
        MPI_Wait(&req, &st);
    }

    // 力の書き戻し
    int f_index = 0;
    for (auto& atom : vars->atoms) {
        Force f = recv_force.at(f_index);
        if (atom.id == f.id) {
            atom.vx += f.vx;
            atom.vy += f.vy;
            f_index++;
        }
    assert(f_index == recv_force.size());
    }
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