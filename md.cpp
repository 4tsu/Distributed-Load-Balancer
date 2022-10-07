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
        for (int i=0; i<sum_recv; i++) {
            Atom a;
            a.id = recvbuf_id[i];
            a.x = recvbuf[i*4+0];
            a.y = recvbuf[i*4+1];
            a.vx = recvbuf[i*4+2];
            a.vy = recvbuf[i*4+3];
            vars->other_atoms.push_back(a);
        }
    }
    // ローカルにペアリスト作成
    pl->make_pair(vars, sysp);
}



// ペアリスト作成・更新
void MD::check_pairlist(void) {

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
    // std::cout << mi.rank << " pairlist " << pl->list.size() << std::endl;
    // std::cout << mi.rank << " pairlist other" << pl->other_list.size() << std::endl;
    // std::cout << mi.rank << " other_atoms " << vars->other_atoms.size() << std::endl;
    for (auto &pl : pl->list) {
        std::cout << pl.idi << " " << pl.idj << std::endl;
    }
    for (auto &pl : pl->other_list) {
        std::cout << std::min(pl.idi, pl.idj) << " " << std::max(pl.idi, pl.idj) << std::endl;
    }
    // step 0 情報の出力
    obs->export_cdview(vars->atoms, *sysp, mi);

    // 計算ループ
    for (int step=1; step<=steps; step++) {
        vars->time += dt;
    }
}



// =====================================================