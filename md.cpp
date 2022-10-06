#include "md.hpp"
#include "systemparam.hpp"
#include "observer.hpp"
#include "variables.hpp"

// =====================================================

MD::MD(MPIinfo mi){
    vars = new Variables();
    sysp = new Systemparam();
    obs = new Observer();
    this->mi = mi;
}

MD::~MD(void){
    delete vars;
    delete obs;
    delete sysp;
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

    // 等間隔分割
    for (int i=0; i<N; i++) {
        static int my_id = 0;
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
            int id = myN*mi.rank + my_id;
            vars->add_atoms(id,x,y);
            assert(x_min<=x && x<=x_max);
            assert(y_min<=y && y<=y_max);
            my_id++;
        }
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
    
    // step 0 情報の出力
    obs->export_cdview(vars->atoms, *sysp, mi);

    // 計算ループ
    for (int step=1; step<=steps; step++) {
        vars->time += dt;
    }
}



// =====================================================