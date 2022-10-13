#include "md.hpp"

// ========================================

int main(int argc, char **argv) {
    // 実行時間計測
    clock_t start = clock();

    MPI_Init(&argc, &argv);
    MPIinfo mi;
    setup_info(mi);

    MD *md = new MD(mi);
    md->set_params(1000, 1, 0.0020); // int steps, int ob_interval, double dt
    md->set_box(100, 10, 10, 3.5); // int N, double xl, double yl, double cutoff
    md->set_margin(0.5);
    md->run();

    MPI_Finalize();

    // 実行時間計測
    delete md;
    clock_t end = clock();
    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC*1000.0;
    fprintf(stderr, "time %lf[ms]\n", time);
}
// ========================================