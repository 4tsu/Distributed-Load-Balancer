#include "md.hpp"

// ========================================

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPIinfo mi;
    setup_info(mi);

    // 実行時間計測
    clock_t start = clock();

    MD *md = new MD(mi);
    md->set_params(2000, 20, 0.0020); // int steps, int ob_interval, double dt
    md->set_box(1600, 40, 40); // unsigned long N, double xl, double yl
    md->set_cutoff(3.5); // double cutoff
    md->set_margin(0.5);
    md->set_config("make");
    md->set_sdd(0);
    md->run();
    delete md;

    // 実行時間計測
    clock_t end = clock();
    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC*1000.0;
    double avg_time;
    MPI_Allreduce(&time, &avg_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (mi.rank == 0)
        fprintf(stderr, "avg.time %lf[ms] (%d procs)\n", avg_time/static_cast<double>(mi.procs), mi.procs);

    MPI_Finalize();
}
// ========================================