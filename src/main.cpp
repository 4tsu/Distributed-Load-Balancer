#include "md.hpp"
// #include "2d/md.hpp"

// ========================================

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPIinfo mi;
    setup_info(mi);

    // 実行時間計測
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    MD *md = new MD(mi);
    md->set_params(20, 20, 0.0020); // int steps, int ob_interval, double dt
    md->set_box(8000, 20, 20, 20); // unsigned long N, double xl, double yl(, double zl)
    md->set_cutoff(3.5); // double cutoff
    md->set_margin(0.5); // double margin
    md->set_config("make"); // configuration file name
    md->set_sdd(0); // Load-Balancer type
    md->run();
    delete md;

    // 実行時間計測
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    const double time = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000.0 / 1000.0);
    double avg_time;
    MPI_Allreduce(&time, &avg_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (mi.rank == 0)
        std::fprintf(stderr, "avg.time %lf[ms] (%d procs)\n", avg_time/static_cast<double>(mi.procs), mi.procs);

    MPI_Finalize();
}
// ========================================