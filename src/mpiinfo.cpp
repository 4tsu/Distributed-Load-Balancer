#include "mpiinfo.hpp"

// -----------------------------------------------------

// MPIプロセス情報を管理する構造体のセットアップ
void setup_info(MPIinfo &mi){
    int rank = 0;
    int procs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    mi.rank = rank;
    mi.procs = procs;
    
    int max_divisor1, max_divisor2;
    for (int i=1; i<=std::floor(cbrt(procs)); i++) {
        if (procs%i==0) {
            max_divisor1 = i;
        }
    }
    for (int i=1; i<=std::floor(sqrt(procs/max_divisor1)); i++) {
        if ((procs/max_divisor1)%i==0) {
            max_divisor2 = i;
        }
    }
    mi.npx = max_divisor1;
    mi.npy = max_divisor2;
    mi.npz = procs/(max_divisor1*max_divisor2);
}
// ==================================