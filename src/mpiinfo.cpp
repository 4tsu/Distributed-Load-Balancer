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
    
    // int d2[2] = {};
    // MPI_Dims_create(procs, 2, d2);
    // mi.npx = d2[0];
    // mi.npy = d2[1];

    int max_divisor;
    for (int i=1; i<=std::floor(sqrt(procs)); i++) {
        if (procs%i==0) {
            max_divisor = i;
        }
    }
    mi.npx = max_divisor;
    mi.npy = procs/max_divisor;
}
// ==================================