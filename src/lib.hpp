#pragma once

#include <cassert>
#include <cstdio>
#include <math.h>
#include <unistd.h>
#include <cassert>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <limits>
#include <chrono>

#include <mpi.h>

#ifdef FS
#include <filesystem>
#endif

// ----------------------

//デバッグ用
inline void ckpt() {
    static int c = 1;
    fprintf(stderr, "c.%d\n", c);
    MPI_Barrier(MPI_COMM_WORLD);
    c++;
}
inline void ckpt2() {
    static int d = 1;
    fprintf(stderr, "d.%d\n", d);
    MPI_Barrier(MPI_COMM_WORLD);
    d++;
}
inline void ckpt3() {
    static int e = 1;
    fprintf(stderr, "e.%d\n", e);
    MPI_Barrier(MPI_COMM_WORLD);
    e++;
}
inline int getrank() {
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    return my_rank;
}
inline void slprnk() {
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    sleep(my_rank*2);
}

inline double input_double() {
    double d;
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank==0) {
        std::cin >> d;
    }
    MPI_Bcast(&d, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return d;
}

// ----------------------