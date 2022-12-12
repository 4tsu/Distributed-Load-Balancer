#pragma once

// ==========================================
struct MPIinfo {
    int rank;
    int procs;
    int npx, npy, npz;
};

// ------------------------------------------

void setup_info(MPIinfo &mi);

// ==========================================