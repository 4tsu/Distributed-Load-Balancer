#pragma once

// ==========================================
struct MPIinfo {
    int rank;
    int procs;
    int npx, npy;
};

// ------------------------------------------

void setup_info(MPIinfo &mi);

// ==========================================