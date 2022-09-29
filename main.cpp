#include <time.h>
#include <stdio.h>
#include "md.hpp"

// ========================================

int main(int argc, char **argv) {
    clock_t start = clock();

    MD *md = new MD();
    md->set_params(1000, 10, 0.0020);
    md->run();

    delete md;
    clock_t end = clock();
    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC*1000.0;
    printf("time %lf[ms]\n", time);
}
// ========================================