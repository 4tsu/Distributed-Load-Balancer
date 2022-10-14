#pragma once

#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <cassert>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <filesystem>
#include <string>

#include <mpi.h>

// ----------------------

//デバッグ用
inline void ckpt() {
    static int c = 1;
    fprintf(stderr, "c.%d\n", c);
    c++;
}
inline void ckpt2() {
    static int d = 1;
    fprintf(stderr, "d.%d\n", d);
    d++;
}
inline void ckpt3() {
    static int e = 1;
    fprintf(stderr, "e.%d\n", e);
    e++;
}

// ----------------------