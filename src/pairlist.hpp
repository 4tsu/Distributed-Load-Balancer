#pragma once
#include "variables.hpp"
#include "systemparam.hpp"

// ============================================

// ペアリストに使う
struct Pair {
    unsigned long i, j;
    unsigned long idi, idj;
};

// ---------------------------------------

void set_pair(Pair &pair, unsigned long i, unsigned long j, unsigned long idi, unsigned long idj);

// ---------------------------------------

// ペアリストクラス
class PairList {
public:
    std::vector<Pair> list;
    std::vector<std::vector<Pair>> other_list;
    void make_pair(Variables* vars, Systemparam*);
};

// ============================================