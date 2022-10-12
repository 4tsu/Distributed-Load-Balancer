#pragma once
#include "variables.hpp"
#include "systemparam.hpp"
#include "domainpair.hpp"

// ============================================

// ペアリストに使う
struct Pair {
    int i, j;
    int idi, idj;
};

// ---------------------------------------

void set_pair(Pair &pair, int i, int j, int idi, int idj);

// ---------------------------------------

// ペアリストクラス
class PairList {
public:
    std::vector<Pair> list;
    std::vector<std::vector<Pair>> other_list;
    void make_pair(Variables* &vars, Systemparam*, DomainPairList*);
};

// ============================================