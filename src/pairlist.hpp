#pragma once
#include "variables.hpp"
#include "systemparam.hpp"
#include "subregion.hpp"

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
private:
    void set_mesh(Variables*, Systemparam*);
    void set_index(Variables*, Systemparam*);
    void search(int, Variables*, Systemparam*);
    void search_neighbor(int, int, Variables*, Systemparam*);
    void mesh_search(Variables*, Systemparam*);
    void clear_all(void);
    double lmx, lmy;
    int nmx, nmy, num_mesh;
    std::vector<unsigned long> counter, head_index, sorted_index;
    std::vector<double> limits;
public:
    std::vector<Pair> list;
    std::vector<std::vector<Pair>> other_list;
    void make_pair(Variables* vars, Systemparam*);
};

// ============================================