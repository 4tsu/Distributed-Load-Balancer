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
 
    void set_mesh_ext(const std::vector<Atom> &, Systemparam*);
    void set_index_ext(const std::vector<Atom> &, Systemparam*,
         std::vector<unsigned long> &, std::vector<unsigned long> &, std::vector<unsigned long> &, std::vector<unsigned long> &, std::vector<unsigned long> &);
    void search_ext(int, int, const std::vector<Atom> &, const std::vector<Atom> &, Systemparam*,
         const std::vector<unsigned long> &, const std::vector<unsigned long> &, const std::vector<unsigned long> &, const std::vector<unsigned long> &, const std::vector<unsigned long> &);
    void mesh_search_ext(const std::vector<Atom> &, Systemparam*);
    void clear_ext(void);
    int nmx_ext, nmy_ext, num_mesh_ext;
   
public:
    std::vector<Pair> list;
    std::vector<std::vector<Pair>> other_list;
    void make_pair(Variables* vars, Systemparam*);
};

// ============================================