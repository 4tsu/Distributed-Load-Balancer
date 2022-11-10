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
    void arrange_pairs(unsigned long);
    void clear_all(void);
    double lmx, lmy;
    int nmx, nmy, num_mesh;
    std::vector<unsigned long> counter, head_index, sorted_index;
    std::vector<double> limits;
 
    void set_mesh_ext(Systemparam*);
    void set_index_ext(const std::vector<Atom> &, Systemparam*,
         std::vector<unsigned long> &, std::vector<unsigned long> &, std::vector<unsigned long> &, std::vector<unsigned long> &, std::vector<unsigned long> &);
    void search_ext(int, int, const std::vector<Atom> &, const std::vector<Atom> &, Systemparam*,
         const std::vector<unsigned long> &, const std::vector<unsigned long> &, const std::vector<unsigned long> &);
    void mesh_search_ext(const std::vector<Atom> &, std::vector<Atom> &, Systemparam*);
    void arrange_pairs_ext(const std::vector<unsigned long> &new_index, unsigned long new_pn);
    void pick_atoms(std::vector<Atom> &, std::vector<unsigned long> &, std::vector<unsigned long> &);
    void clear_ext(void);
    int nmx_ext, nmy_ext, num_mesh_ext;
    std::vector<Pair> one_other_list;
    std::vector<bool> across_border;
    std::vector<bool> survivor_list;
   
public:
    std::vector<Pair> list;
    std::vector<std::vector<Pair>> other_list;
    void make_pair(Variables* vars, Systemparam*);
};

// ============================================