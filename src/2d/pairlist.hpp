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
    void set_mesh(Variables*);
    void set_index(Variables*);
    void search(long, Variables*);
    void search_neighbor(long, long, long, Variables*);
    void mesh_search(Variables*);
    void search_all(Variables*);
    void arrange_pairs(unsigned long);
    void clear_all(void);
    double lmx, lmy;
    long nmx, nmy, num_mesh;
    std::vector<unsigned long> counter, head_index;
    std::vector<Atom> sorted_atoms;
    std::vector<double> limits;
 
    void set_mesh_ext(void);
    void set_index_ext(const std::vector<Atom> &);
    void search_ext(long, long, long, const std::vector<Atom> &, const std::vector<Atom> &);
    void mesh_search_ext(const std::vector<Atom> &, std::vector<Atom> &);
    void arrange_pairs_ext(const std::vector<unsigned long> &new_index, unsigned long new_pn);
    void pick_atoms(std::vector<Atom> &, std::vector<unsigned long> &, std::vector<unsigned long> &);
    void search_all_ext(const std::vector<Atom> &, const std::vector<Atom> &);
    void clear_ext(void);
    long nmx_ext, nmy_ext, num_mesh_ext;
    std::vector<Pair> one_other_list;
    std::vector<unsigned long> counter_ext, head_index_ext, sorted_index_ext;
    std::vector<bool> survivor_list;
   
public:
    std::vector<Pair> list;
    std::vector<std::vector<Pair>> other_list;
    void make_pair(Variables* vars);
    void make_pair_ext(Variables* vars);
};

// ============================================