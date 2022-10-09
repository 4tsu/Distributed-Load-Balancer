#include "pairlist.hpp"

// ============================================

void set_pair(Pair &pair, int i, int j, int idi, int idj) {
    pair.i = i;
    pair.j = j;
    pair.idi = idi;
    pair.idj = idj;
}

// ---------------------------------------

// ペアリスト作成・更新
// ローカルにペアリスト作成
void PairList::make_pair(Variables *vars, Systemparam *sysp) {
    // 自領域内ペア
    list.clear();
    Atom *atoms = vars->atoms.data();
    int pn = vars->number_of_atoms();
    for (int i=0; i<pn; i++) {
        double ix = atoms[i].x;
        double iy = atoms[i].y;
        for (int j=i+1; j<pn; j++) {
            double dx = atoms[j].x - ix;
            double dy = atoms[j].y - iy;
            periodic_distance(dx, dy, sysp);
            double r = sqrt(dx*dx + dy*dy);
            if (r > sysp->co_margin) {
                continue;
            }
            Pair pair;
            set_pair(pair, i, j, atoms[i].id, atoms[j].id);
            this->list.push_back(pair);
        }
    }
    // 他領域粒子とのペア
    other_list.clear();
    std::vector<Pair> one_other_list;
    // ペアリストに載らない粒子情報はother_atomsから削除する
    std::vector<std::vector<Atom>> new_other_atom;
    for (auto &one_other_atoms : vars->other_atoms) {
        Atom *other_atoms = one_other_atoms.data();
        const int other_pn = one_other_atoms.size();
        std::vector<Atom> one_new_other_atom;
        for (int i=0; i<other_pn; i++) {
            double ix = other_atoms[i].x;
            double iy = other_atoms[i].y;
            bool survive = false;
            for (int j=0; j<pn; j++) {
                double dx = atoms[j].x - ix;
                double dy = atoms[j].y - iy;
                periodic_distance(dx, dy, sysp);
                double r = sqrt(dx*dx + dy*dy);
                if (r > sysp->co_margin) {
                    continue;
                }
                Pair pair;
                set_pair(pair, j, i, atoms[j].id, other_atoms[i].id);
                one_other_list.push_back(pair);
                survive = true;
            }
            if (survive)
                one_new_other_atom.push_back(other_atoms[i]);
        }
        this->other_list.push_back(one_other_list);
        new_other_atom.push_back(one_new_other_atom);
    }
    vars->other_atoms = new_other_atom;
    assert(this->other_list.size() == vars->other_atoms.size());
}

// ============================================