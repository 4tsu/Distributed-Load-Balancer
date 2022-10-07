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
    Atom *other_atoms = vars->other_atoms.data();
    int other_pn = vars->number_of_other_atoms();
    for (int i=0; i<pn; i++) {
        double ix = atoms[i].x;
        double iy = atoms[i].y;
        for (int j=0; j<other_pn; j++) {
            double dx = other_atoms[j].x - ix;
            double dy = other_atoms[j].y - iy;
            periodic_distance(dx, dy, sysp);
            double r = sqrt(dx*dx + dy*dy);
            if (r > sysp->co_margin) {
                continue;
            }
            Pair pair;
            set_pair(pair, i, j, atoms[i].id, other_atoms[j].id);
            this->other_list.push_back(pair);
        }
    }
}

// ============================================