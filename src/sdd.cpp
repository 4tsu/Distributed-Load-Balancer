#include <sdd.hpp>



// ============================================

void Sdd::init(const int sdd_type, Variables* vars, Systemparam* sysp
               , const MPIinfo &mi, SubRegion* sr) {
    this->sdd_type = sdd_type;

    if        (sdd_type==0) {
        calc_bounds(sysp, mi);
        return;

    } else if (sdd_type==1) {
        global_sort(vars, sysp, mi, sr);

    } else if (sdd_type==2) {
        calc_bounds(sysp, mi);
        simple(vars, sysp, mi, sr);
        voronoi_init(vars, sysp, mi, sr);
        return;
    }
}



void Sdd::run(Variables* vars, Systemparam* sysp, const MPIinfo &mi, SubRegion* sr) {

    if        (sdd_type==0) {
        simple(vars, sysp, mi, sr);

    } else if (sdd_type==1) {
        global_sort(vars, sysp, mi, sr);

    } else if (sdd_type==2) {
        voronoi(vars, sysp, mi, sr);
    }
}



unsigned long Sdd::ideal(Systemparam* sysp, const MPIinfo &mi) {
    return static_cast<unsigned long>(sysp->N/mi.procs);
}



void Sdd::calc_bounds(Systemparam* sysp, const MPIinfo &mi) {
    int ix = mi.rank%mi.npx;
    int iy = mi.rank/mi.npx;
    const double lpx = sysp->xl/static_cast<double>(mi.npx);
    const double lpy = sysp->yl/static_cast<double>(mi.npy);
    this->bottom = lpy*static_cast<double>(iy);
    this->top    = lpy*static_cast<double>(iy+1);
    this->left   = lpx*static_cast<double>(ix);
    this->right  = lpx*static_cast<double>(ix+1);
}



void Sdd::simple(Variables* vars, Systemparam* sysp, const MPIinfo &mi, SubRegion* sr) {
    std::vector<std::vector<Atom>> migration_atoms(mi.procs);
    std::vector<Atom> new_atoms;
    const double lpx = sysp->xl/static_cast<double>(mi.npx);
    const double lpy = sysp->yl/static_cast<double>(mi.npy);
    for (const Atom atom : vars->atoms) {
        if (atom.x>right || atom.x<left || atom.y>top || atom.y<bottom) {
            int new_rank = static_cast<int>(floor(atom.y/lpy)*mi.npx + floor(atom.x/lpx));
            migration_atoms.at(new_rank).push_back(atom);
        } else {
            new_atoms.push_back(atom);
        }
    }
    vars->atoms = new_atoms;

    // 移動する原子の通信
    // 最初にAllgatherで移動粒子数のテーブルを共有（1次元で送って戻す）
    // 次に必要なペアでだけ通信をする

}
// ============================================