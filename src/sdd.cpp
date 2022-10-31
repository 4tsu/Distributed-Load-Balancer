#include <sdd.hpp>



// ============================================

void Sdd::init(const int sdd_type, Variables* vars, Systemparam* sysp
               , const MPIinfo &mi, SubRegion* sr) {
    this->sdd_type = sdd_type;

    if        (sdd_type==0) {
        return;

    } else if (sdd_type==1) {
        global_sort(vars, sysp, mi, sr);

    } else if (sdd_type==2) {
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

// ============================================