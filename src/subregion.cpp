#include "subregion.hpp"

// ============================================

void set_dp(DomainPair &dp, int ip, int jp) {
    dp.i = ip;
    dp.j = jp;
}



std::vector<double> calc_limit(Variables* vars) {
    if (vars->number_of_atoms() > 0) {
        double x_min = vars->atoms.at(0).x;
        double x_max = vars->atoms.at(0).x;
        double y_min = vars->atoms.at(0).y;
        double y_max = vars->atoms.at(0).y;
        for (auto atom : vars->atoms) {
            if      (atom.x < x_min)
                x_min = atom.x;
            else if (atom.x > x_max)
                x_max = atom.x;
            if      (atom.y < y_min)
                y_min = atom.y;
            else if (atom.y > y_max)
                y_max = atom.y;
        }
        std::vector<double> v{x_min, x_max, y_min, y_max};
        return v;
    } else {
        std::vector<double> v{0, 0, 0, 0};
        return v;
    }
}

// --------------------------------------------

// 相互作用する可能性のある領域ぺア検出
void SubRegion::make_dplist(MPIinfo mi, Variables* vars, Systemparam* sysp) {
    /// まずは中心と半径を計算して共有する
    this->calc_center(vars, sysp);
    this->calc_radius(vars, sysp);
    this->communicate_centradi(mi);
    this->dplist.clear();
    this->dplist_reverse.clear();
    DomainPair dp;
    if (mi.procs <= 1)
        return;
    
    /// O(np^2)で全探索
    /// 作用反作用を考慮して半分に落とすのを、計算が効率よくなるようにやる
    /// なるべく全ての領域で同じくらいの長さの領域ペアリストを持ちたい
    //// dplistとdplist_reverseを分けて、二重ループを展開できるはず。
    for (int i=0; i<(mi.procs/2); i++) {
        for (int j=i+1; j<(mi.procs/2+i+1); j++) {
            if (i!=mi.rank && j!=mi.rank)
                continue;
            if (this->judge(i,j,sysp))
                continue;
            if (i==mi.rank) {
                set_dp(dp, i, j);
                this->dplist.push_back(dp);
            } else if (j==mi.rank){
                set_dp(dp, j, i);
                this->dplist_reverse.push_back(dp);
            }
        }
    }
    for (int i=mi.procs/2; i<mi.procs-1; i++) {
        for (int j=0; j<(i-(mi.procs/2)); j++) {
            if (i!=mi.rank && j!=mi.rank)
                continue;
            if (this->judge(i,j,sysp))
                continue;
            if (i==mi.rank) {
                set_dp(dp, i, j);
                this->dplist.push_back(dp);
            } else if (j==mi.rank){
                set_dp(dp, j, i);
                this->dplist_reverse.push_back(dp);
            }
        }
        for (int j=i+1; j<mi.procs; j++) {
            if (i!=mi.rank && j!=mi.rank)
                continue;
            if (this->judge(i,j,sysp))
                continue;
            if (i==mi.rank) {
                set_dp(dp, i, j);
                this->dplist.push_back(dp);
            } else if (j==mi.rank){
                set_dp(dp, j, i);
                this->dplist_reverse.push_back(dp);
            }
        }
    }
    for (int i=mi.procs-1; i<mi.procs; i++) {
        for (int j=0; j<(i-mi.procs/2); j++) {
            if (i!=mi.rank && j!=mi.rank)
                continue;
            if (this->judge(i,j,sysp))
                continue;
            if (i==mi.rank) {
                set_dp(dp, i, j);
                this->dplist.push_back(dp);
            } else if (j==mi.rank){
                set_dp(dp, j, i);
                this->dplist_reverse.push_back(dp);
            }
        }
    }
}



bool SubRegion::judge(int i, int j, Systemparam* sysp) {
    if (i==j)
        return true;
    if (is_empties.at(i) || is_empties.at(j))
        return true;
    double dx = centers.at(i).at(0) - centers.at(j).at(0);
    double dy = centers.at(i).at(1) - centers.at(j).at(1);
    periodic_distance(dx, dy, sysp);
    double gap = sqrt(dx*dx + dy*dy) - radii.at(i) - radii.at(j);
    if (gap > sysp->co_margin)
        return true;
    return false;
}



void SubRegion::communicate_centradi(const MPIinfo &mi) {
    this->centers.clear();
    this->radii.clear();
    this->is_empties.clear();
    std::vector<double> recvbuf_x(mi.procs);
    std::vector<double> recvbuf_y(mi.procs);
    std::vector<double> recvbuf_r(mi.procs);
    std::vector<int> recvbuf_e(mi.procs);
    MPI_Allgather(&center[0], 1, MPI_DOUBLE, recvbuf_x.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&center[1], 1, MPI_DOUBLE, recvbuf_y.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&radius, 1, MPI_DOUBLE, recvbuf_r.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&is_empty, 1, MPI_INT, recvbuf_e.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for (int i=0; i<mi.procs; i++) {
        std::vector<double> one_center = {recvbuf_x.at(i), recvbuf_y.at(i)};
        this->centers.push_back(one_center);
        this->radii.push_back(recvbuf_r.at(i));
        this->is_empties.push_back(recvbuf_e.at(i));
    }
}



void SubRegion::calc_center(Variables* vars, Systemparam* sysp) {
    const unsigned long pn = vars->atoms.size();

    if (pn == 0) {
        this->center[0] = 0.0;
        this->center[1] = 0.0;
        this->is_empty = true;
        return;
    }
    
    is_empty = false;
    const double origin_ax = vars->atoms.at(0).x;
    const double origin_ay = vars->atoms.at(0).y;
    double sx = 0;
    double sy = 0;
    for (auto atom: vars->atoms) {
        double dx = atom.x - origin_ax;
        double dy = atom.y - origin_ay;
        periodic_distance(dx, dy, sysp);
        sx += dx;
        sy += dy;
    }
    sx /= static_cast<double>(pn);
    sy /= static_cast<double>(pn);
    sx += origin_ax;
    sy += origin_ay;
    periodic_coordinate(sx, sy, sysp);
    this->center[0] = sx;
    this->center[1] = sy;
}



void SubRegion::calc_radius(Variables* vars, Systemparam* sysp) {
    double r_max = 0;
    for (auto atom: vars->atoms) {
        double dx = atom.x - this->center[0];
        double dy = atom.y - this->center[1];
        periodic_distance(dx, dy, sysp);
        double r = sqrt(dx*dx + dy*dy);
        if (r>r_max) {
            r_max = r;
        }
    }
    this->radius = r_max;
}

// ============================================