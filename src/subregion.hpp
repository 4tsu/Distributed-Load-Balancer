#pragma once
#include "mpiinfo.hpp"
#include "variables.hpp"
#include "systemparam.hpp"

// ============================================

// 粒子情報をやり取りするべきプロセス（SubRegion）のペア
struct DomainPair{
    int i,j;   // ランク番号のペア
};

// --------------------------------------------

// 空間分割し、1プロセスに1つ割り当てられた領域
class SubRegion {
private:
    bool judge(int, int, Systemparam*);
    void communicate_centradi(const MPIinfo &mi);
    std::vector<std::vector<double>> centers;
    std::vector<double> radii;
    std::vector<bool> is_empties;
    bool is_empty;

public:
    std::vector<DomainPair> dplist;
    std::vector<DomainPair> dplist_reverse;
    void make_dplist(MPIinfo, Variables*, Systemparam*);

    double center[2];
    double radius;
    void calc_center(Variables*, Systemparam*);
    void calc_radius(Variables*, Systemparam*);
    void set_bias(double);
};

// --------------------------------------------

std::vector<double> calc_limit(Variables*);

// ============================================