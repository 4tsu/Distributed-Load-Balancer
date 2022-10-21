#include "subdomain.hpp"

// ============================================

void set_dp(DomainPair &dp, int ip, int jp) {
    dp.i = ip;
    dp.j = jp;
}

// --------------------------------------------

// 4プロセス並列を仮定、手動でリストを構築している
void SubDomain::make_list(MPIinfo mi) {
    dplist.clear();
    dplist_reverse.clear();
    DomainPair dp;
    if (mi.procs <= 1)
        return;
    if (mi.rank == 0) {
        set_dp(dp, 0, 1);
        this->dplist.push_back(dp);
        set_dp(dp, 0, 2);
        this->dplist.push_back(dp);
        set_dp(dp, 0, 3);
        this->dplist_reverse.push_back(dp);
    } else if (mi.rank == 1) {
        set_dp(dp, 1, 2);
        this->dplist.push_back(dp);
        set_dp(dp, 1, 3);
        this->dplist.push_back(dp);
        set_dp(dp, 1, 0);
        this->dplist_reverse.push_back(dp);
    } else if (mi.rank == 2) {
        set_dp(dp, 2, 3);
        this->dplist.push_back(dp);
        set_dp(dp, 2, 0);
        this->dplist_reverse.push_back(dp);
        set_dp(dp, 2, 1);
        this->dplist_reverse.push_back(dp);
    } else {
        set_dp(dp, 3, 0);
        this->dplist.push_back(dp);
        set_dp(dp, 3, 1);
        this->dplist_reverse.push_back(dp);
        set_dp(dp, 3, 2);
        this->dplist_reverse.push_back(dp);
    }
}

// ============================================