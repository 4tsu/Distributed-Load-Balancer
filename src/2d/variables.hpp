#pragma once
#include "mpiinfo.hpp"
#include "systemparam.hpp"

// ========================================================================

struct Atom {
    unsigned long id;
    double  x,  y;
    double vx, vy;
};




struct Force {
    unsigned long id;
    double vx, vy;
};



class Variables {
public:
    std::vector<Atom> atoms;
    std::vector<std::vector<Atom>> other_atoms;

    // 自分が受け取りたい他領域粒子の情報（MPI_Irecvで使用）
    std::vector<unsigned long> recv_size;   // バイト単位、Atomベクターの大きさ
    std::vector<std::vector<unsigned long>> recv_list;   // IDのみ格納
    // 他領域が受け取りたい自領域粒子の情報（MPI_Isendで使用）
    std::vector<unsigned long> send_size;   // バイト単位、ペアリスト構築時に他領域から受け取っておく
    std::vector<std::vector<unsigned long>> send_list;   // 受け取っておく
    std::vector<std::vector<Atom*>> send_atoms;   // 上の受け取りを受けて作成する
    std::vector<std::vector<Force>> sending_force;   // 力積の書き戻しに使用
    
    double time;
    Variables(void) {time = 0.0;}
    void add_atoms(unsigned long id, double x, double y);
    void export_cdview(void);
    unsigned long number_of_atoms(void) {return atoms.size();}
    void set_initial_velocity(const double, MPIinfo);
    double margin_life;
    void set_margin_life(double);
    double max_velocity(void);
    void pack_send_atoms(void);

    double xp_max;
    double xp_min;
    double yp_max;
    double yp_min;

};

// ========================================================================