#include "pairlist.hpp"

// ============================================

void set_pair(Pair &pair, unsigned long i, unsigned long j, unsigned long idi, unsigned long idj) {
    pair.i = i;
    pair.j = j;
    pair.idi = idi;
    pair.idj = idj;
}

// ---------------------------------------

// ペアリスト作成・更新
// ローカルにペアリスト作成
void PairList::make_pair(Variables* vars, Systemparam* sysp) {
    // 自領域内ペア
    list.clear();
    this->mesh_search(vars, sysp);
    // 他領域粒子とのペア
    // ペアリストに載らない粒子情報はother_atomsから削除する
    // 通信がいらなくなった他領域も、DomainPairListから削除する
    // 他領域から情報を受け取るべき粒子のリストcomm_recv_listも作成
    this->other_list.clear();
    vars->recv_list.clear();
    vars->recv_size.clear();
    unsigned long pn = vars->number_of_atoms();
    Atom* atoms = vars->atoms.data();
    std::vector<std::vector<Atom>> new_other_atom;
    for (auto &one_other_atoms : vars->other_atoms) {
        Atom *other_atoms = one_other_atoms.data();
        const unsigned long other_pn = one_other_atoms.size();
        std::vector<Pair> one_other_list;
        std::vector<Atom> one_new_other_atom;
        std::vector<unsigned long> one_recv_list;
        unsigned long pair_j_index = 0;   // other_atomを間引いてしまうので、ペアのインデックスは新しく用意
        for (unsigned long i=0; i<other_pn; i++) {
            double ix = other_atoms[i].x;
            double iy = other_atoms[i].y;
            bool survive = false;
            for (unsigned long j=0; j<pn; j++) {
                double dx = atoms[j].x - ix;
                double dy = atoms[j].y - iy;
                periodic_distance(dx, dy, sysp);
                double r = sqrt(dx*dx + dy*dy);
                if (r > sysp->co_margin) {
                    continue;
                }
                Pair pair;
                set_pair(pair, j, pair_j_index, atoms[j].id, other_atoms[i].id);
                one_other_list.push_back(pair);
                survive = true;
            }
            if (survive){
                one_new_other_atom.push_back(other_atoms[i]);
                one_recv_list.push_back(other_atoms[i].id);
                pair_j_index++;
            }
        }
        if (one_other_list.size() != 0)
            this->other_list.push_back(one_other_list);
        new_other_atom.push_back(one_new_other_atom);
        vars->recv_list.push_back(one_recv_list);
        vars->recv_size.push_back(one_recv_list.size()*sizeof(Atom));
    }
    vars->other_atoms = new_other_atom;

for (auto l : list) {
    printf("%lu-%lu\n", l.idi, l.idj);
}
}



void PairList::set_mesh(Variables* vars, Systemparam* sysp) {
    this->limits = calc_limit(vars);
    this->nmx = static_cast<int>((limits.at(1) - limits.at(0))/sysp->co_margin) - 1;
    this->lmx = static_cast<double>(limits.at(1) - limits.at(0))/nmx;
    this->nmy = static_cast<int>((limits.at(3) - limits.at(2))/sysp->co_margin) - 1;
    this->lmy = static_cast<double>(limits.at(3) - limits.at(2))/nmy;
    assert(nmx>2 && nmy>2);
    assert(lmx>sysp->co_margin && lmy>sysp->co_margin);
    this->num_mesh = nmx*nmy;
    assert(num_mesh < std::numeric_limits<int>::max());
    this->counter.resize(num_mesh);
    this->head_index.resize(num_mesh);
    this->sorted_index.resize(vars->number_of_atoms());
}



void PairList::set_index(Variables* vars, Systemparam* sysp) {
    unsigned long pn = vars->number_of_atoms();
    std::vector<unsigned long> position_buffer(pn);
    for (unsigned long i=0; i<pn; i++) {
        int ix = static_cast<int>((vars->atoms.at(i).x-limits.at(0))/lmx);
        int iy = static_cast<int>((vars->atoms.at(i).y-limits.at(2))/lmy);
        if (vars->atoms.at(i).x==limits.at(1))
            ix = nmx-1;
        if (vars->atoms.at(i).y==limits.at(3))
            iy = nmy-1;
        assert(0<=ix && ix<nmx);
        assert(0<=iy && iy<nmy);
        int im = ix + iy*nmx;
        assert(0<=im && im<num_mesh);
        counter.at(im)++;
        position_buffer.at(i) = im;
    }
    unsigned long total = std::accumulate(counter.begin(), counter.end(), 0);
    assert(total==pn);

    head_index.at(0) = 0;
    int sum = 0;
    for (int i=1; i<num_mesh; i++) {
        sum += counter.at(i-1);
        head_index.at(i) = sum;
    }

    std::vector<unsigned long> indexes(num_mesh);
    std::fill(indexes.begin(), indexes.end(), 0);
    for (unsigned long i=0; i<pn; i++) {
        int im = position_buffer.at(i);
        unsigned long j = head_index.at(im) + indexes.at(im);
        assert(sorted_index.at(j) == 0);
        sorted_index.at(j) = i;
        indexes.at(im)++;
    }
}



void PairList::search(int im, Variables* vars, Systemparam* sysp) {
    int ih = head_index.at(im);
    unsigned long pnm = counter.at(im);
    Atom *atoms = vars->atoms.data();
    for (unsigned long m=ih; m<ih+pnm-1; m++) {
        for (unsigned long n=m+1; n<ih+pnm; n++) {
            int i = sorted_index.at(m);
            int j = sorted_index.at(n);
            double dx = atoms[j].x - atoms[i].x;
            double dy = atoms[j].y - atoms[i].y;
            periodic_distance(dx, dy, sysp);
            double r = sqrt(dx*dx + dy*dy);
            if (r > sysp->co_margin)
                continue;
            Pair p;
            set_pair(p, i, j, atoms[i].id, atoms[j].id);
            this->list.push_back(p);
        }
    }
}



void PairList::search_neighbor(int im, int jm, Variables* vars, Systemparam* sysp) {
    if (jm<0) {
        jm += num_mesh;
    } else if (jm>=num_mesh) {
        jm -= num_mesh;
    }
    Atom *atoms = vars->atoms.data();
    for (unsigned long m=head_index[im]; m<head_index[im]+counter[im]; m++) {
        for (unsigned long n=head_index[jm]; n<head_index[jm]+counter[jm]; n++) {
            int i = sorted_index.at(m);
            int j = sorted_index.at(n);
            double dx = atoms[j].x - atoms[i].x;
            double dy = atoms[j].y - atoms[i].y;
            periodic_distance(dx, dy, sysp);
            double r = sqrt(dx*dx + dy*dy);
            if (r>sysp->co_margin)
                continue;
            Pair p;
            set_pair(p, i, j, atoms[i].id, atoms[j].id);
            this->list.push_back(p);
        }
    }
}



void PairList::mesh_search(Variables* vars, Systemparam* sysp) {
    clear_all();
    set_mesh(vars, sysp);
    set_index(vars, sysp);
    for (int i=0; i<num_mesh; i++) {
        search(i, vars, sysp);
        search_neighbor(i, i+1, vars, sysp);
        search_neighbor(i, i+nmx, vars, sysp);
        search_neighbor(i, i+nmx+1, vars, sysp);
        search_neighbor(i, i-nmx+1, vars, sysp);
    }
}



void PairList::clear_all(void) {
    lmx = 0;
    lmy = 0;
    nmx = 0;
    nmy = 0;
    num_mesh = 0;
    counter.clear();
    head_index.clear();
    sorted_index.clear();
    limits.clear();
}
// ============================================