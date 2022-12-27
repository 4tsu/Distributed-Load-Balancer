#include "pairlist.hpp"

namespace sysp = systemparam;

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
void PairList::make_pair(Variables* vars) {
    // 自領域内ペア
    list.clear();
    this->mesh_search(vars);
    this->arrange_pairs(vars->number_of_atoms());
}


// 他領域粒子とのペア
// ペアリストに載らない粒子情報はother_atomsから削除する
// 通信がいらなくなった他領域も、DomainPairListから削除する
// 他領域から情報を受け取るべき粒子のリストcomm_recv_listも作成
void PairList::make_pair_ext(Variables* vars) {
    this->other_list.clear();
    vars->recv_list.clear();
    vars->recv_size.clear();
    for (auto &one_other_atoms : vars->other_atoms) {
        mesh_search_ext(vars->atoms, one_other_atoms);
        std::vector<unsigned long> one_recv_list;
        std::vector<unsigned long> new_index;
        pick_atoms(one_other_atoms, one_recv_list, new_index);
        vars->recv_list.push_back(one_recv_list);
        vars->recv_size.push_back(one_recv_list.size()*sizeof(Atom));
        arrange_pairs_ext(new_index, one_other_atoms.size());
        if (one_other_list.size() != 0)
            this->other_list.push_back(one_other_list);
    }
}



void PairList::set_mesh(Variables* vars) {
    this->limits = calc_limit(vars);
    this->nmx = static_cast<long>((limits.at(1) - limits.at(0))/sysp::co_margin) - 1;
    this->nmy = static_cast<long>((limits.at(3) - limits.at(2))/sysp::co_margin) - 1;
    if (nmx<3 || nmy<3)
        return;
    this->lmx = (limits.at(1) - limits.at(0))/static_cast<double>(nmx);
    this->lmy = (limits.at(3) - limits.at(2))/static_cast<double>(nmy);
    assert(lmx>sysp::co_margin && lmy>sysp::co_margin);
    this->num_mesh = nmx*nmy;
    assert((nmx+2)*(nmy+2)< std::numeric_limits<long>::max());
    this->counter.resize(num_mesh);
    this->head_index.resize(num_mesh);
}



void PairList::set_index(Variables* vars) {
    unsigned long pn = vars->number_of_atoms();
    std::vector<unsigned long> position_buffer(pn);
    for (unsigned long i=0; i<pn; i++) {
        long ix = static_cast<long>((vars->atoms.at(i).x-limits.at(0))/lmx);
        long iy = static_cast<long>((vars->atoms.at(i).y-limits.at(2))/lmy);
        if (vars->atoms.at(i).x==limits.at(1))
            ix = nmx-1;
        if (vars->atoms.at(i).y==limits.at(3))
            iy = nmy-1;
        assert(0<=ix && ix<nmx);
        assert(0<=iy && iy<nmy);
        long im = ix + iy*nmx;
        assert(0<=im && im<num_mesh);
        counter.at(im)++;
        position_buffer.at(i) = im;
    }
    unsigned long total = std::accumulate(counter.begin(), counter.end(), 0);
    assert(total==pn);

    head_index.at(0) = 0;
    unsigned long sum = 0;
    for (long i=1; i<num_mesh; i++) {
        sum += counter.at(i-1);
        head_index.at(i) = sum;
    }

    std::vector<unsigned long> indexes(num_mesh);
    std::fill(indexes.begin(), indexes.end(), 0);
    sorted_atoms.resize(pn);
    for (unsigned long i=0; i<pn; i++) {
        long im = position_buffer.at(i);
        unsigned long j = head_index.at(im) + indexes.at(im);
        sorted_atoms.at(j) = vars->atoms.at(i);
        indexes.at(im)++;
    }
    std::copy(sorted_atoms.begin(), sorted_atoms.end(), vars->atoms.begin());
}



void PairList::search(long im, Variables* vars) {
    unsigned long ih = head_index.at(im);
    unsigned long pnm = counter.at(im);
    Atom *atoms = vars->atoms.data();
    for (unsigned long i=ih; i<ih+pnm; i++) {
        for (unsigned long j=i+1; j<ih+pnm; j++) {
            if (i==j)
                continue;
            double dx = atoms[j].x - atoms[i].x;
            double dy = atoms[j].y - atoms[i].y;
            periodic_distance(dx, dy);
            double r = std::sqrt(dx*dx + dy*dy);
            if (r > sysp::co_margin)
                continue;
            Pair p;
            set_pair(p, i, j, atoms[i].id, atoms[j].id);
            this->list.push_back(p);

        }
    }
}



void PairList::search_neighbor(long im, long jmx, long jmy, Variables* vars) {
    if (jmx<0) {
        jmx += nmx;
    } else if (jmx>=nmx) {
        jmx -= nmx;
    }
    if (jmy<0) {
        jmy += nmy;
    } else if (jmy>=nmy) {
        jmy -= nmy;
    }
    long jm = jmx + jmy*nmx;
    Atom *atoms = vars->atoms.data();
    for (unsigned long i=head_index[im]; i<head_index[im]+counter[im]; i++) {
        for (unsigned long j=head_index[jm]; j<head_index[jm]+counter[jm]; j++) {
            double dx = atoms[j].x - atoms[i].x;
            double dy = atoms[j].y - atoms[i].y;
            periodic_distance(dx, dy);
            double r = std::sqrt(dx*dx + dy*dy);
            if (r>sysp::co_margin)
                continue;
            Pair p;
            set_pair(p, i, j, atoms[i].id, atoms[j].id);
            this->list.push_back(p);
       }
    }
}



void PairList::mesh_search(Variables* vars) {
    clear_all();
    set_mesh(vars);
    if (nmx>2 && nmy>2) {
        set_index(vars);
        for (long i=0; i<num_mesh; i++) {
            search(i, vars);
            long ix = i%nmx;
            long iy = i/nmx;
            search_neighbor(i, ix+1, iy  , vars);
            search_neighbor(i, ix, iy+1  , vars);
            search_neighbor(i, ix+1, iy+1, vars);
            search_neighbor(i, ix+1, iy-1, vars);
        }
    } else {
        search_all(vars);
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
    sorted_atoms.clear();
    limits.clear();
}



void PairList::arrange_pairs(unsigned long pn) {
    unsigned long ln = list.size();
    std::vector<Pair> new_list(ln);
    std::vector<unsigned long> pair_counter(pn);
    std::fill(pair_counter.begin(), pair_counter.end(), 0);
    for (auto p : list) {
        pair_counter.at(p.i)++;
    }

    std::vector<unsigned long> head_index_pair(pn);
    unsigned long sum = 0;
    for (unsigned long i=1; i<pn; i++) {
        sum += pair_counter.at(i-1);
        head_index_pair.at(i) = sum;
    }

    std::vector<unsigned long> indexes(pn);
    std::fill(indexes.begin(), indexes.end(), 0);
    for (std::size_t i=0; i<ln; i++) {
        Pair p = list.at(i);
        unsigned long j = head_index_pair.at(p.i) + indexes.at(p.i);
        new_list.at(j) = p;
        indexes.at(p.i)++;
    }

    list.clear();
    list.resize(ln);
    std::copy(new_list.begin(), new_list.end(), list.begin());
}



void PairList::search_all(Variables* vars) {
    list.clear();
    Atom *atoms = vars->atoms.data();
    const unsigned long pn = vars->number_of_atoms();
    for (unsigned long i=0; i<pn; i++) {
        double ix = atoms[i].x;
        double iy = atoms[i].y;
        for (unsigned long j=i+1; j<pn; j++) {
            double dx = atoms[j].x - ix;
            double dy = atoms[j].y - iy;
            periodic_distance(dx, dy);
            double r = std::sqrt(dx*dx + dy*dy);
            if (r > sysp::co_margin) {
                continue;
            }
            Pair pair;
            set_pair(pair, i, j, atoms[i].id, atoms[j].id);
            this->list.push_back(pair);
        }
    }
}


void PairList::set_mesh_ext(void) {
    nmx_ext = nmx+2;
    nmy_ext = nmy+2;
    num_mesh_ext = nmx_ext*nmy_ext;

    counter_ext.resize(num_mesh_ext);
    head_index_ext.resize(num_mesh_ext);
}



void PairList::set_index_ext(const std::vector<Atom> &atoms) {
    std::vector<unsigned long> inside_mesh_index;
    std::vector<unsigned long> position_buffer;
    unsigned long pn = atoms.size();
    for (unsigned long i=0; i<pn; i++) {
        double dx = atoms.at(i).x - limits.at(0);
        double dy = atoms.at(i).y - limits.at(2);
        long ix = std::floor(dx/lmx)+1;
        long iy = std::floor(dy/lmy)+1;
        long im = -1;
        if (0<=ix && ix<nmx_ext && 0<=iy && iy<nmy_ext) {
            im = ix+iy*nmx_ext;
        } else {
            bool flag = false;
            for (int xc=0; xc<2; xc++) {
                for (int yc=2; yc<4; yc++) {
                    dx = atoms.at(i).x - limits.at(xc);
                    dy = atoms.at(i).y - limits.at(yc);
                    periodic_distance(dx, dy);
                    ix = std::floor(dx/lmx)+nmx*xc+1;
                    iy = std::floor(dy/lmy)+nmy*(yc-2)+1;
                    if (0<=ix && ix<nmx_ext && 0<=iy && iy<nmy_ext) {
                        im = ix+iy*nmx_ext;
                        flag = true;
                        break;
                    }
                }
                if (flag) break;
            }
        }
        
        if (im>=0) {
            inside_mesh_index.push_back(i);
            counter_ext.at(im)++;
            position_buffer.push_back(im);
        }
    }
    head_index_ext.at(0) = 0;
    unsigned long sum = 0;
    for (long i=1; i<num_mesh_ext; i++) {
        sum += counter_ext.at(i-1);
        head_index_ext.at(i) = sum;
    }

    sorted_index_ext.resize(inside_mesh_index.size());
    std::vector<unsigned long> indexes(num_mesh_ext);
    std::fill(indexes.begin(), indexes.end(), 0);
    for (unsigned long i=0; i<inside_mesh_index.size(); i++) {
        long im = position_buffer.at(i);
        unsigned long j = head_index_ext.at(im) + indexes.at(im);
        assert(sorted_index_ext.at(j) == 0);
        sorted_index_ext.at(j) = inside_mesh_index.at(i);
        indexes.at(im)++;
    }
}



void PairList::search_ext(long im, long jmex, long jmey, const std::vector<Atom> &my_atoms, const std::vector<Atom> &ext_atoms) {
    if (jmex<0) {
        jmex += nmx_ext;
    } else if (jmex>=nmx_ext) {
        jmex -= nmx_ext;
    }
    if (jmey<0) {
        jmey += nmy_ext;
    } else if (jmey>=nmy_ext) {
        jmey -= nmy_ext;
    }
    long jme = jmex + jmey*nmx_ext;
    for (unsigned long n=head_index_ext.at(jme); n<head_index_ext.at(jme)+counter_ext.at(jme); n++) {
        bool survive = false;
        unsigned long j = sorted_index_ext.at(n);
        for (unsigned long i=head_index.at(im); i<head_index.at(im)+counter.at(im); i++) {
            double dx = ext_atoms.at(j).x - my_atoms.at(i).x;
            double dy = ext_atoms.at(j).y - my_atoms.at(i).y;
            periodic_distance(dx, dy);
            double r = std::sqrt(dx*dx + dy*dy);
            if (r > sysp::co_margin)
                continue;
            Pair p;
            set_pair(p, i, j, my_atoms.at(i).id, ext_atoms.at(j).id);
            this->one_other_list.push_back(p);
            survive = true;
        }
        if (survive) {
            this->survivor_list.at(j) = true;
        }
    }
}



// 実行はmesh_search()よりも後でなければならない。
void PairList::mesh_search_ext(const std::vector<Atom> &my_atoms, std::vector<Atom> &ext_atoms) {
    this->clear_ext();
    this->set_mesh_ext();
    this->set_index_ext(ext_atoms);
    survivor_list.resize(ext_atoms.size());
    std::fill(survivor_list.begin(), survivor_list.end(), false);
    if (nmx>2 && nmy>2) {
        for (long i=0; i<num_mesh; i++) {
            long jmex = i%nmx + 1;
            long jmey = i/nmx + 1;
            for (int m=-1; m<2; m++) {
                for (int n=-1; n<2; n++) {
                    this->search_ext(i, jmex+m, jmey+n, my_atoms, ext_atoms);
                }
            }
            
            // 拡張メッシュ領域が周期境界によってつながっている場合を考慮。
            if (i%nmx<2) {
                for (int m=-1; m<2; m++) {
                    this->search_ext(i, jmex-2, jmey+m , my_atoms, ext_atoms);
                    this->search_ext(i, jmex-3, jmey+m , my_atoms, ext_atoms);
                }
            }
            if (i%nmx>nmx-3) {
                for (int m=-1; m<2; m++) {
                    this->search_ext(i, jmex+2, jmey+m , my_atoms, ext_atoms);
                    this->search_ext(i, jmex+3, jmey+m , my_atoms, ext_atoms);
                }
            }
            if (i/nmx<2) {
                for (int m=-1; m<2; m++) {
                    this->search_ext(i, jmex+m, jmey-2, my_atoms, ext_atoms);
                    this->search_ext(i, jmex+m, jmey-3, my_atoms, ext_atoms);
                }
            }
            if (i/nmx>nmy-3) {
                for (int m=-1; m<2; m++) {
                    this->search_ext(i, jmex+m, jmey+2, my_atoms, ext_atoms);
                    this->search_ext(i, jmex+m, jmey+3, my_atoms, ext_atoms);
                }
            }

            if(i==0) {
                this->search_ext(i, nmx_ext-1 , nmy_ext-1, my_atoms, ext_atoms);
                this->search_ext(i, nmx_ext-2 , nmy_ext-1, my_atoms, ext_atoms);
                this->search_ext(i, nmx_ext-1 , nmy_ext-2, my_atoms, ext_atoms);
                this->search_ext(i, nmx_ext-2 , nmy_ext-2, my_atoms, ext_atoms);
            }
            if(i==nmx-1) {
                this->search_ext(i, 0, nmy_ext-1, my_atoms, ext_atoms);
                this->search_ext(i, 0, nmy_ext-2, my_atoms, ext_atoms);
                this->search_ext(i, 1, nmy_ext-1, my_atoms, ext_atoms);
                this->search_ext(i, 1, nmy_ext-2, my_atoms, ext_atoms);
            }
            if(i==num_mesh-nmx) {
                this->search_ext(i, nmx_ext-1, 0, my_atoms, ext_atoms);
                this->search_ext(i, nmx_ext-2, 0, my_atoms, ext_atoms);
                this->search_ext(i, nmx_ext-1, 1, my_atoms, ext_atoms);
                this->search_ext(i, nmx_ext-2, 1, my_atoms, ext_atoms);
            }
            if(i==num_mesh-1) {
                this->search_ext(i, 0, 0, my_atoms, ext_atoms);
                this->search_ext(i, 0, 1, my_atoms, ext_atoms);
                this->search_ext(i, 1, 0, my_atoms, ext_atoms);
                this->search_ext(i, 1, 1, my_atoms, ext_atoms);
            }
        }
    } else {
        search_all_ext(my_atoms, ext_atoms);
    }
}



void PairList::clear_ext(void) {
    nmx_ext = 0;
    nmy_ext = 0;
    num_mesh_ext = 0;
    one_other_list.clear();
    survivor_list.clear();
    counter_ext.clear();
    head_index_ext.clear();
    sorted_index_ext.clear();
}



void PairList::search_all_ext(const std::vector<Atom> &my_atoms, const std::vector<Atom> &ext_atoms) {
    const unsigned long pn = ext_atoms.size();
    for (unsigned long i=0; i<pn; i++) {
        double ix = ext_atoms[i].x;
        double iy = ext_atoms[i].y;
        bool survive = false;
        for (unsigned long j=0; j<my_atoms.size(); j++) {
            double dx = my_atoms[j].x - ix;
            double dy = my_atoms[j].y - iy;
            periodic_distance(dx, dy);
            double r = std::sqrt(dx*dx + dy*dy);
            if (r > sysp::co_margin) {
                continue;
            }
            Pair pair;
            set_pair(pair, j, i, my_atoms[j].id, ext_atoms[i].id);
            one_other_list.push_back(pair);
            survive = true;
        }
        if (survive){
            this->survivor_list.at(i) = true;
        }
    }
}



void PairList::pick_atoms(std::vector<Atom> &ext_atoms, std::vector<unsigned long> &one_recv_list, std::vector<unsigned long> &new_index) {
    std::vector<Atom> new_ext_atoms;
    unsigned long pn = survivor_list.size();
    new_index.resize(pn);
    unsigned long index = 0;
    for (unsigned long i=0; i<pn; i++) {
        if (survivor_list.at(i)) {
            new_ext_atoms.push_back(ext_atoms.at(i));
            one_recv_list.push_back(ext_atoms.at(i).id);
            new_index.at(i) = index;
            index++;
        }
    }
}



void PairList::arrange_pairs_ext(const std::vector<unsigned long> &new_index, unsigned long new_pn) {
    unsigned long ln = one_other_list.size();
    std::vector<Pair> new_list(ln);
    std::vector<unsigned long> pair_counter(new_pn);
    std::fill(pair_counter.begin(), pair_counter.end(), 0);
    for (std::size_t i=0; i<ln; i++) {
        Pair& p = one_other_list.at(i);
        p.j = new_index.at(p.j);
        pair_counter.at(p.j)++;
    }

    std::vector<unsigned long> head_index_pair(new_pn);
    unsigned long sum = 0;
    for (unsigned long i=1; i<new_pn; i++) {
        sum += pair_counter.at(i-1);
        head_index_pair.at(i) = sum;
    }

    std::vector<unsigned long> indexes(new_pn);
    std::fill(indexes.begin(), indexes.end(), 0);
    for (std::size_t i=0; i<ln; i++) {
        Pair p = one_other_list.at(i);
        unsigned long j = head_index_pair.at(p.j) + indexes.at(p.j);
        new_list.at(j) = p;
        indexes.at(p.j)++;
    }

    one_other_list.clear();
    one_other_list.resize(ln);
    std::copy(new_list.begin(), new_list.end(), one_other_list.begin());
}

// ============================================
