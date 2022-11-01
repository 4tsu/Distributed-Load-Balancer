#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>



void read(std::string filename) {
    std::ifstream reading_file;
    reading_file.open(filename, std::ios::in);
    std::string line;
    bool is_step = false;
    bool is_init = false;
    bool is_num = false;
    int is_bounds = 0;
    bool is_atoms = false;
    int begin_step = 0;
    unsigned int N;
    double x_max, x_min, y_max, y_min;
    while(std::getline(reading_file, line)) {
        // LAMMPS出力ではじめに出てくる、初期配置のdumpであれば無視する。
        if (line=="ITEM: TIMESTEP") {
            is_step = true;
            is_init = false;
            continue;
        } else if (is_step) {
            is_step = false;
            if (line=="0") {
                is_init = true;
            } else {
                begin_step = std::stoi(line);
            }
            continue;
        }
        if (is_init) {
            continue;
        }
        
        // ユーザーが定義したdump、もしくはLAMMPSの最終結果のdumpを読み込む
        if (line=="ITEM: NUMBER OF ATOMS") {
            is_num = true;
            continue;
        } else if (is_num) {
            is_num = false;
            N = std::stoul(line);
            continue;
        }

        if (std::equal(line.begin(), line.begin()+16, "ITEM: BOX BOUNDS")) {
            is_bounds = true;
            continue;
        } else if (std::equal(line.begin(), line.begin()+10, "ITEM: ATOM")) {
            is_bounds = 0;
            is_atoms = true;
            continue;
        } else if (is_bounds) {
            std::string var1, var2;
            auto space_pos = line.find(" ", 0);
            auto length = line.size();
            for (int i=0; i<space_pos; i++)
                var1 += line[i];
            for (int i=space_pos; i<length; i++)
                var2 += line[i];
            if (is_bounds==1) {
                x_min = std::stod(var1);
                x_max = std::stod(var2);
            } else if (is_bounds==2) {
                y_min = std::stod(var1);
                y_max = std::stod(var2);
            }
            is_bounds++;
            continue;
        }

        unsigned int id;
        double x, y, vx, vy;
        std::vector<std::string> vars;
        auto offset = std::string::size_type(0);
        while (true) {
            auto pos = line.find(" ", offset);
            if (pos==std::string::npos) {
                vars.push_back(line.substr(offset));
                break;
            }
            vars.push_back(line.substr(offset, pos-offset));
            offset = pos + 1;
        }
        id = std::stoi(vars.at(0));
        x  = std::stod(vars.at(1));
        y  = std::stod(vars.at(2));
        vx = std::stod(vars.at(3));
        vy = std::stod(vars.at(4));
        std::cout << id << " ";
        std::cout << x << " ";
        std::cout << y << " ";
        std::cout << vx << " ";
        std::cout << vy << std::endl;
    }
    std::cout << begin_step << std::endl;
    std::cout << N << std::endl;
    std::cout << x_min << ", " << x_max << std::endl;
    std::cout << y_min << ", " << y_max << std::endl;
}

int main(int argc, char *argv[]){
    read("droplet.dump");
    return 0;
}