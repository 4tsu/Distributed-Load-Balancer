#include "calctimer.hpp"

// ==========================================

void CalcTimer::start(void) {
    assert(!this->is_measuring);
    this->start_time = std::chrono::system_clock::now();
    this->is_measuring = true;
}



void CalcTimer::stop(void) {
    assert(this->is_measuring);
    end = std::chrono::system_clock::now();
    const double rap = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - this->start_time).count() / 1000.0 / 1000.0);
    this->result += rap;
    this->is_measuring = false;
}



double CalcTimer::get_result(void) {
    return result;
}



void CalcTimer::reset(void) {
    this->result = 0;
}



void CalcTimer::share(double* results){
    assert(!this->is_measuring);
    MPI_Gather(&this->result, 1, MPI_DOUBLE, results, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
} 
// ==========================================
