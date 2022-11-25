#pragma once

// ==========================================

class CalcTimer {
    private:
        std::chrono::system_clock::time_point start_time, end;
        double result;
        bool is_measuring = false;
    
    public:
        void start(void);
        void stop(void);
        double get_result(void);
        void reset(void);
        void share(double*);
};

// ==========================================