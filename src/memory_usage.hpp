#pragma once

// ==========================================

class MemoryUsage {
    private:
				struct rusage r_usage;
				unsigned long max_mem;
    
    public:
        void check(void);
        unsigned long get_result(void);
        void reset(void);
};

// ==========================================
