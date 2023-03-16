#include "memory_usage.hpp"

// ==========================================

void MemoryUsage::check(void) {
    getrusage(RUSAGE_SELF, &r_usage);
    unsigned long mem_usage = r_usage.ru_maxrss;
		std::cout << mem_usage << std::endl;
		if (max_mem < mem_usage) {
				max_mem = mem_usage;
		}
}



unsigned long MemoryUsage::get_result(void) {
    return max_mem;
}



void MemoryUsage::reset(void) {
    max_mem = 0;
}

// ==========================================
