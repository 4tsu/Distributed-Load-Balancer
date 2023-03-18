#include "memory_usage.hpp"

// ==========================================

namespace memory_usage {
		struct rusage r_usage;
		unsigned long max_mem = 0;
		
		void check(void) {
				getrusage(RUSAGE_SELF, &r_usage);
				unsigned long mem_usage = r_usage.ru_maxrss;
				if (max_mem < mem_usage) {
						max_mem = mem_usage;
				}
		}

		unsigned long get_result(void) {
				return max_mem;
		}

		void reset(void) {
				max_mem = 0;
		}
}

// ==========================================
