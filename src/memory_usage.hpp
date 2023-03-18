#pragma once

// ==========================================

namespace memory_usage {
		extern struct rusage r_usage;
		extern unsigned long max_mem;
		
		void check(void);
		unsigned long get_result(void);
		void reset(void);
}

// ==========================================
