#ifndef tbb_libs_h
#define tbb_libs_h
/**
   @file
   @brief Specifies tbb libraries.
   @ingroup common
   @author Ole Kristian Ekseth (oekseth)
   @date 26.06.2012 by oekseth (initial)
**/ 
#include "concurrent_hash_map.h"
#include "concurrent_queue.h"
#include "concurrent_vector.h"
#include "pipeline.h"
#include "tick_count.h"
#include "task_scheduler_init.h"
#include "tbb_allocator.h"
#include "tick_count.h"
#include "pipeline.h"
#include "parallel_for.h"
#include "task_scheduler_init.h"
#include "atomic.h" // for atomic varaibles
#include "spin_mutex.h"
/* #include "concurrent_hash_map.h" */
/* #include "concurrent_queue.h" */
/* #include "concurrent_vector.h" */
/* #include "pipeline.h" */
/* #include "tick_count.h" */
/* #include "task_scheduler_init.h" */
/* #include "tbb_allocator.h" */
/* #include "tick_count.h" */
/* #include "pipeline.h" */
/* #include "parallel_for.h" */
/* #include "task_scheduler_init.h" */
/* #include "atomic.h" // for atomic varaibles */
/* #include "spin_mutex.h" */
using namespace tbb; 
//! Usage: Holds the co-orthologs for every protein in a stack
typedef concurrent_queue<o_rel_t> stack_rel; 
//! A spin_lock, used when setting- and releasing an unique id in the pipes.
typedef spin_mutex slock_t; 

#endif
