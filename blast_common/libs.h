#ifndef libs_h
#define libs_h
/**
   @file
   @brief Specifies commonly used libraries.
   @ingroup common
   @author Ole Kristian Ekseth (oekseth)
   @date 29.12.2011 by oekseth (initial)
**/ 
#include "../configure.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <map> // the c++ hash map to use
#include <queue> // a FIFO to retrieve the file locations
#include <syscall.h>
#include <unistd.h>
#include <sys/types.h>
#include <cmph.h>
//#include "tbb.h"
#include <list>
// Degubber
//#include "libfence.a"
using namespace std;

#endif
