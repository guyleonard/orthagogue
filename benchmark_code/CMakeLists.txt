#cmake_minimum_required(VERSION 2.6)
cmake_minimum_required(VERSION 2.6)


#
# Conditonal variables:
#cmakedefine BUILD_MAN_PAGE
#cmakedefine assert_code
#cmakedefine NDEBUG
#option (BUILD_MAN_PAGE "Build a man page for the TurboOrtho Software" ON) 
#option (assert_code "Exectues logical software tests on the TurboOrtho Software" ON) 
#option (assert_code "Exectues logical software tests on the TurboOrtho Software" ON) 
#option (NDEBUG "Exectues in-line software tests on the TurboOrtho Software" ON) 


project(Benchmark)
ADD_DEFINITIONS(
-O2
# "-W -Wall -ansi -g
#
# -std=c++0x 
#-pedantic -Wno-write-strings"
)
include_directories("${PROJECT_BINARY_DIR}")
#add_library(bench_mark	  benchmark.cxx)

#add_library(internal_blast internal_blast.cxx)
add_executable(x_bench_mark benchmark.cxx)


#target_link_libraries(x_bench_mark)