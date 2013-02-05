CMAKE_MINIMUM_REQUIRED(VERSION 2.6) 

#
# Enable CMakefile to be used both as subpart of of other config file, and run seperately:
get_directory_property(not_the_dad_dir PARENT_DIRECTORY)
if(not_the_dad_dir)
  set(PRODUCE_RUNNABLE 0)
  set(suffix "")
else() # Decleares a lot of stuff:
  set(suffix "../")
  set(VERSION 0.9.7)
  PROJECT(blast_common)
  set(PRODUCE_RUNNABLE 1) 
  configure_file("${PROJECT_BINARY_DIR}/../project_config.h.in" "${PROJECT_BINARY_DIR}/../project_config.h")
  set(DEFENITIONS_ADDED 0)
  set(COMPILATION_DEBUG 1)  
  set(AVOID_STRING_WARNING 1)
  include(../config_files/project_properties.cmake)
  include(../config_files/compilespec.cmake)
endif()

#
# Defines the other directories to include source from:
set(FOLDERLOC_LOGGBUILDER "${suffix}log_builder")
set(FOLDERLOC_TERMINAL "${suffix}terminal_input")
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}) 
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER})  
include(${PROJECT_BINARY_DIR}/${suffix}config_files/cmph_paths.cmake)
include(${PROJECT_BINARY_DIR}/${suffix}config_files/tbb_paths.cmake)
include(${PROJECT_BINARY_DIR}/${suffix}cmake_scripts/findcmdlinelibs.cmake)
include(${PROJECT_BINARY_DIR}/${suffix}cmake_scripts/findloglibs.cmake)

#
# Sets the specific files to be used for this build:
include(files_for_libs.cmake)

#
# Attaches all the files to libraries iot make them accessable:
foreach(name ${CXXFILES_BLAST_COMMON})
  get_filename_component(comp "${name}" NAME_WE)
  add_library("${comp}" STATIC "${name}")
endforeach(name)

add_library(blast_common STATIC ${CXXFILES_BLAST_COMMON} ${CXXFILES_BLAST_COMMON_LOCAL})

target_link_libraries(blast_common ${TBB_LIBRARY}
  ${CMPH_LIBS}
  ${LOG_LIBS}
  ${CMD_LIBS}
  ) 
target_link_libraries(blast_memory ${TBB_LIBRARY}
  ${CMD_LIBS}
  ${CMPH_LIBS}
  ${LOG_LIBS}
  ) 

#
# Builds the exectuable:
add_executable(x_blast_common blast_common_main.cxx)
target_link_libraries(x_blast_common
  blast_common ${TBB_LIBRARY} 
  prot_list
  ${CMPH_LIBS}
  ${LOG_LIBS}
  ${CMD_LIBS}
  )

#
# If set prints the status of (important variable) configuration used:
if(PRODUCE_RUNNABLE) 
  include(../config_files/settings_and_status.cmake)
  set(EXECUTABLE_NAME "orthaGogue")
  set(TURBO_LIB_NAME "${EXECUTABLE_NAME}-lib")
endif()

#
# Initiates- and builds the packages:
set(TURBO_MEMORY_LIB_NAME "${EXECUTABLE_NAME}-memory-lib")
add_library(${TURBO_MEMORY_LIB_NAME} STATIC blast_memory.cxx blast_memory.h)
target_link_libraries(${TURBO_MEMORY_LIB_NAME} blast_memory ${TBB_LIBRARY}
  ${CMD_LIBS}
  ${CMPH_LIBS}
  ${LOG_LIBS}
  ) 
# install(FILES blast_memory.h DESTINATION "include/${EXECUTABLE_NAME}")
# install(TARGETS blast_memory DESTINATION "lib/${EXECUTABLE_NAME}" EXPORT "${TURBO_MEMORY_LIB_NAME}")
# install(TARGETS ${TURBO_MEMORY_LIB_NAME} DESTINATION "lib/${EXECUTABLE_NAME}" EXPORT "${TURBO_MEMORY_LIB_NAME}")
# install(EXPORT "${TURBO_MEMORY_LIB_NAME}" DESTINATION "lib/${TURBO_MEMORY_LIB_NAME}")


if(PRODUCE_RUNNABLE) 
  # install(FILES "../orthAgogue-conf.cmake" DESTINATION "lib/${EXECUTABLE_NAME}")
  include("../orthAgogue_description.cmake")
endif()
