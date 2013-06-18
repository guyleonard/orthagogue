CMAKE_MINIMUM_REQUIRED(VERSION 2.6) 

#
# Enable CMakefile to be used both as subpart of of other config file, and run seperately:
get_directory_property(not_the_dad_dir PARENT_DIRECTORY)
if(not_the_dad_dir)
  set(suffix "")
  set(PRODUCE_RUNNABLE 0)
else() # Decleares a lot of stuff:
  set(suffix "../")
  set(VERSION 0.9.7)
  PROJECT(blast_filtering)  
  set(PRODUCE_RUNNABLE 1) 
  set(DEFENITIONS_ADDED 0)
  set(COMPILATION_DEBUG 1)  
  set(AVOID_STRING_WARNING 1)
  include(../config_files/project_properties.cmake)
  include(${PROJECT_BINARY_DIR}/${suffix}config_files/compilespec.cmake)
  configure_file("${PROJECT_BINARY_DIR}/../project_config.h.in" "${PROJECT_BINARY_DIR}/../project_config.h") 
endif()

#
# Defines the other directories to include source from:
set(FOLDERLOC_COMMON "${suffix}blast_common")
set(FOLDERLOC_LOGGBUILDER "${suffix}log_builder")
set(FOLDERLOC_TERMINAL "${suffix}terminal_input")
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}) 
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER})  
include(${PROJECT_BINARY_DIR}/${suffix}config_files/cmph_paths.cmake)
include(${PROJECT_BINARY_DIR}/${suffix}config_files/tbb_paths.cmake)
include(${PROJECT_BINARY_DIR}/${suffix}cmake_scripts/findcmdlinelibs.cmake)
include(${PROJECT_BINARY_DIR}/${suffix}cmake_scripts/findloglibs.cmake)
include(${PROJECT_BINARY_DIR}/${suffix}cmake_scripts/get_blast_common_libs.cmake)

#
# Sets the specific files to be used for this build:
include(files_for_libs.cmake)

#
# Attaches all the files to libraries iot make them accessable:
#foreach(name ${CXXFILES_BLAST_FILTERING})
add_library(blast_filtering_m STATIC ${CXXFILES_BLAST_FILTERING})

foreach(name ${CXXFILES_BLAST_FILTERING_AS_LIBRARY})
  get_filename_component(comp "${name}" NAME_WE)
  add_library(${comp} STATIC ${name})
  add_dependencies(blast_filtering_m ${comp})
endforeach()


# Builds the library paths:
target_link_libraries(blast_filtering_m 
  ${TBB_LIBRARY} 
  ${CMPH_LIBS}
  "${CMPH_LIBS_ARGUMENT}"
  ${LOG_LIBS}
  #"${BLAST_COMMON_LIBS}"
  ${CMD_LIBS}
  ${LIBS_COMMON}
  ) 
target_link_libraries(blast_filtering ${TBB_LIBRARY} 
  ${CMD_LIBS}
  ${CMPH_LIBS}
  ${LOG_LIBS}
  ) 
#
# Builds the exectuable:
add_executable(x_blast_filtering  blast_filtering_main.cxx)
target_link_libraries(x_blast_filtering 
  ${TBB_LIBRARY} 
  blast_filtering_m
  ${LIBS_COMMON}
  ) 

#
# If set prints the status of (important variable) configuration used:
if(PRODUCE_RUNNABLE) 
  include(../config_files/settings_and_status.cmake)
  set(EXECUTABLE_NAME "orthAgogue")
  set(TURBO_LIB_NAME "${EXECUTABLE_NAME}-lib")
endif()

#
# Initiates- and builds the packages:
set(TURBO_FILTER_LIB_NAME "${EXECUTABLE_NAME}-filter-lib")
add_library(${TURBO_FILTER_LIB_NAME} STATIC blast_filtering.cxx blast_filtering.h)
target_link_libraries(${TURBO_FILTER_LIB_NAME} blast_filtering ${TBB_LIBRARY} 
  ${CMD_LIBS}
  ${CMPH_LIBS}
  ${LOG_LIBS}
  ) 
# install(FILES blast_filtering.h DESTINATION "include/${EXECUTABLE_NAME}")
# install(TARGETS blast_filtering DESTINATION "lib/${EXECUTABLE_NAME}" EXPORT "${TURBO_FILTER_LIB_NAME}")
# install(TARGETS ${TURBO_FILTER_LIB_NAME} DESTINATION "lib/${EXECUTABLE_NAME}" EXPORT "${TURBO_FILTER_LIB_NAME}")
# install(EXPORT "${TURBO_FILTER_LIB_NAME}" DESTINATION "lib/${TURBO_FILTER_LIB_NAME}")
if(PRODUCE_RUNNABLE) 
  # install(FILES "../orthAgogue-conf.cmake" DESTINATION "lib/${EXECUTABLE_NAME}")
  include("../orthAgogue_description.cmake")
endif()
