CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#
# Enable CMakefile to be used both as subpart of of other config file, and run seperately:
get_directory_property(not_the_dad_dir PARENT_DIRECTORY)
if(not_the_dad_dir)
  set(PRODUCE_RUNNABLE 0)
  set(suffix "")
  #include(../config_files/cmph_paths.cmake)
else() # Decleares a lot of stuff:
  set(VERSION 0.7.0)
  set(suffix "../")
  PROJECT(logBuilder)
  set(PRODUCE_RUNNABLE 1)
  # configure a header file to pass some of the CMake settings to the source code
  configure_file("${PROJECT_BINARY_DIR}/../project_config.h.in" "${PROJECT_BINARY_DIR}/../project_config.h")
  set(AVOID_STRING_WARNING 1)
  set(COMPILATION_DEBUG 1)  
  include(../config_files/project_properties.cmake)
  include(../config_files/compilespec.cmake)
endif()

#
# Defines the other directories to include source from:
set(FOLDERLOC_COMMON "${suffix}blast_common")
set(FOLDERLOC_TERMINAL "${suffix}terminal_input")
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}) 
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_COMMON}) 
include_directories (${PROJECT_BINARY_DIR}/${PATH_SUFFIX}) 
include(${PROJECT_BINARY_DIR}/${suffix}config_files/cmph_paths.cmake)
include(${PROJECT_BINARY_DIR}/${suffix}cmake_scripts/findcmdlinelibs.cmake)
include(${PROJECT_BINARY_DIR}/${suffix}config_files/tbb_paths.cmake)

#
# Sets the specific files to be used for this build:
include(files_for_libs.cmake)

#
# Attaches all the files to libraries iot make them accessable:
foreach(name ${CXXFILES_LOGGBUILDER})
  get_filename_component(comp "${name}" NAME_WE)
  add_library(${comp} ${name})
  target_link_libraries(${comp} ${CMD_LIBS} ${TBB_LIBRARY})
endforeach()

#
# Builds the exectuable:
add_executable(x_log_builder log_builder_main.cxx)
target_link_libraries(x_log_builder log_builder  ${TBB_LIBRARY} ${CMD_LIBS})

#
# If set prints the status of (important variable) configuration used:
if(PRODUCE_RUNNABLE) 
  include(../config_files/settings_and_status.cmake)
endif()