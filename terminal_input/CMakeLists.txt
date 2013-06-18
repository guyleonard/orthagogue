CMAKE_MINIMUM_REQUIRED(VERSION 2.6) 

#
# Enable CMakefile to be used both as subpart of of other config file, and run seperately:
get_directory_property(not_the_dad_dir PARENT_DIRECTORY)
if(not_the_dad_dir)
  set(PRODUCE_RUNNABLE 0)
  include(../config_files/compilespec.cmake)
  set(suffix "")
else() # Decleares a lot of stuff:
  set(suffix "../")
  set(VERSION 0.9.7)
  PROJECT(cmd_line)
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
set(FOLDERLOC_COMMON "${suffix}blast_common")
set(FOLDERLOC_LOGGBUILDER "${suffix}log_builder")
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER})  
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_COMMON}) 
include(${PROJECT_BINARY_DIR}/${suffix}cmake_scripts/get_blast_common_libs.cmake)
include(../config_files/cmph_paths.cmake)

#
# Sets the specific files to be used for this build:
include(files_for_libs.cmake)

add_library(blast_terminal cmd_list.cxx cmd_argument.cxx)


#
# Attaches all the files to libraries iot make them accessable:
foreach(name ${CXXFILES_BLAST_TERMINAL})
  get_filename_component(comp "${name}" NAME_WE)
  add_library(${comp} ${name})
  add_dependencies(blast_terminal ${comp})
endforeach()

add_library(cmd_line cmd_list.cxx cmd_argument.cxx)
add_dependencies(cmd_line blast_terminal)

target_link_libraries(cmd_line ${TBB_LIBRARY})
target_link_libraries(cmd_argument ${TBB_LIBRARY})

#
# Builds the exectuable:
add_executable(x_terminal main_cmd_line.cxx)
target_link_libraries(x_terminal cmd_line) 

#
# If set prints the status of (important variable) configuration used:
if(PRODUCE_RUNNABLE) 
  include(../config_files/settings_and_status.cmake)
endif()