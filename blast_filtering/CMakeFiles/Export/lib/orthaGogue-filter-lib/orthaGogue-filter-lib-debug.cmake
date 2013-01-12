#----------------------------------------------------------------
# Generated CMake target import file for configuration "DEBUG".
#----------------------------------------------------------------

# Commands may need to know the format version.
SET(CMAKE_IMPORT_FILE_VERSION 1)

# Compute the installation prefix relative to this file.
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)

# Import target "blast_filtering" for configuration "DEBUG"
SET_PROPERTY(TARGET blast_filtering APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
SET_TARGET_PROPERTIES(blast_filtering PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "/gpfs/shareapps/apps/modulessoftware/intel/compilers/12.1.0/composer_xe_2011_sp1.6.233/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so;/home/olekrie/orthaGogue_0.9.9.6/src/terminal_input/libcmd_argument.a;/home/olekrie/orthaGogue_0.9.9.6/src/terminal_input/libcmd_list.a;/home/olekrie/bin/cmph/src/.libs/libcmph.so;/home/olekrie/orthaGogue_0.9.9.6/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/orthaGogue/libblast_filtering.a"
  )

# Import target "orthaGogue-filter-lib" for configuration "DEBUG"
SET_PROPERTY(TARGET orthaGogue-filter-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
SET_TARGET_PROPERTIES(orthaGogue-filter-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "blast_filtering;/gpfs/shareapps/apps/modulessoftware/intel/compilers/12.1.0/composer_xe_2011_sp1.6.233/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so;/home/olekrie/orthaGogue_0.9.9.6/src/terminal_input/libcmd_argument.a;/home/olekrie/orthaGogue_0.9.9.6/src/terminal_input/libcmd_list.a;/home/olekrie/bin/cmph/src/.libs/libcmph.so;/home/olekrie/orthaGogue_0.9.9.6/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/orthaGogue/liborthaGogue-filter-lib.a"
  )

# Cleanup temporary variables.
SET(_IMPORT_PREFIX)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
