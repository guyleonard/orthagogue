#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
SET(CMAKE_IMPORT_FILE_VERSION 1)

# Compute the installation prefix relative to this file.
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)

# Import target "blast_filtering" for configuration "Release"
SET_PROPERTY(TARGET blast_filtering APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(blast_filtering PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/usr/local/intel/compilers/12.1.0/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/terminal_input/libcmd_argument.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/terminal_input/libcmd_list.a;/usr/local/lib/libcmph.so;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthaGogue/libblast_filtering.a"
  )

# Import target "orthaGogue-filter-lib" for configuration "Release"
SET_PROPERTY(TARGET orthaGogue-filter-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(orthaGogue-filter-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "blast_filtering;/usr/local/intel/compilers/12.1.0/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/terminal_input/libcmd_argument.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/terminal_input/libcmd_list.a;/usr/local/lib/libcmph.so;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthaGogue/liborthaGogue-filter-lib.a"
  )

# Cleanup temporary variables.
SET(_IMPORT_PREFIX)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
