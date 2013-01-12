#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
SET(CMAKE_IMPORT_FILE_VERSION 1)

# Compute the installation prefix relative to this file.
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)

# Import target "blast_parsing" for configuration "Release"
SET_PROPERTY(TARGET blast_parsing APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(blast_parsing PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/usr/local/intel/compilers/12.1.0/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so;/usr/local/lib/libcmph.so;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libtaxa.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libparse.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/liblist_norm.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libnorm_t.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libprot_list.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libalgo_overlap.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libid_simil_list.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/liblist_file_parse.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libfile_parse.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/librelations_list.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libmcl_format.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libmcl_bunch.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libprotein_relation.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libprotein_vector.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libblast_memory.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libindex_list.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libindex.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/librel.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libp_rel.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/blast_common/libmpi_read.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/log_builder/liblog_builder.a;/usr/local/lib/libcmph.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthaGogue/libblast_parsing.a"
  )

# Import target "orthaGogue-parse-lib" for configuration "Release"
SET_PROPERTY(TARGET orthaGogue-parse-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(orthaGogue-parse-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "blast_parsing;/usr/local/intel/compilers/12.1.0/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/terminal_input/libcmd_argument.a;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/terminal_input/libcmd_list.a;/usr/local/lib/libcmph.so;/norstore/user/olekrie/orthAgogue_1.0.0.0/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthaGogue/liborthaGogue-parse-lib.a"
  )

# Cleanup temporary variables.
SET(_IMPORT_PREFIX)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
