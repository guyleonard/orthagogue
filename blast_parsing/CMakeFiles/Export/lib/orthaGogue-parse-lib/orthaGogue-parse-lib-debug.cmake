#----------------------------------------------------------------
# Generated CMake target import file for configuration "DEBUG".
#----------------------------------------------------------------

# Commands may need to know the format version.
SET(CMAKE_IMPORT_FILE_VERSION 1)

# Compute the installation prefix relative to this file.
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)

# Import target "blast_parsing" for configuration "DEBUG"
SET_PROPERTY(TARGET blast_parsing APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
SET_TARGET_PROPERTIES(blast_parsing PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "/gpfs/shareapps/apps/modulessoftware/intel/compilers/12.1.0/composer_xe_2011_sp1.6.233/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so;/home/olekrie/bin/cmph/src/.libs/libcmph.so;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libtaxa.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libparse.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/liblist_norm.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libnorm_t.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libprot_list.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libalgo_overlap.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libid_simil_list.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/liblist_file_parse.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libfile_parse.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/librelations_list.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libmcl_format.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libmcl_bunch.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libprotein_relation.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libprotein_vector.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libblast_memory.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libindex_list.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libindex.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/librel.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libp_rel.a;/home/olekrie/orthaGogue_0.9.9.6/src/blast_common/libmpi_read.a;/home/olekrie/orthaGogue_0.9.9.6/src/log_builder/liblog_builder.a;/home/olekrie/bin/cmph/src/.libs/libcmph.so"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/orthaGogue/libblast_parsing.a"
  )

# Import target "orthaGogue-parse-lib" for configuration "DEBUG"
SET_PROPERTY(TARGET orthaGogue-parse-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
SET_TARGET_PROPERTIES(orthaGogue-parse-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "blast_parsing;/gpfs/shareapps/apps/modulessoftware/intel/compilers/12.1.0/composer_xe_2011_sp1.6.233/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so;/home/olekrie/orthaGogue_0.9.9.6/src/terminal_input/libcmd_argument.a;/home/olekrie/orthaGogue_0.9.9.6/src/terminal_input/libcmd_list.a;/home/olekrie/bin/cmph/src/.libs/libcmph.so;/home/olekrie/orthaGogue_0.9.9.6/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/orthaGogue/liborthaGogue-parse-lib.a"
  )

# Cleanup temporary variables.
SET(_IMPORT_PREFIX)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
