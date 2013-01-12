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
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/usr/lib/libtbb.so;/usr/local/lib/libcmph.so;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libtaxa.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libparse.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/liblist_norm.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libnorm_t.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libprot_list.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libalgo_overlap.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libid_simil_list.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/liblist_file_parse.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libfile_parse.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/librelations_list.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libmcl_format.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libmcl_bunch.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/liblist_file_chunk.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libprotein_relation.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libprotein_vector.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libblast_memory.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libindex_list.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libindex.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/librel.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libp_rel.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/libmpi_read.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/log_builder/liblog_builder.a;/usr/local/lib/libcmph.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthAgogue/libblast_parsing.a"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS blast_parsing )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_blast_parsing "${_IMPORT_PREFIX}/lib/orthAgogue/libblast_parsing.a" )

# Import target "orthAgogue-parse-lib" for configuration "Release"
SET_PROPERTY(TARGET orthAgogue-parse-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(orthAgogue-parse-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "blast_parsing;/usr/lib/libtbb.so;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/terminal_input/libcmd_argument.a;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/terminal_input/libcmd_list.a;/usr/local/lib/libcmph.so;/home/klatremus/Dokumenter/Work/code/orthAgogue/src/log_builder/liblog_builder.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/orthAgogue/liborthAgogue-parse-lib.a"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS orthAgogue-parse-lib )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_orthAgogue-parse-lib "${_IMPORT_PREFIX}/lib/orthAgogue/liborthAgogue-parse-lib.a" )

# Loop over all imported files and verify that they actually exist
FOREACH(target ${_IMPORT_CHECK_TARGETS} )
  FOREACH(file ${_IMPORT_CHECK_FILES_FOR_${target}} )
    IF(NOT EXISTS "${file}" )
      MESSAGE(FATAL_ERROR "The imported target \"${target}\" references the file
   \"${file}\"
but this file does not exist.  Possible reasons include:
* The file was deleted, renamed, or moved to another location.
* An install or uninstall procedure did not complete successfully.
* The installation package was faulty and contained
   \"${CMAKE_CURRENT_LIST_FILE}\"
but not all the files it references.
")
    ENDIF()
  ENDFOREACH()
  UNSET(_IMPORT_CHECK_FILES_FOR_${target})
ENDFOREACH()
UNSET(_IMPORT_CHECK_TARGETS)

# Cleanup temporary variables.
SET(_IMPORT_PREFIX)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
