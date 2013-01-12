# The log-builder:
find_library(INTERNAL_BLAST_LIBS
  NAMES internal_blast
  PATHS ${PROJECT_BINARY_DIR}/../internal_blast/
  )
if(NOT INTERNAL_BLAST_LIBS) 
include_directories ("${PROJECT_SOURCE_DIR}/../internal_blast/")
include_directories ("${PROJECT_SOURCE_DIR}/internal_blast/")
#  message("PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")
  find_library(INTERNAL_BLAST_LIBS
    NAMES internal_blast
    PATHS ${PROJECT_BINARY_DIR}/internal_blast/
    )
else() 
#  include_directories ("${PROJECT_SOURCE_DIR}/internal_blast/")
  include_directories ("${PROJECT_SOURCE_DIR}/../internal_blast/")
endif()