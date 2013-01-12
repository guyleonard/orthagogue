if(NOT DEFINED TURBO_LIBS) # If not already set:
  find_library(TURBO_LIBS
    NAMES internal
    PATHS ${PROJECT_BINARY_DIR}../internal_blast/
    )
  find_path(TURBO_INCLUDES internal_blast.h
#    PATHS ${PROJECT_BINARY_DIR}../internal_blast/
    PATHS ../internal_blast/
    )
  if (TURBO_LIBS AND TURBO_INCLUDES AND EXISTS "${TURBO_INCLUDES}/cmph.h")
    include_directories(${TURBO_INCLUDES})
  else()
  endif()
else() 
endif()

