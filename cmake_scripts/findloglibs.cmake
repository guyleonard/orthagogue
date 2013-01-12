# The log-builder:
find_library(LOG_LIBS
  NAMES log_builder
  PATHS ${PROJECT_BINARY_DIR}/../log_builder/
  )
if(NOT LOG_LIBS) 
  set(LOG_LIBS "")
  include_directories ("${PROJECT_SOURCE_DIR}/../log_builder/")
#  message("PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")
  #message("LOG_LIBS = ${LOG_LIBS}")
  find_library(LOG_LIBS
    NAMES log_builder
    PATHS ${PROJECT_BINARY_DIR}/log_builder/
    )
  if(NOT LOG_LIBS) 
    set(LOG_LIBS "")
#message("tries including:\n${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}/files_for_libs.cmake")
    include ("${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}/files_for_libs.cmake") 
    # Adds them to the variable.
#    include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}) 
    foreach(name ${CXXFILES_LOGGBUILDER})
      get_filename_component(comp "${name}" NAME_WE)
      set(LOG_LIBS ${LOG_LIBS}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}/lib${comp}.a")
    endforeach(name)
#    message("i findloglibs: LOG_LIBS=${LOG_LIBS}...\n")
  endif()
else() 
  include_directories ("${PROJECT_SOURCE_DIR}/log_builder")
endif()