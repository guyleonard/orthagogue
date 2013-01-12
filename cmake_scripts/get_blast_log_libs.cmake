# Gets the files made libraries of:
include ("${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}/files_for_libs.cmake") 
# Adds them to the variable.
include_directories("${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}")
set(LOG_LIBS "")
set(LIBS_LOG "")
set(FILES_LOG "")
foreach(name ${CXXFILES_LOGGBUILDER})
  get_filename_component(comp "${name}" NAME_WE)
  set(LOG_LIBS ${LOG_LIBS}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}/lib${comp}.a")
  set(LIBS_LOG ${LIBS_LOG}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}/lib${comp}.a")
  set(FILES_LOG ${FILES_LOG} "${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}/${name}")
#  message("${PROJECT_BINARY_DIR}/${FOLDERLOC_LOGGBUILDER}/${name}")
endforeach(name)

if(NOT LOG_LIBS) 
  message("an error in the script building the log_builder libraries- and absolute paths.")
endif()
