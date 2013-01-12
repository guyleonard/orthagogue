# Gets the files made libraries of:
include ("${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}/files_for_libs.cmake") 
# Adds them to the variable.
include_directories("${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}")
set(CMD_LIBS "")
set(LIBS_CMD "")
set(FILES_CMD "")
foreach(name ${CXXFILES_BLAST_TERMINAL})
  get_filename_component(comp "${name}" NAME_WE)
  set(LIBS_CMD ${LIBS_CMD}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}/lib${comp}.a")
  set(CMD_LIBS ${CMD_LIBS}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}/lib${comp}.a")
  set(FILES_CMD ${FILES_CMD} "${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}/${name}")
#  message("${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}/${name}")
endforeach(name)

if(NOT CMD_LIBS) 
  message("an error in the script building the cmd libraries- and absolute paths.")
endif()
