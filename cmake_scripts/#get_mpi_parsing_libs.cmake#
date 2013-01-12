# Gets the files made libraries of:
include ("${PROJECT_BINARY_DIR}/${FOLDERLOC_PARSING}/files_for_libs.cmake") 
# Adds them to the variable.
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_PARSING}) 
set(LIBS_PARSING "")
set(FILES_PARSING "")
foreach(name ${CXXFILES_BLAST_PARSING})
  get_filename_component(comp "${name}" NAME_WE)
  set(LIBS_PARSING ${LIBS_PARSING}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_PARSING}/lib${comp}.a")
  set(FILES_PARSING ${FILES_PARSING}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_PARSING}/${name}")
endforeach(name)

