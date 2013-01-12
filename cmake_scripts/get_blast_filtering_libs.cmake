# Gets the files made libraries of:
include ("${PROJECT_BINARY_DIR}/${FOLDERLOC_FILTERING}/files_for_libs.cmake") 
# Adds them to the variable.
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_FILTERING}) 
set(LIBS_FILTERING "")
set(FILES_FILTERING "")
foreach(name ${CXXFILES_BLAST_FILTERING})
  get_filename_component(comp "${name}" NAME_WE)
  set(LIBS_FILTERING ${LIBS_FILTERING}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_FILTERING}/lib${comp}.a")
  set(FILES_FILTERING ${FILES_FILTERING}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_FILTERING}/${name}")
endforeach(name)
