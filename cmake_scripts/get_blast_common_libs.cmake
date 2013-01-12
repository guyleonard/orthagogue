# Gets the files made libraries of:
include ("${PROJECT_BINARY_DIR}/${FOLDERLOC_COMMON}/files_for_libs.cmake") 
# Adds them to the variable.
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_COMMON}) 
set(LIBS_COMMON  "")
set(FILES_COMMON "")
foreach(name ${CXXFILES_BLAST_COMMON})
  get_filename_component(comp "${name}" NAME_WE)
  set(LIBS_COMMON  ${LIBS_COMMON}   "${PROJECT_BINARY_DIR}/${FOLDERLOC_COMMON}/lib${comp}.a")
  set(FILES_COMMON ${FILES_COMMON}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_COMMON}/${name}")
endforeach(name)