# Gets the files made libraries of:
include ("${PROJECT_BINARY_DIR}/${FOLDERLOC_MPI}/files_for_libs.cmake") 
# Adds them to the variable.
include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_MPI}) 
set(LIBS_MPI "")
set(FILES_MPI "")
foreach(name ${CXXFILES_BLAST_MPI})
  get_filename_component(comp "${name}" NAME_WE)
  set(LIBS_MPI ${LIBS_MPI}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_MPI}/lib${comp}.a")
  set(FILES_MPI ${FILES_MPI}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_MPI}/${name}")
endforeach(name)
