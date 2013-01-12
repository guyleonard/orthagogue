# The log-builder:
find_library(CMD_LIBS_ARGUMENT
  NAMES cmd_argument
  PATHS ${PROJECT_BINARY_DIR}/../terminal_input/
  )
find_library(CMD_LIBS
  NAMES cmd_line
  PATHS ${PROJECT_BINARY_DIR}/../terminal_input
  )
if(NOT CMD_LIBS) 
  set(CMD_LIBS "")
  include_directories ("${PROJECT_SOURCE_DIR}/../terminal_input")
# message("Not found at path: ${PROJECT_BINARY_DIR}/../terminal_input")
#   message("PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")
#   message("LOG_LIBS = ${LOG_LIBS}")
  find_library(CMD_LIBS
    NAMES cmd_line
    PATHS ${PROJECT_BINARY_DIR}/terminal_input/
    )
  find_library(CMD_LIBS_ARGUMENT
    NAMES cmd_argument
    PATHS ${PROJECT_BINARY_DIR}/terminal_input/
    )
  if(NOT CMD_LIBS) 
    set(CMD_LIBS "")
    include_directories("${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}")
    include ("${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}/files_for_libs.cmake") 
    # Adds them to the variable.
    include_directories (${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}) 
    foreach(name ${CXXFILES_BLAST_TERMINAL})
      get_filename_component(comp "${name}" NAME_WE)
      set(CMD_LIBS ${CMD_LIBS}  "${PROJECT_BINARY_DIR}/${FOLDERLOC_TERMINAL}/lib${comp}.a")
    endforeach(name)
  endif()
else() 
  include_directories ("${PROJECT_SOURCE_DIR}/terminal_input")
endif()


# if(NOT ${CMD_LIBS})
#   message(".....an error\n")
# else()
#   message(".....ok..found CMD_LIBS\n")
# endif()
