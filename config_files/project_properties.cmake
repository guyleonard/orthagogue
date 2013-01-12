if(NOT PROJECT_PROP_SET) # User section
  MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
  # allow more human readable "if then else" constructs
  SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )  

  #SET(CMAKE_BUILD_TYPE distribution)
  #SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
  #SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")
 
  #                                        #
  # GLOBAL VARIABLES FOR THE CONFIG FILE:  #
  #                                        #
  # The project version number:
  SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
  SET( ${PROJECT_NAME}_MINOR_VERSION 1 )
  SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )
  add_definitions(-DVERSION="${VERSION}")
  set(PROJECT_PROP_SET 1)  

  ## DISK BLOCK SIZES ##
  MATH(EXPR def_size_file_read "1024 * 1024 * 5")
  set(BLOCK_FILE_READ_SIZE ${def_size_file_read} CACHE STRING "The number of bytes to read from the disk on each call to fread(..)")
  MATH(EXPR def_size "1024 * 1024 * 20")
  set(DISK_BUFFER_SIZE ${def_size}  CACHE STRING "The size in Bytes to send from one pipe to another. Value would affect time- and memory usage.") # The default size

  ## DOCUMENATION ##
  set(BUILD_LOG_FILES 1 CACHE STRING "Builds the log files documenting the running time.")
  set(BUILD_MAN_PAGE 0 CACHE STRING "Builds the man page for the softwares usage")
  set(LOG_FOLDER_NAME "report_orthaGogue" CACHE STRING "Location storing the log files in.")

  ## DEBUG PARAMETERS ##
  set(NDEBUG 1 CACHE STRING "Includes extra parameters to the software for debugging purposes (named DEUBG_*).")
  set(INCLUDE_CMD_DEBUG_PARAMS 1 CACHE STRING "")
  set(MEMORY_CONSUMPTION_LEVEL 1 CACHE STRING "If set uses the maximises the speed by using most of the avilable memory on your system. By default set to 1")
  set(MAN_PATH_USER 0 CACHE STRING "") # If set stores the manpath describing the usage of this tool at this location. Else the root dir of the cmake installation is used.

  
  #  set(assert_code 1 CACHE STRING "Set if asserts of the functions is to be compiled in.")
  
  #                                            #
  # GLOBAL VARIABLES ACCESSABLE FOR THE USER:  #
  #                                            #
  
  option(USE_MPI "Use several nodes wrapping the processes with 'mpic++' " OFF)
  ## DISK BLOCK SIZES ##
  option(BLOCK_FILE_READ_SIZE "The number of bytes to read from the disk on each call to fread(..)" ${def_size_file_read})
  option(DISK_BUFFER_SIZE  "The size in Bytes to send from one pipe to another. Value would affect time- and memory usage." ${def_size})

  ## DOCUMENATION ##
  option(BUILD_LOG_FILES, "Builds the log files documenting the running time." ON)
  option(BUILD_MAN_PAGE, "Builds the man page for the softwares usage." OFF)

  ## DEBUG PARAMETERS ##
  option(assert_code, "Set if asserts of the functions is to be compiled in." ON)
  option(NDEBUG, "Enables the extra asserts to be run during execution: practical if software behaves in an unexpected manner." ON)
  option(INCLUDE_CMD_DEBUG_PARAMS "Includes extra parameters to the software for debugging purposes (named DEUBG_*)." ON)
  option(MEMORY_CONSUMPTION_LEVEL
   "If set uses the maximises the speed by using most of the avilable memory on your system. By default set to 1" 1)
 option(LOG_FOLDER_NAME "Specififies the location to store the generated logfiles." report_orthaGogue)
#  option(MAN_PATH_USER     "If set stores the manpath describing the usage of this tool at this location. Else the root dir of the cmake installation is used. Default set to 0" 0)

#message("in project_properties the DISK_BUFFER_SIZE = ${DISK_BUFFER_SIZE}")
endif()