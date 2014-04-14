# Finds the correct include path:



get_directory_property(not_the_dad_dir PARENT_DIRECTORY)
set(CMAKE_INCLUDE_PATH     "$ENV{HOME}/bin/")

if(${USE_MPI} STREQUAL "1")
    SET(CMAKE_CXX_COMPILER mpic++)
endif()

if("${USE_MPI}" STREQUAL ON)
  SET(CMAKE_CXX_COMPILER mpic++)
endif()

 
# Owerrites settings if other specifications given on input
if("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
  set(COMPILATION_DEBUG 1)
  set(NDEBUG OFF)
else()	
  # Enable releasesymbols by default
  #message("is_set:\tCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")	
  set(CMAKE_BUILD_TYPE Release)
  set(COMPILATION_DEBUG 0)  
  set(NDEBUG ON) # Discards the 'asserts' in the code.
  #  message("in compilerspec, type set to Release")
endif()

if(NOT "${not_the_dad_dir}" STREQUAL "")
  set(PATH_ONE "${PROJECT_BINARY_DIR}/config_files/CheckCXXCompilerFlag.cmake")
  if(EXISTS "${PATH_ONE}")
    include("${PATH_ONE}")
  else()
    set(PATH_ONE "${PROJECT_BINARY_DIR}/../config_files/CheckCXXCompilerFlag.cmake")
    if(EXISTS "${PATH_ONE}")
      include("${PATH_ONE}")
    endif()
  endif()
else()
  set(PATH_ONE "${PROJECT_BINARY_DIR}/config_files/CheckCXXCompilerFlag.cmake")
  if(EXISTS "${PATH_ONE}")
    include("${PATH_ONE}")
  else()
    set(PATH_ONE "${PROJECT_BINARY_DIR}/../config_files/CheckCXXCompilerFlag.cmake")
    if(EXISTS "${PATH_ONE}")
      include("${PATH_ONE}")
    endif()
  endif()
endif(NOT "${not_the_dad_dir}" STREQUAL "")

# Sets the compiler-defentions:

if(NOT DEFENITIONS_ADDED) 
  CHECK_CXX_COMPILER_FLAG(-Wno-write-strings HAS_NO_WRITE_STRINGS_FLAG)    
  #
  # Intiiates the list, to be used prinitng the settings:
  set(COMPILESPEC_SET  "")  

  if(COMPILATION_DEBUG) 
    #
    # Sets some specific openmpi debug variables:
    if(CMAKE_BUILD_TYPE STREQUAL DEBUG)
      if(USE_MPI)
	if("${USE_MPI}" STREQUAL ON)
	  CHECK_CXX_COMPILER_FLAG(--enable-debug ENABLE_DEBUG)
	  CHECK_CXX_COMPILER_FLAG(--enable-memchecker ENABLE_MEMCHECK)
	  # TODO: Path below is hardcoded: If preferable, make it dynamic, or use an NEV-path.
	  CHECK_CXX_COMPILER_FLAG(--with-valgrind=/usr/bin/valgrind PATH_VALGRIND) 
	  if(ENABLE_DEBUG)
            add_definitions(--enable-debug)
	    set(COMPILESPEC_SET  --enable-debug " " ${COMPILESPEC_SET})
	  endif()
	  if(ENABLE_MEMCHECK)
            add_definitions(--enable-memchecker)
	    set(COMPILESPEC_SET  --enable-memchecker " " ${COMPILESPEC_SET})
	  endif()

	  if(PATH_VALGRIND)
            add_definitions(--with-valgrind=/usr/bin/valgrind)
	    set(COMPILESPEC_SET  --with-valgrind=/usr/bin/valgrind " " ${COMPILESPEC_SET})
	  endif()
	  
	endif()
      endif()
    endif()

    #
    # Sets some standard variables the the ordinary c++-compiler:_
    CHECK_CXX_COMPILER_FLAG(-g HAS_DEBUG_FLAG)
    CHECK_CXX_COMPILER_FLAG(-Wall HAS_WALL_FLAG)    
    CHECK_CXX_COMPILER_FLAG(-ansi HAS_ANSI_FLAG)    

#    message("in compilespec_set '${CMAKE_BUILD_TYPE}'")
    if(CMAKE_BUILD_TYPE STREQUAL DEBUG)
      if(HAS_DEBUG_FLAG)
        add_definitions(-g)
	set(COMPILESPEC_SET  -g " " ${COMPILESPEC_SET})
      endif()

      if(HAS_WALL_FLAG)
        add_definitions(-Wall)
	set(COMPILESPEC_SET  -Wall " " ${COMPILESPEC_SET})
      endif()

      if(HAS_ANSI_FLAG)
        add_definitions(-ansi)
	set(COMPILESPEC_SET  -ansi " " ${COMPILESPEC_SET})
      endif()

      if(HAS_CXX0X_FLAG)
        add_definitions(-std=c++0x)
	set(COMPILESPEC_SET  -std=c++0x " " ${COMPILESPEC_SET})
      endif()

      if(HAS_PEDANTIC_FLAG)
        add_definitions(-pedantic)
	set(COMPILESPEC_SET  -pedantic " " ${COMPILESPEC_SET})
      endif()

      #    add_definitions(-D__DEBUG__)
    else()
      # For future options.
    endif()
  endif()
  #
  # General properties
  if(HAS_NO_WRITE_STRINGS_FLAG)
    add_definitions(-Wno-write-strings)
    set(COMPILESPEC_SET  -Wno-write-strings " " ${COMPILESPEC_SET})
    message("--\t added -Wno-write-string") 
  else()
  endif()
  set(DEFENITIONS_ADDED 1)
endif()

