#message("ENV{HOME}=$ENV{HOME}")
#message("ENV{HOME}=$ENV{TBBROOT}")
if(NOT DEFINED TBB_LIBRARY OR NOT EXISTS "${TBB_LIBRARY}") # If not already set:  
  find_library(TBB_LIBRARY
    NAMES tbb
    PATHS 
    "$ENV{TBB_PATH_LIB}"
    "/usr/src/redhat/BUILD/tbb-1.1/src/.libs/"
    "$ENV{HOME}/bin/"
    "$ENV{HOME}/bin/tbb/"
    "/usr/local/intel/compilers/12.1.0/tbb/"
    "$ENV{TBBROOT}/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/"
    "/usr/local/intel/compilers/12.1.0/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/"
    )
#message("-- $ENV{TBBROOT}")
#message(    "$ENV{TBBROOT}/lib/intel64/cc4.1.0_libc2.4_kernel2.6.1") 
  find_path(TBB_INCLUDES tbb.h
    $ENV{TBB_PATH_INCLUDES}
    /usr/local/include/tbb/
    /usr/include/tbb/
    /usr/src/redhat/BUILD/tbb-1.1/src/.libs/
    "$ENV{HOME}/bin/src/"
    "$ENV{HOME}/bin/tbb/src/"
    $ENV{TBBROOT}/include/tbb/
    "/usr/local/intel/compilers/12.1.0/tbb/"
    "/usr/local/intel/compilers/12.1.0/tbb/include/tbb/"
    )
  if (TBB_LIBRARY AND TBB_INCLUDES AND EXISTS "${TBB_INCLUDES}/tbb.h")
    include_directories(${TBB_INCLUDES})
#    message("TBB_LIBRARY=${TBB_LIBRARY} AND TBB_INCLUDES=${TBB_INCLUDES}")
  else()
    message("1.\tInternal variables: TBB_LIBRARY=${TBB_LIBRARY} AND TBB_INCLUDES=${TBB_INCLUDES}")
    message("2.\tSystem variables:   TBB_PATH_LIB=${TBB_PATH_LIB} AND TBB_PATH_INCLUDES=${TBB_PATH_INCLUDES}")
    message(FATAL_ERROR "The tbb library (path) not found in your systems library. Please first assure that the tbb library has been installed, try setting the system variables 'TBB_PATH_LIB'=<directory (folder) of your tbb library> and 'TBB_PATH_INCLUDES'=<directory (folder) of your tbb header files>,  and if this does not bring the nuts to the squirrel, then give the developer Ole Kristian Ekseth (oekseth@gmail.com) an update if you found it, but the software didn't, as it would imply you are running an architecture- or update not tested. In sum, a gentle nod to oekseth@gmail.com would be velcomly received.")
  endif()
else() 
endif()




