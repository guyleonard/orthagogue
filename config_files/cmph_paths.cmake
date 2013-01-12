#message("ENV{HOME}=$ENV{HOME}")
if(NOT DEFINED CMPH_LIBS) # If not already set:
  find_library(CMPH_LIBS
    NAMES cmph
    PATHS
    "$ENV{CMPH_PATH_LIB}"
    /usr/local/lib/
    /usr/src/redhat/BUILD/cmph-1.1/src/.libs/
    "$ENV{HOME}/bin/src/.libs/"
    "$ENV{HOME}/bin/cmph/src/.libs/"
    )
  find_path(CMPH_INCLUDES cmph.h
    "$ENV{CMPH_PATH_INCLUDES}"
    /usr/local/include
    /usr/include
    /usr/src/redhat/BUILD/cmph-1.1/src/.libs/
    "$ENV{HOME}/bin/src/"
    "$ENV{HOME}/bin/cmph/src/"
    )
  if (CMPH_LIBS AND CMPH_INCLUDES AND EXISTS "${CMPH_INCLUDES}/cmph.h")
    include_directories(${CMPH_INCLUDES})
#    message("CMPH_LIBS=${CMPH_LIBS} AND CMPH_INCLUDES=${CMPH_INCLUDES}")
  else()
    message("1.\tInternal variables: CMPH_LIBS=${CMPH_LIBS} AND CMPH_INCLUDES=${CMPH_INCLUDES}")
    message("2.\tSystem variables:   CMPH_PATH_LIB=${CMPH_PATH_LIB} AND CMPH_PATH_INCLUDES=${CMPH_PATH_INCLUDES}")
    message(FATAL_ERROR "The cmph library (path) not found in your systems library. Please first assure that the cmph library has been installed, try setting the system variables 'CMPH_PATH_LIB'=<directory (folder) of your cmph library> and 'CMPH_PATH_INCLUDES'=<directory (folder) of your cmph header files>,  and if this does not bring the nuts to the squirrel, then give the developer Ole Kristian Ekseth (oekseth@gmail.com) an update if you found it, but the software didn't, as it would imply you are running an architecture- or update not tested. In sum, a gentle nod to oekseth@gmail.com would be velcomly received.")
  endif()
else() 
endif()




