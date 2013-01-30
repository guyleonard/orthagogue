# Set version information in a config.h file
#set(myproject_VERSION_MAJOR |) 
#set(myproject_VERSION_MINOR 0)
include(InstallRequiredSystemLibraries)

set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
if(NOT CPACK_PACKAGE_VERSION_MAJOR)
  SET(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR_VERSION}")
endif()
if(NOT CPACK_PACKAGE_VERSION_MINOR)
  SET(CPACK_PACKAGE_VERSION_MINOR "${MINOR_VERSION}")
endif()
if(NOT CPACK_PACKAGE_VERSION_PATCH)
  SET(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")
endif()
set(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")


#SET(CPACK_BINARY_RPM "ON")
#SET(CPACK_RPM_COMPONENT_INSTALL "ON" )

# SET(CPACK_BINARY_BUNDLE "ON") # Only for "OsX"
# SET(CPACK_BINARY_ZIP "ON")
# SET(CPACK_SOURCE_ZIP "ON")
message("CMake_VERSION_MAJOR=${CMake_VERSION_MAJOR}")
message("CMake_VERSION_MINOR=${CMake_VERSION_MINOR}")
#TODO: Define the usage of the below variables:
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
SET(CPACK_PACKAGE_EXECUTABLES "${EXECUTABLE_NAME}" "My ${EXECUTABLE_NAME}")
SET(CPACK_PACKAGE_CONTACT "Ole Kristian Ekseth (oekseth@gmail.com)")

# FIXME: Validate the below settings:
SET(CPACK_SYSTEM_NAME "Linux_${CMAKE_SYSTEM}")
SET(CPACK_TOPLEVEL_TAG "Linux-Source_${CMAKE_SYSTEM}")

if(NOT DEFINED PACKAGE)
  set(PACKAGE "DEB")
endif()

message("\t PACKAGE=${PACKAGE}")

if(DEFINED PACKAGE)

  set(CPACK_PACKAGE_NAME "ortAgogue")
  set(CPACK_PACKAGE_VENDOR "Ole Kristian Ekseth <oekseth@gmail.com>, The Systems Biology Group, NTNU, Norway")
#SET(CPACK_PACKAGE_VENDOR  "O.K. Ekseth (email: oekseth@gmail.com), High Performance Computing Group, NTNU, Norway")
  SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ortAgogue - A parallel filtering tool for orthology estimation.")
#  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ortAgogue - orthAgogue: a tool for high speed estimation of homology relations within and between species in massive data sets.")
  set(CPACK_PACKAGE_DESCRIPTION "orthAgogue, a tool for high-speed orthology estimation: The comparison of genes and gene products across species depends on high quality tools to determine the relationships between gene and protein sequences from various species. Although some excellent applications are available and widely used, there is room for further improvement, especially concerning performance and precision. We therefore developed orthAgogue: a tool for high speed estimation of homology relations within and between species in massive data sets. orthAgogue is easy to use and offers flexibility through a range of optional parameters.\n Homepage: https://code.google.com/p/orthagogue/")
  set(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
  set(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
  set(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
  set(CPACK_PACKAGE_INSTALL_DIRECTORY "orthAgogue-${VERSION}")
  SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
  SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
  # set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/ReadMe.txt")
  # set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/Copyright.txt")
  set(CPACK_STRIP_FILES TRUE)
  set(CPACK_TOPLEVEL_TAG "package")
  set(CPACK_PACKAGE_FILE_NAME "orthAgogue-${VERSION}")

  if(PACKAGE STREQUAL "RPM")
    set(CPACK_GENERATOR "RPM")
    if(NOT DEFINED CPACK_RPM_PACKAGE_RELEASE)
      set(CPACK_RPM_PACKAGE_RELEASE "1.fc12")
    endif(NOT DEFINED CPACK_RPM_PACKAGE_RELEASE)
    set(CPACK_RPM_PACKAGE_LICENSE "GPLv2+")
    #! Dependencies:
    set(CPACK_RPM_PACKAGE_REQUIRES "linux-vdso")
    set(CPACK_RPM_PACKAGE_REQUIRES "libcmph")
    set(CPACK_RPM_PACKAGE_REQUIRES "libtbb")
    set(CPACK_RPM_PACKAGE_REQUIRES "libcmph")
    set(CPACK_RPM_PACKAGE_REQUIRES "libstdc++")
    set(CPACK_RPM_PACKAGE_REQUIRES "libm")
    set(CPACK_RPM_PACKAGE_REQUIRES "libgcc_s")
    set(CPACK_RPM_PACKAGE_REQUIRES "libc")
    set(CPACK_RPM_PACKAGE_REQUIRES "libdl")
    set(CPACK_RPM_PACKAGE_REQUIRES "libpthread")

    set(CPACK_PACKAGE_FILE_NAME "ortAgogue-${VERSION}-${CPACK_RPM_PACKAGE_RELEASE}")
  endif(PACKAGE STREQUAL "RPM")

  if(PACKAGE STREQUAL "DEB") # Debian pacakage generation:
    SET(CPACK_BINARY_DEB "ON")

# -------------------
    set(CPACK_GENERATOR "DEB")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Ole Kristian Ekseth <oekseth@gmail.com>")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.4), libgcc1 (>= 1:4.1.1), libstdc++6 (>= 4.1.1) libcmph0 (>=0.9) libtbb2 (>=4.0)")
#    set(CPACK_DEBIAN_PACKAGE_SECTION "liblo7")
#SET(CPACK_DEBIAN_PACKAGE_DEPENDS "kdebase-runtime (>= 4:4.2.1), kdelibs5 (>= 4:4.2.1), libc6 (>= 2.1.3), libgcc1 (>= 1:4.1.1), libplasma3, libqt4-dbus (>= 4.5.0), libqtcore4 (>= 4.5.0), libqtgui4 (>= 4.5.0), libstdc++6 (>= 4.2.1)")
    set(CPACK_DEBIAN_PACKAGE_VERSION "${VERSION}-ubuntu")
  endif(PACKAGE STREQUAL "DEB")

  set(CPACK_SOURCE_GENERATOR "ZIP")
#  set(CPACK_SOURCE_GENERATOR "TGZ")
  set(CPACK_SOURCE_PACKAGE_FILE_NAME "orthAgogue-${VERSION}" CACHE INTERNAL "tarball basename")
  set(CPACK_IGNORE_FILES ";/\\\\.git/;/\\\\.hg/;CMakeLists.txt.user;/build/;/website/;/reference/;/artwork/;/heuristical_orthology_building/;/control_set/;/*.a/")
  set(CPACK_SOURCE_IGNORE_FILES ${CPACK_IGNORE_FILES})

  include(CPack)

endif(DEFINED PACKAGE)

#INCLUDE(CPack) 