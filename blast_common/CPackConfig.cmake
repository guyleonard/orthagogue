# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. Example variables are:
#   CPACK_GENERATOR                     - Generator used to create package
#   CPACK_INSTALL_CMAKE_PROJECTS        - For each project (path, name, component)
#   CPACK_CMAKE_GENERATOR               - CMake Generator used for the projects
#   CPACK_INSTALL_COMMANDS              - Extra commands to install components
#   CPACK_INSTALLED_DIRECTORIES           - Extra directories to install
#   CPACK_PACKAGE_DESCRIPTION_FILE      - Description file for the package
#   CPACK_PACKAGE_DESCRIPTION_SUMMARY   - Summary of the package
#   CPACK_PACKAGE_EXECUTABLES           - List of pairs of executables and labels
#   CPACK_PACKAGE_FILE_NAME             - Name of the package generated
#   CPACK_PACKAGE_ICON                  - Icon used for the package
#   CPACK_PACKAGE_INSTALL_DIRECTORY     - Name of directory for the installer
#   CPACK_PACKAGE_NAME                  - Package project name
#   CPACK_PACKAGE_VENDOR                - Package project vendor
#   CPACK_PACKAGE_VERSION               - Package project version
#   CPACK_PACKAGE_VERSION_MAJOR         - Package project version (major)
#   CPACK_PACKAGE_VERSION_MINOR         - Package project version (minor)
#   CPACK_PACKAGE_VERSION_PATCH         - Package project version (patch)

# There are certain generator specific ones

# NSIS Generator:
#   CPACK_PACKAGE_INSTALL_REGISTRY_KEY  - Name of the registry key for the installer
#   CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS - Extra commands used during uninstall
#   CPACK_NSIS_EXTRA_INSTALL_COMMANDS   - Extra commands used during install


SET(CPACK_BINARY_BUNDLE "")
SET(CPACK_BINARY_CYGWIN "")
SET(CPACK_BINARY_DEB "ON")
SET(CPACK_BINARY_DRAGNDROP "")
SET(CPACK_BINARY_NSIS "")
SET(CPACK_BINARY_OSXX11 "")
SET(CPACK_BINARY_PACKAGEMAKER "")
SET(CPACK_BINARY_RPM "")
SET(CPACK_BINARY_STGZ "")
SET(CPACK_BINARY_TBZ2 "")
SET(CPACK_BINARY_TGZ "")
SET(CPACK_BINARY_TZ "")
SET(CPACK_BINARY_ZIP "")
SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")
SET(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
SET(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.4), libgcc1 (>= 1:4.1.1), libstdc++6 (>= 4.1.1) libcmph0 (>=0.9) libtbb2 (>=4.0)")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Ole Kristian Ekseth <oekseth@gmail.com>")
SET(CPACK_DEBIAN_PACKAGE_VERSION "1.0.0-ubuntu")
SET(CPACK_GENERATOR "DEB")
SET(CPACK_IGNORE_FILES ";/\\.git/;/\\.hg/;CMakeLists.txt.user;/build/;/website/;/reference/;/artwork/;/heuristical_orthology_building/;/control_set/;/*.a/")
SET(CPACK_INSTALL_CMAKE_PROJECTS "/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common;blast_common;ALL;/")
SET(CPACK_INSTALL_PREFIX "/usr/local")
SET(CPACK_MODULE_PATH "")
SET(CPACK_NSIS_DISPLAY_NAME "orthAgogue-1.0.0")
SET(CPACK_NSIS_INSTALLER_ICON_CODE "")
SET(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
SET(CPACK_NSIS_PACKAGE_NAME "orthAgogue-1.0.0")
SET(CPACK_OUTPUT_CONFIG_FILE "/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/CPackConfig.cmake")
SET(CPACK_PACKAGE_CONTACT "Ole Kristian Ekseth (oekseth@gmail.com)")
SET(CPACK_PACKAGE_DEFAULT_LOCATION "/")
SET(CPACK_PACKAGE_DESCRIPTION "orthAgogue, a tool for high-speed orthology estimation: The comparison of genes and gene products across species depends on high quality tools to determine the relationships between gene and protein sequences from various species. Although some excellent applications are available and widely used, there is room for further improvement, especially concerning performance and precision. We therefore developed orthAgogue: a tool for high speed estimation of homology relations within and between species in massive data sets. orthAgogue is easy to use and offers flexibility through a range of optional parameters.
 Homepage: https://code.google.com/p/orthagogue/")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/ReadMe.txt")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ortAgogue - A parallel filtering tool for orthology estimation.")
SET(CPACK_PACKAGE_EXECUTABLES "orthaGogue;My orthaGogue")
SET(CPACK_PACKAGE_FILE_NAME "orthAgogue-1.0.0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "orthAgogue-1.0.0")
SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "orthAgogue-1.0.0")
SET(CPACK_PACKAGE_NAME "ortAgogue")
SET(CPACK_PACKAGE_RELOCATABLE "true")
SET(CPACK_PACKAGE_VENDOR "Ole Kristian Ekseth <oekseth@gmail.com>, The Systems Biology Group, NTNU, Norway")
SET(CPACK_PACKAGE_VERSION "1.0.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_RESOURCE_FILE_LICENSE "/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/Copyright.txt")
SET(CPACK_RESOURCE_FILE_README "/usr/share/cmake-2.8/Templates/CPack.GenericDescription.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "/usr/share/cmake-2.8/Templates/CPack.GenericWelcome.txt")
SET(CPACK_SET_DESTDIR "OFF")
SET(CPACK_SOURCE_CYGWIN "")
SET(CPACK_SOURCE_GENERATOR "ZIP")
SET(CPACK_SOURCE_IGNORE_FILES "/\\.git/;/\\.hg/;CMakeLists.txt.user;/build/;/website/;/reference/;/artwork/;/heuristical_orthology_building/;/control_set/;/*.a/")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/klatremus/Dokumenter/Work/code/orthAgogue/src/blast_common/CPackSourceConfig.cmake")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "orthAgogue-1.0.0")
SET(CPACK_SOURCE_TBZ2 "")
SET(CPACK_SOURCE_TGZ "")
SET(CPACK_SOURCE_TZ "")
SET(CPACK_SOURCE_ZIP "")
SET(CPACK_STRIP_FILES "TRUE")
SET(CPACK_SYSTEM_NAME "Linux_Linux-3.2.0-36-generic")
SET(CPACK_TOPLEVEL_TAG "package")
