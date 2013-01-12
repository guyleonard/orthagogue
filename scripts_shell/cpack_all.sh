#!/bin/sh
#------------------------------------------------------------
# THE CPACK CONFIGURATION chmod
#------------------------------------------------------------
#echo "Runs the cpack configuration."

# FIXME: inkluder nedenfor!

# cmake .;
# # Binary package:
cpack -C CPackConfig.cmake;
# cpack -C CPackConfig.cmake -G DEB;
# cpack -C CPackConfig.cmake -G NSIS;
# cpack -C CPackConfig.cmake -G RPM;
# # Source package:
# cpack -C CPackSourceConfig.cmake;
