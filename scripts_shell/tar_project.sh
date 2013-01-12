#!/bin/sh

echo ""
echo "Runs the cpack config."
# Binary package:
cpack -C CPackConfig.cmake
# Source package:
cpack -C CPackSourceConfig.cmake

exit;