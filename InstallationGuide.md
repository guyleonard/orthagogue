## Software dependences ##
  * Intel TBB (http://threadingbuildingblocks.org/)
    * To install on Linux Ubuntu run “sudo apt-get install libtbb-dev”.
  * The perfect hash library cmph (http://cmph.sourceforge.net/).
    * To install on Linux Ubuntu run “sudo apt-get install libcmph-dev”.

## Installation using pre-compiled packages for Linux ##
Follow the instructions provided on the download site  https://code.google.com/p/orthagogue/downloads/list . As listed on the page we provide binaries for two Linux distributions (i.e. "RPM" and "DEB"). If the computations are performed on a cluster, we advise the use of the MPI version, which combines MPI with Intel TBB.

## Building binaries from the source code ##
Additional requirements:
  * CMake (OS independent)
  * make (OS independent)
  * C++ compiler (OS dependent)
    * Tested with the Intel  and gnu g++ compilers. To install the latter on Linux Ubuntu run “sudo apt-get install g++”.
  * OpenMPI compiler (if  the MPI extension is used)
    * Tested with the mpic++ compiler. To install the latter on Linux Ubuntu run “sudo apt-get install mpic++”.

A quick guide for installing locally on a Linux machine:
  1. Clone the project from the Mercurial repository https://code.google.com/p/orthagogue/source/checkout and change to the directory 'orthagogue'
  1. Type ./install.bash in the source code folder, or "./mpi\_install.bash" if MPI is used.
  1. To confirm the success of installation type “./orthAgogue” . You should see a help message explaining the use of the software.