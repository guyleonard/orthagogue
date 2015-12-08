The comparison of genes and gene products across species depends on high quality tools to determine the relationships between gene and protein sequences from various species:
  * Although some excellent applications are available and widely used, there is room for further improvement, especially concerning performance and precision.
  * We therefore developed orthAgogue: a tool for high speed estimation of homology relations within and between species in massive data sets.
  * orthAgogue is easy to use and offers flexibility through a range of optional parameters.

## Getting Started ##
After downloading orthAgogue, unpack it, read the documentation wiki pages (listed on the right side of this front page). The downloadables are:
  1. Debian for large clusters (ombines Open MPI with Intels TBB): http://orthagogue.googlecode.com/hg/orthAgogue-1.0.3.mpi.deb
  1. Debian (e.g. Ubuntu): http://orthagogue.googlecode.com/hg/orthAgogue-1.0.3.deb
  1. RPM for large clusters (combines Open MPI with Intels TBB): http://orthagogue.googlecode.com/hg/orthAgogue-1.0.3.mpi.rpm
  1. RPM: http://orthagogue.googlecode.com/hg/orthAgogue-1.0.3.rpm

The packages are tested on Ubuntu 12.04.1 LTS, using the amd64 architecture. The minimum-software requirements are: libc6 (>= 2.4), libgcc1 (>= 1:4.1.1), libstdc++6 (>= 4.1.1), libcmph-dev (>=0.9), libtbb-dev (>=4.0), mpic++ (>=4.6.3). To install the packages, follow the 'classical appraoch'. For Debian (e.g. Ubuntu), write in tour terminal/bash: sudo dpkg -i orthAgogue-1.0.2.deb

**Note:** The recent version (1.0.3) solves an important issue: sometimes Blast manages to give a higher overlap-score to different proteins than self-comparison. To exemplify, for two arbitrary proteins "A" and "B", the overlap(A, B) (assigned by Blast) are sometimes greater than the maximum overlap of overlap(A,A) or overlap(B,B). In earlier versions, these cases have been asserted to be both rare and insignificant, though thanks to the feedback from a user, the issue was identified. In brief, if your calculations does not match your expectations (which is an unlikely case using orthAgogue, and therefore of importance), please do not hesitate contacting us, either through the 'issue' tab, or directly to the main-developer (oekseth@gmail.com).

## Who are using orthAgogue? ##

The Semantic-Systems Biology Group at NTNU, Norway.

## What are the dependencies? ##

All dependencies are freely available under the GNU public licence.
  1. Parallelization is enabled through Intel TBB and OPEN-MPI.
  1. Text-hashing is made effective by the cmph hashing library.
  1. Building- and compilation through the Cmake and g++ tool chains.

The dependencies may either be installed locally (in $ENV{HOME}/bin/) or at the default locations (if installed using root).