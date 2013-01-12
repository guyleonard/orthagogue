--------------------------------------------------
         R E A D M E  
---------------------------------------------------

Installation depends on the type you downloaded. E.g. unpacking the
folders for the types *.tar.gz file is: orthAgogue-.*Linux.tar.gz  | tar xvf -

Different releases are found at: http://folk.ntnu.no/olekrie/TurboOrtho/release_turbo/

Dependencies of the code is: "intel tbb" and "cmph" found in the Linux
repository.

Documentation of the code is located at: http://folk.ntnu.no/olekrie/TurboOrtho/Documentation/html/index.html

Usage of the code is presented when the software is run without any
params: ./orthaGogue (or "orthaGogue" directly if installed using
the root).

If the developement version were used in the installation, specific
updates using the cmake configuration tools is easy. Examples of CMake
parameters set from terminal: 
cmake . -DMEMORY_CONSUMPTION_LEVEL=1 -Dassert_code=1
-DBUILD_LOG_FILES=1 -DDISK_BUFFER_SIZE=1024  
---
cmake -DMEMORY_CONSUMPTION_LEVEL=0 . 1>out_smake.txt; make
1>out_make.txt; \ 
valgrind --leak-check=full ./orthaGogue -i all.blast -s '_' -p 0 \
-t 1 -c 1 1>out.txt;


Questions to be forwarded to Ole Kristian Ekseth: oekseth@gmail.com
