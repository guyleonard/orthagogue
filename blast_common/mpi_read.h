#ifndef mpi_file_read_h
#define mpi_file_read_h
/**
   @file
   @brief Maps a protein-id to a set of other proteins (i.e. identifies the protein pairs).
   @ingroup blastfile_container
   @date 04.01.2012 by oekseth (cleanup).
   @author Ole Kristian Ekseth (oekseth).
**/ 
//#include "libs.h"
#include "types.h"
#include "log_builder.h"
/*
 * Copyright 2012 Ole Kristian Ekseth (oekseth@gmail.com)
 *
 * This file is part of orthAgogue.
 *
 * orthAgogue is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * orthAgogue is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with orthAgogue
. If not, see <http://www.gnu.org/licenses/>.
 */

/**
   @class mpi_read
   @brief Procedures for reading a file using MPI procedures.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth).
   @date 05.07.2012 by oekseth (initial).
**/ 
class mpi_read {
#ifdef USE_MPI
 private:
  int myrank;
  MPI_File fh;
  MPI_Datatype filetype; 
  bool file_is_set;
  loint reading_file_start_position;
  loint reading_file_end_position;
  loint reading_file_length;
  loint block_length;
  //! @return the file chunk:
  char *get_file_chunk(loint &string_length);
 public:  
/**
   @brief Reads the data
   @param <str>    The char buffer to put things into.
   @param <str_len>  The maximum number of elements to read.
   @return the number of elements read.
**/
  loint read_chars(char *string, size_t string_length);

  //! @returns a chunk of a block:
  char *read_a_block(loint &string_length) {
    return get_file_chunk(string_length);
  }
//! @return true if there is more data to read.
  bool has_more_data_in_file_to_read();
   
 static void close(mpi_read *&obj) {
   if(obj) {obj->free_mem(), delete obj, obj = NULL;}
 }
  //! Deallocates the memory.
 void free_mem();

  //! The constructor:
  mpi_read(MPI_Comm comm, int _myrank, char *file_name, loint _reading_file_start_position, loint _reading_file_end_position, loint _reading_file_length, loint _read_block_length);
  

  //! The constructor:
 mpi_read() : 
   myrank(0), file_is_set(false), reading_file_start_position(0), reading_file_end_position(0), reading_file_length(0), block_length(0)
   {}
#endif
};


/**
   @class mpi_read
   @brief Procedures for reading a file using MPI procedures.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth).
   @date 05.07.2012 by oekseth (initial).
**/ 
//typedef class mpi_read mpi_read_t;
#endif
