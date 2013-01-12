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
   @file
   @brief Read a block of data from the blast file, producing an object of tyoe string_section consisting of a char-block with extra meta-data.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
**/
#ifndef read_file_h
#define read_file_h
#include "../configure.h"
#include "log_builder.h"
#include "mpi_read.h"
#include "parse_read_blocks.h"
#include "string_section.h"
/**
   @class read_file
   @brief Read a block of data from the blast file, producing an object consisting of a char-block with extra meta-data.
   @ingroup parsing_ops
   @remarks Task is procing an oject of type string_section".
   -# Validates reading of data by generating an empty test-set with a predefined set of newline-chars.
   -# Uses either MPI og ordinary c-routines for reading, dependent of the constructuro used initiating this class.
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 25.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a library)
   @date 05.07.2012 by oekseth (added MPI functionality and improved of the internals asserts in current methods.)
**/
class read_file {
 private:
  uint disk_buffer_size;
  int myrank;
  mpi_read *mpi_obj;

  uint buffer_rest_total_size;
  char *BUFFER_REST; 
  char *buffer_rest;
  // Gives the current location in the rest buffer
  long int buffer_rest_size;
  mem_loc chars_sendt_to_parsing;

  //! Allocates and initializes the read-buffer to be sent downward the pipe:
  void allocate(char *&buffer_main, char *&buffer_end, const long int buffer_in_mem_size);
  //! Reads the file:
  void read_from_file(char *buffer_main, const long int buffer_in_mem_size);  
  //! Reads the file:
  bool read_from_file(char *buffer_main, long int buffer_in_mem_size, long int buffer_rest_size, long int &size_current_buffer);
  //! Finding the last new_line char:
  char *get_logical_end(char *buffer_main,char *buffer_main_end, const bool end_of_file, long int &buffer_rest_size);
  //! Copies the owerflowing data from the previous buffer to this.
  void copy_overflow(char *buffer_main, char *buffer_rest, long int buffer_rest_size);
  //! Updates the buffer with data to be included on the next run.
  void update_overflow_buffer(char *buffer_main_start, char *buffer_rest, char *buffer_rest_start, char *buffer_main_logical_end, long int buffer_rest_size,
 long int buffer_rest_total_size, long int size_current_buffer);
//! Verifies the correctness of the block
  void verify_correctness_of_block(int myrank, char *start, char *buffer_end, char seperator, const int line_caller, const char *file_caller);
  //! Prints the data in the buffer:
  void print_buffer(char *buffer_main, char *buffer_end);
public:
  //! Pointer to the file.
  FILE *file_input;
  /**
     Reads the file and dump a block of chars to the next element in the pipe
     @return An object holding the chars read.
  **/
  struct string_section *read_file_blocks(char seperator, uint &block_cnt, const bool perform_tests);
  //! Positions the file-pointer at the start of the file.
  void rewind_file() {rewind(file_input);}
  /**
     @brief Sets the position of the file pointer
     @return true if operation was possible
  **/
  bool set_file_pointer_position(loint file_position);
  //! Gets the real position in the file using ftell
  long ftell_get_real_position_in_file(){
    if(file_input) return ftell(file_input);
    else return 0;
  }
  /**
     @brief Remove content of overflow buffer.
     @remarks Useful when jumpring around in a file, only retrieving the content, ie, not merging it     
  **/
  void set_overflow_buffer_to_emtpy();

  //! @return the section of the data using the parse-blocks as a basis:
  string_section *get_section();
  /**
     @brief De-allocates the memory.
     @remark To be called after the operation is finalized
  **/
  void free_mem();
  //! Returns an initiated object of this class.
  static read_file *init(uint _disk_buffer_size, char *f_name) {return new read_file(_disk_buffer_size, f_name);}
#ifdef USE_MPI
  //! Returns an initiated object of this class using MPI procedures
  static read_file *init(uint _disk_buffer_size, MPI_Comm comm, int myrank, char *f_name, loint reading_file_start_position, loint reading_file_end_position, loint reading_file_length) {
    return new read_file(_disk_buffer_size, comm, myrank, f_name, reading_file_start_position, reading_file_end_position, reading_file_length);
  }
#endif

  //! De-allocates memory of the object given as input.
  static void close(read_file *&tmp) {
    if(tmp) tmp->free_mem(), delete tmp, tmp = NULL;
  }
  //! The method of parallisation.
  /*overrride*/void* operator()(void* item);
  //! The constructor.
  read_file(uint disk_buffer_size, char *f_name);

#ifdef USE_MPI
  //! The constructor for MPI-procedure-initiation.
  read_file(uint disk_buffer_size, MPI_Comm comm,int myrank, char *f_name, loint reading_file_start_position, loint reading_file_end_position, loint reading_file_length);
#endif

  //! The test function for the private parts of this class.
  void assert_private_parts();
  //! The main test function for this class  
  static void assert_class(const bool print_info);
};

/**
   @brief Read a block of data from the blast file, producing an object consisting of a char-block with extra meta-data.
   @ingroup parsing_ops
   @author: Ole Kristian Ekseth (oekseth)
   @date 25.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a library)
**/
typedef class read_file read_file_t;
#endif
