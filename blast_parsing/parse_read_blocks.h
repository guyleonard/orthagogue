/**
   @file
   @brief Assures correct reading of blast files
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
**/
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
#ifndef parse_reead_blocks_h
#define parse_reead_blocks_h
#include "string_section.h"
#include "thread_reading.h"
#include "../configure.h"
#include "log_builder.h"
/**
   @class parse_read_blocks
   @brief Assures correct reading of blast files
   @ingroup parsing_container
   @remarks The file is read twise: Collecting break-point in the first run,
   enables efficeint reading in the second run. This implies the the class is
   used during the second read in order to ensure the correct order of the
   blocks written to the buffer for the current taxon by avoiding the danger of
   two similar protein vectors overlapping each other
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 15.09.2011 by oekseth (asserts)
   @date 19.11.2011 by oekseth (each p-thread reading their from their own file)
   @date 03.12.2011 by oekseth (interaction with class 'buffer_string_list' caused some modifications.)
**/
class parse_read_blocks {
 private:
  uint disk_buffer_size;
 public:
  /**
     @brief Holds file positions.
     @remarks Understanding this is of importance if the procedure fo file reading in class pipe_parse_parse is to be understood.
     - Each value in this list refers to a specific position in the file used producing the data
     - The length of each block is found taking arr[index]-arr[index-1];
  **/
  loint *arr; 
  //! The size set for the list.
  loint size;
  //! The last inserted position in the array.
  loint position_in_arr;
  //! During the process of reading, the next index position to retrieve.
  loint read_postion_to_retrieve;
 private:
  thread_reading_t *thread_list; // the list of threads used for specific reading
  uint n_threads;
  loint last_file_pos_in_mem;
  //! Allocates and initializes the read-buffer to be sent downward the pipe:
  void allocate(char *&buffer_main, const long int buffer_in_mem_size);
  //! Allocates and initializes the read-buffer to be sent downward the pipe:
  void allocate(char *&buffer_main, char *&buffer_end, const long int buffer_in_mem_size);

  //! Reads from the given file-pointer
  void read_file(char *file_name, long int pos_in_file, char *buffer_main, const long int buffer_in_mem_size);

  /**
     @brief Reads from the given file-pointer
     @param <file_input> The file pointer tor ead from.
     @param <buffer_main> The reserved char array to read the data into.
     @param <buffer_in_mem_size> The variable to get number of chars to read.
     @return true if the there are no more chars to read.
     @author Ole Kristian Ekseth (oekseth)
  **/
  bool read_file(FILE *file_input, char *buffer_main, const long int buffer_in_mem_size);
 public:
  /**
     @brief A possible start position for the reading
     @param <numb> The last position in the file stored in local memory.
     @remarks Task to know where to read the first chars from when reading the second run starts, i.e. where to move the file pointer before reading starts. In the context where the first n chars are stored in memory, and the rest read from the file.
   **/
  void update_last_file_pos_im_mem(loint numb) {
    if (numb > last_file_pos_in_mem) last_file_pos_in_mem = numb;
  }  
  //! Asserts that the block ends with a new-line, and if not prints the debug-info
  void assert_newline_end(char *buffer_main, long int size);
  //! Gets the block using fseek:
  char *get_block(char *file_name, long int pos_in_file, long int size);
  /**
     @brief Inserts the data inside the structure:
     @remarks Size of each block will (on average) be proportional to the size of the data used for input.
  **/
  void insert(loint file_pos);
  //! Returns the section of the data using the parse-blocks as a basis:
  string_section *get_string_section_from_file(FILE *file_input);
  //! Writing of the blocks in memory to stdout.
  void printBlocks();
  /**
     @return The length of the block, incrementing the poiner for the next read.
     @remarks '0' if no more blocks to read
  **/
  long int getNextReadLength();

  /**
     @brief Generate (writes) a log file describing the details of the memory signature for this 'unit'
     @remarks Useful for analyzing- optimizing the memory allocation procedures.
  **/
  void log_generate_memory_allocation_overview(uint n_threads);

/**
   @brief Initiates the file reading, using the data set. Thread safe.
   @param <_n_threads> The number of threads that weill be calling this function.
   @param <file_name> The name of the file to be read. (Should be the same file as block data are stored for.) 
   @remarks
   - To be called when the file is read- and the 'blocks' stored in the list named 'arr' are consistent.
   - Sets the index for the positions to retrieve in the objects of type thread_reading.
   @author Ole Kristian Ekseth (oekseth)
**/
  void init_thread_file_reading(uint _n_threads, char *file_name);
  //! Returns true if all the data in the fiels are read:
  bool all_data_read(bool print_result);
  //! @return True if all data in the file is read.
  bool assert_file_is_read() {return all_data_read(false);} 
  //! Returns the section of the data using the parse-blocks as a basis:
  string_section *get_string_section(uint thread_id);
  //! De-allocates memory stored for this object.
  void free_mem();
  //! The constructor.
  parse_read_blocks(uint _disk_buffer_size);

#ifdef assert_code
  //! Asserts the thread reading.
  void assert_thread_file_reading();
  //! Asserts the private functions of this class.
  void assert_private_parts();
#endif
  //! The main test function for this class.
  static void assert_class(const bool print_info);
};

/**
   @brief Assures correct reading of blast files
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
   @date 03.12.2011 by oekseth (interaction with class 'buffer_string_list' caused some modifications.)
**/
typedef class parse_read_blocks parse_read_blocks_t;

#endif
