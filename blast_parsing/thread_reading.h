/**
   @file
   @brief  Reads a file for å specific thread.
   @ingroup parsing_ops
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
#ifndef thread_reading_h
#define thread_reading_h
#include "string_section.h"
#include "../configure.h"
/**
   @class thread_reading
   @brief  Reads a file for å specific thread.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 19.11.2011 by oekseth (initial)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/
class thread_reading {
 private:
  uint disk_buffer_size;
  char *file_name;
  FILE *file_input;
  // Below values refers to indexes used in the list found in class 'parse_read_blocks':
  uint block_index_start;
  uint block_index_end;
  uint current_index; // in range [block_index_start, block_index_end-1]
  loint global_file_reading_start; // USed if its the first file pointer to get data.
 public:
  //! @return the location of where in the file the object will start its reading.
  uint get_index_start(){return block_index_start;}
  //! @return the location of where in the file the object will end its reading.
  uint get_index_end(){return block_index_end;}
  //! Prints information about the setting regarding the file reading for this object.
  void print_block();
  //! @return the total length of the data to be read by this object.
  loint get_total_length(loint *arr);
  //! @return true if the file is read
  bool file_is_read(uint my_id, bool print_data_);
 private:
  //! @return the length of the block given: assumes the indexes ar in the range of the arr given.
  loint get_block_length(loint *arr, uint start_index, uint end_index);

  //! @return the starting position in the file for the given index:
  loint get_file_pos(loint *arr, uint index);
  //! Allocates and initializes the read-buffer to be sent downward the pipe:
  void allocate(char *&buffer_main, char *&buffer_end, const loint buffer_in_mem_size);
  /**! @return the length of the block, incrementing the poiner for the next read
     @Return: '0' if no more blocks to read
  */
  loint getNextReadLength(loint *arr);
  /**
     @brief Reads the file using the position given as input.
     @param <file_position_array>  The array where the position data is located.
     @param <buffer_main> The reserved char array to read the data into.
     @param <buffer_in_mem_size> The variable to get the length of the chars inserted into the char parameter.
     @return true if the there are no more chars to read.
     @author Ole Kristian Ekseth (oekseth)
  **/
  bool read_file(loint *file_position_array, char *buffer_main, loint &buffer_in_mem_size);
 public:
  /**
     @brief Reads a block from the input file using positions in input.
     @param <my_id> The identifier specifying the if of the file pointer to use.
     @param <file_position_array>  The array where the position data is located.
     @return the section of the data using the parse-blocks as a basis:
     @remarks Is thread-safe. Steps taken are:
     - Get number of chars to read (length), summing up several 'lengths" if they in total constitute a sum of bytes less that the variable disk_buffer_size.
     - On the first call moves the file pointer to the point were this object starts the read. (For an example of this usage, see class pipe_parse_merge in combination with class buffer_string_list.
     @author Ole Kristian Ekseth (oekseth)
  **/
  string_section *get_string_section(uint my_id, loint *file_position_array);

  //! @return the section of the data using the parse-blocks as a basis:
  string_section *get_string_section(loint *file_position_array);
  //! De-allocates memory for this object.
  void free_mem() {if(file_input) fclose(file_input);}
  //! De-allocates memory for the collection of objects, given as input.
  static void free_mem(thread_reading *list, uint size);
  //! The constructor.
  thread_reading(uint _disk_buffer_size, char *_file_name, uint _block_index_start, uint _block_index_end, loint _global);
  //! The constructor.
  thread_reading();
  //! Assert functions for the class.
  static void assert_class(const bool print_info);
};

/**
   @brief  Reads a file for å specific thread.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/
typedef class thread_reading thread_reading_t;
#endif
