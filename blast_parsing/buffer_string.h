/**
   @file
   @brief Holds a block of data read from a text file.
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
   @date 02.11.2011 by oekseth (initial)
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
#ifndef buffer_string_h
#define buffer_string_h
#include "types.h"
#include "../configure.h"
#include "log_builder.h"
/**
   @class buffer_string
   @brief Holds a block of data read from a text file.
   @ingroup parsing_container
   @remark The data block contains (a) the owerflow from the previous process and (b) the block of interest from 'this' process, with the index given by the block number sendt from 'pipe_parse_file'
   @author Ole Kristian Ekseth (oekseth)
   @date 02.11.2011 by oekseth (initial)
**/
class buffer_string {
 private:
  char *first_buffer;
  loint size_first_buffer;
  char *second_buffer;
  loint size_second_buffer;
  void print_segment(char *m_start, loint length);
 public:
  //! @return the combined size of both blocks:
  loint get_size(){return (size_first_buffer + size_second_buffer);}
  //! @return the size of the first block
  loint get_size_first_buffer() const {return size_first_buffer;}
  //! @return the size of the second block
  loint get_size_second_buffer() const {return size_second_buffer;}
  //! Prints info about the object.
  void print_info();
  //! @return true if second char buffer is set.
  bool has_first_block() {return (NULL != first_buffer);}
  //! @return true if second char buffer is set.
  bool has_main_block() {return (NULL != second_buffer);}
  //! @return true if eiter first- or second char buffer is set.
  bool has_data() {return (first_buffer || second_buffer);}
  /**
     @param <buff> The data to copy into the first buffer.
     @param <size> The length of the data to copy into the first buffer.
   **/
  void copy_into_first_buffer(char *buff, loint size);

  /**
     @brief Sets the second buffer with the params given.
     @remark Uses the given pointer, without any copying.
  **/
  void set_second_buffer(char *&buff, loint size);
  //! Printing the content of the first buffer
  void print_first_buffer();
  //! Printing the content of the second buffer 
  void print_second_buffer();
  //! @return the merged result of the two internal buffers.
  char *get_string(loint &length) {
    char *str = NULL; get_string(length, str);return str;
  }    
  //! Sets the input param as the merged result of the two internal buffers.
  void get_string(loint &length, char *&str);
  //! Returns an allocated object of the class.
  static buffer_string *init(loint size) {return new buffer_string[size]();} 
  /**
     @brief De-allocates the list of objects given as input, defined
     @param <tmp> The list of objects to de-allocate.
     @param <size_index> The number of objects provided to de-allocate.
   **/
  static void close(buffer_string *&tmp, loint &size_index);
  //! De-allocates the object given as input.
  static void close(buffer_string &tmp) {tmp.free_mem(0);}
  //! De-allocates the memory for this object.
  void free_mem(uint block_id);
  //! De-allocates the memory for this object.
  void free_mem();
  //! The constructor.
 buffer_string() :
  first_buffer(NULL), size_first_buffer(0), second_buffer(NULL), size_second_buffer(0) {};
  //! The test function for this class.
  static void assert_class(const bool print_info);
};

/**
   @brief Holds a block of data read from a text file.
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
   @date 02.11.2011 by oekseth (initial)
**/
typedef buffer_string buffer_string_t;
#endif
