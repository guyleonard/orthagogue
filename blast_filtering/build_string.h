#ifndef build_string_h
#define build_string_h
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
#include "../configure.h"
#include "types.h"
#include "tbb_libs.h"
#include "log_builder.h"
/**
   @file
   @class build_string
   @brief Produces a row of chars for the protein given.
   @ingroup filter_container
   @author Ole Kristian Ekseth (oekseth
   @date 21.12.2010 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 31.12.2011 by oekseth (cleanup)
**/
class build_string {
 private:
  long int last_index;
  uint size_blast_all;
  //! The position of the strings logical start.
  char* logical_start;
  //! The position of the strings logical end.
  char* logical_end;  
  //! Prints the buffer stored in this object.
  void print();
 public:
  //! Writes the string to the file, and clears the memory
  void write_data_to_file(FILE  *out_file, const bool PIPE);
  //! Sets the header for the string to be written to the output file.
  void setHeader(uint world_index_in);
  //! Called when no data shall be printed to the stream
  void avoidPrinting();
  //! If set to true, the data in the poiners are written to a file.
  void writeRow(bool _has_data);
  //! @return true if data shall be printed
  bool hasData();
  //! @return pointer to beginning of the sequence the string represents.
  char* begin();
  //! Appends the chars to the buffer
  void append(char* first, char* last);
  /**
     @brief adds the content of the input to the string located in this object.
     @param <start> Pointer to the first char in the sequence of chars.
     @param <end> Pointer to the last char in the sequence of chars.
  */
  void copy(char *start, char *end);

  //! adds the content of the input to the string located in this object.
  void copy_struct(build_string row);
  
  //! Appends the line-end characters specified for the 'blast' file format
  void end_blast_line();
  //! Returns the size of the index
  long int getIndex();

  //! Appends the line end ('\0' character)
  void finalize();
  //! Length of sequence
  uint size() {return(logical_end-logical_start);}

  /**
     @brief Increases the size of the buffer if it does not fit.
     @param <length_argument> the length of the input.
     @return the size of the argument to be copied.
     @remarks To be used at insertion of a new part of the char array
  */ 
  mem_loc resize_if_to_large(mem_loc length_argument);

  //! Allocates a list of 'this' class, but do not initialize it.
  static build_string *allocate_class(uint size);

  //! Frees the memory allocated for this object
  void free_mem();
  //! The constructor.
  build_string(uint _size_blast_all);
#ifdef assert_code
  //! Asserts of private pars of this class.
  void assert_private_parts();
  //! The main test function for this class  
  static void assert_class(const bool print_info);
#endif
};

/**
   @brief Produces a row of chars for the protein given.
   @ingroup filter_container
   @author Ole Kristian Ekseth (oekseth
   @date 16.09.2011 by oekseth (asserts)
   @date 31.12.2011 by oekseth (cleanup)
**/
typedef build_string build_string_t;
#endif
