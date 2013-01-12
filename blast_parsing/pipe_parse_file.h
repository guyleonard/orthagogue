#ifndef pipe_parse_file_h
#define pipe_parse_file_h
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
#include "parse_read_blocks.h"
#include "string_section.h"
#include "read_file.h"
#include "list_file_parse.h"
#include "log_builder.h"
#include "buffer_string_list.h"
#include "../configure.h"
/**
   @file
   @class pipe_parse_file
   @brief Reads the blast file sending it down the pipe.
   @ingroup pipe_parsing
   @remark The two reading procedures are:
   - Reads a block of lines for each call to the operator or
   - reads either to be called by the 'read-file-blocks' method or via the tbbb operator using pipeline architecture.
   - Validates reading of data by generating an empty test-set with a predefined set of newline-chars.
   @return A "string_section_t" block of text data, starting on a line start, and ending on a line end
   with a sze of approx the value set at the variable "disk_buffer_size"
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
**/
class pipe_parse_file: public tbb::filter {
 private:
  char seperator; // Used for validation of the reading
  uint disk_buffer_size;
  log_builder_t *log;
  uint block_cnt; // counts the blocks used
  parse_read_blocks *parseBlocks; // holds the blocks to read
  read_file_t *fileRead;
 public:
  //! Set to true in the constructor, and 'false' in 'rewind_file(..)'
  bool FIRST_RUN; 
 private:
  char *FILE_INPUT_NAME;
  loint chars_in_file_read_first_run;
  int myrank;
 public:
  /**
     @return An object consisting of a char-block with extra meta-data.
     @remark Using method found in object of type read_file.
  **/
  struct string_section *read_file_blocks();
  //  int curr_taxa_in;   int curr_taxa_out;

  //! Returns the section of the data using the parse-blocks as a basis:
  string_section *get_section();
  //! Sets the file to the start_pos of the file
  void rewind_file(parse_read_blocks *_parseBlocks);
  //! To be called after the operation is finalized
  void close_file(int n_threads);
  //! The method of parallisation.
  /*overrride*/void* operator()(void* item);
  //! Constructor
  pipe_parse_file(char *f_name, log_builder_t *_log);
  //! Constructor
  pipe_parse_file(char seperator, uint _disk_buffer_size, char *f_name, log_builder_t *_log);
#ifdef USE_MPI
  //! Constructor
  pipe_parse_file(char seperator, uint _disk_buffer_size, MPI_Comm comm, int myrank, char *f_name, log_builder_t *_log, loint reading_file_start_position, loint reading_file_end_position, loint reading_file_length);
  //  pipe_parse_file(MPI_File &fh, int myrank, char *f_name, log_builder_t *_log, loint reading_file_start_position, loint reading_file_end_position, loint reading_file_length);
#endif
#ifdef assert_code
  //! The test function for the private parts of this class.
  void assert_private_parts();
#endif 
  //! The main test function for this class  
  static void assert_class(const bool print_info);

  
};

#endif
