#ifndef string_section_h
#define string_section_h
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
#include "blast_extractors.h"
#include "types.h"
#include "tsettings_input.h"
#include "taxon_list_t.h"
#include "parse.h"
/**
   @file
   @struct string_section
   @brief A transport container for a block of chars, holding metadata about it.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 15.09.2011 by oekseth (asserts)
**/
struct string_section {
  //! The start position of the buffer.
  char *buffer_main_start;
  //! The end pos of the buffer.
  char *buffer_main_logical_end;
  //! The previous index data were inserted at.
  long int prev_index;
  //! The relative postion among the blocks (in order to ensure order_preserving)
  int block_cnt;
  //! Set to true if reading is on the end of the file.
  bool end_of_file;
  //! @return the length of the string
  mem_loc getStringLength();
  /**
     @brief An assert function returning true if input equals
     @return true if input equals the content of this object's 'buffer_main_start' variable.
  **/
  bool is_equal(char *inp, long int inp_length);
  //! @return the number of chars containing data.
  long int get_data_length();
  //! Prints the data block given the argument.
  void print_data_block(long int end_pos);
  //! Prints the data block
  void print_data_block();
  //! Prints the segment given the inputs
  static void print_segment(char *m_start, char *end);
  //! @return true if object contains data by comparing the memory addresses.
  bool is_set();
private:
  /**
     @return the file position in the file where the next new line in the blast file starts.
     @remarks Return 0 if a shift from one type of left-most protein to another is not found.
  **/
  loint get_last_line_pos(char *buffer_start, char *last_block_line);
public:
  /**
     @brief Tests if block has a shift from one type of left-most protein to another. 
     @return A value greater than 0 if a position is found.
  **/  
  loint more_than_one_leftmost_protein(taxon_list_t *listProteins, tsettings_input_t blast_settings);
  //! Deallocates the memory
  void finalize();
  //! @return An object of this given the input params.
  static string_section *init(char *buffer_main_start,
			      char *buffer_main_logical_end,
			      long int chars_sendt_to_parsing,
			      int block_cnt, bool end_of_file
			      );
  //! Deallocates the memory given the param.
  static void close(string_section *&obj);
  //! The constructor.
  string_section();
  //! The constructor.
  string_section(char *_buffer_main_start,
		 char *_buffer_main_logical_end,
		 long int _prev_index, int _block_cnt, bool _end_of_file);
  //! The main test function for this class  
  static void assert_class(const bool print_info);  
};
/**
   @brief A transport container for a block of chars, holding metadata about it.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 15.09.2011 by oekseth (asserts)
**/
typedef struct string_section string_section_t;

#endif
