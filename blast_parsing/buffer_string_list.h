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
#ifndef buffer_string_list_h
#define buffer_string_list_h
#include "buffer_string.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "string_section.h"
#include "../configure.h"
/**
   @file
   @class buffer_string_list
   @brief Inserts the data block at the correct location, and retrieves iaw procedure.
   @ingroup parsing_ops
   @remark May be regarded as a wrapper for class 'buffer_string'. 
   - Hopefully all permutations of the inputs are tested.
   @author Ole Kristian Ekseth (oekseth)
   @date 02.11.2011 by oekseth (initial)
**/
class buffer_string_list {
 private:
  loint current_index;
  loint size_index;
  buffer_string *list;
  //! Merges two strings:
  void merge(char *&str, loint &length, char *&temp, loint temp_length);
 public:
  //! @return true if more blocks available to be read
  bool has_data(){
    if(current_index >= size_index) return false;
    else if(!list) return false;
    else if(list[current_index].has_first_block()) return true;
    else if(list[current_index].has_main_block()) return true;
    else return false;
  }
    /**
     @brief Generate (writes) a log file describing the details of the memory signature for strings stored in memory during parsing
     @remarks Useful for analyzing- and optimizing the memory allocation procedures.
  **/
  void log_generate_memory_allocation_overview(uint n_threads);
  //! Prints the properties of the item last inserted.
  void print_item_last_inserted();
  //! Prints the properties of the items in the list
  void print_items();
  //! @return the length (size) of the index.
  loint get_index_length(){return size_index;}
  //! Returns the string given, using the blast file stored in memory.
  void get_string(loint &length, char *&str, loint &found_at_index_pos);
  //! Returns false if no mare data is of the file in memory
  bool get_section(string_section *&section, loint &found_at_index_pos);
  //! Inserts the overflow buffer at the correct location.
  void copy_into_first_buffer(uint block_cnt, char *buff, loint size);
  //! Inserts the main buffer at the correct location.
  void set_second_buffer(uint block_cnt, char *buff, loint size);
  /**
     @brief Includes the data given into the buffer available for all the threads.
     @param <section> The object containing the chars from the blast file to store.
     @param <block_length> The remaining length of data to be merged with the block following this data in the input blast file.
     @return false if data is not inserted into memory.
     @remarks If the index block provided is above the length of the list reserved, data is not stored in this object. This implies that alternative method of getting the data is required. In the gogueOrto project class thread_reading provides this option. For examples of usage, see class pipe_parse_parse.
     @author Ole Kristian Ekseth (oekseth)
     @date 05.01.2012 by oekseth (documentation and debuging).
   **/
  bool update(string_section_t *&section, loint block_length);

  //! Returns the length of the item inserted:
  loint get_length_of_inserted_data(string_section_t *section) {
    if(!section || !list) return 0;
    else {
      loint sum = 0;
      if((loint)(section->block_cnt+1) < size_index) sum = list[section->block_cnt+1].get_size_first_buffer();
      if((loint)(section->block_cnt) < size_index) sum += list[section->block_cnt].get_size_second_buffer();
      return sum;
    }
    
  }
  //! Return the total number of chars represented in the given blocks:
  loint getTotalLengthOfData() {
    if(list) {
      loint sum = 0;
      for(loint i = 0; i < size_index; i++) {
	sum += list[i].get_size();
      }
      return sum;
    } else return 0;
  }
  //! Returns the file size.
  static loint get_file_size(char *file);
  //! Returns an approximation of how mmuch memory may be reserved for this purpose.
  static loint get_max_size_in_buff(char *f_name);
  //! Returns an estimated number of the block_count that may be reserved.
  static loint get_estimate_of_memory_blocks_cnt(uint disk_buffer_size, char *f_name);
  //! De-allocates memory for the object given as input.
  static void close(buffer_string_list &tmp);
  //! De-allocates memory for the object given as input.
  static void close(buffer_string_list *&tmp);
  //! De-allocates memory for this class.
  void free_mem();  
  //! The constructor.
  buffer_string_list(loint size);
  //! The test function for this class.
  static void assert_class(const bool print_info);
};

/**
   @brief Inserts the data block at the correct location, and retrieves iaw procedure.
   @ingroup parsing_ops 
   @author Ole Kristian Ekseth (oekseth)
   @date 02.11.2011 by oekseth (initial)
**/
typedef class buffer_string_list buffer_string_list_t;

#endif
