#ifndef index_h
#define index_h
/**
   @file
   @brief Maps a protein-id to a set of other proteins (i.e. identifies the protein pairs).
   @ingroup parsing_ops
   @date 04.01.2012 by oekseth (cleanup).
   @author Ole Kristian Ekseth (oekseth).
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
#include "types.h"
/**
   @file
   @class index
   @brief Maps a protein-id to a set of other proteins (i.e. identifies the protein pairs).
   @ingroup parsing_ops
   @remark Holds the index used in retrieving data for a spesific protein after read from a binary file.
   -- Each index correponds to a specific protein id.
   @author Ole Kristian Ekseth (oekseth).
   @date 18.03.2011 by oekseth (initial).
   @date 27.08.2011 by oekseth (asserts).
   @date 04.01.2012 by oekseth (cleanup).
**/ 
class index { 
 private:
  uint start; // the start pos in index of type rel_t
  uint length; // the total number of elements from the start of type rel_t
 public:
  //! If outer relations are added, decreases the count of them.
  void decrease_length(){if(length > 0) length--;}
  //! @return true if it has data set
  bool has_data() {if(length > 0) return true; return false;}
  /**
     @brief sets an object of type index;
        @param <block> The object to set values into.
     @remarks Transfering data, sets the last param with it. Performs validation
   **/
  void get_data(index &block) {
    block = index(start, length);
  }
  //! Sets the data for this object:
  void set_data(index obj) {
    start = obj.get_start_pos();
    length = obj.get_length();
  }

  //! @return true if they are equal
  bool is_equal(index obj) {
    if((obj.get_start_pos()==start) && (obj.get_length() == length) ) return true;
    else return false;
  }
  //! @return true if 'this' is greater or equal with respect to the number of objects.
  bool is_greater_or_equal(index obj) {
    if((obj.get_start_pos()==start) && (obj.get_length() <= length) ) return true;
    else return false;
  }

  /**
     @param <obj> The object to test against.
     @param <offset_start_pos> The offset to subtract in the verification.
     @return true if they are equal
  **/
  bool is_equal(index obj, loint offset_start_pos) {
    if(((obj.get_start_pos()+offset_start_pos)==start) && (obj.get_length() == length) ) return true;
    else return false;
  }

  /**
     @param <obj> The object to test against.
     @param <offset_start_pos> The offset to subtract in the verification.
     @return true if they are equal
  **/
  bool is_greater_or_equal(index obj, loint offset_start_pos) {
    if(((obj.get_start_pos()+offset_start_pos)==start) && (obj.get_length() <= length) ) return true;
    else return false;
  }
  /**
     @brief Sets the starting position in the calling object.
     @param <start_pos> The index in the calling object to be set.
  **/
  void set_start_pos(uint start_pos) {start = start_pos;}
  //! Returns the index of where the outer relations for this object starts in the calling object.
  uint get_start_pos() {return start;}
  //! @return the end pos, used in for-loop-iterations:
  uint get_end_pos() {return start + length;}
  //! @return the end point, used in for-loop-iterations:
  uint get_end_point() {return start + length;}
  //! @return the number of outer relations this object is connected with..
  uint get_length() { return length; }
  //! @return the absolute postion for the calling object. Used when merging index structures
  uint get_absolute_start_pos() {return start + length;} // @return the end-point: to get the actual index, substract '-1' from this.
  //! Defines the number of outer relations for this object.
  void set_length(uint l) { length = l;}
  //! Increases the number of outer relations by one.
  void increment_length() { length ++; }
  //! Increases the number of outer relations by the param given.
  void increment_length(uint l) { length += l; }
  //! Prints information of the object
  void print(uint i);
  //! Prints information of the object using the name given as param.
  void print(char *name);
  //! Prints information of the object.
  void print(uint i, char **name);

  
  // TODO: Methods below are rpimary replaces by laternatives, so at some point consider removing them.

  //! Prints information of the list of objects given as input.
  static void print_list(index *list, uint index_size) {
    for(uint i = 0; i < index_size; i++) list[i].print(i);
  }
  //! Init a list of index
  static index *init_list(uint index_size) {
    index *list = (class index*)malloc(sizeof(index)*index_size);
    if(!list) {
      fprintf(stderr, "!!\t Was not able to allucate %u=%.3E elements of class \"index\", each with a size of %u B.\n"
	      "-\t If this error is seen, please contact the developer at [oekseth@gmail.com].\n"
	      "Error at [%s]:%s:%d\n", index_size, (float)index_size, (uint)sizeof(index), __FUNCTION__, __FILE__, __LINE__);
      assert(list);
    }
    for(uint i = 0; i<index_size; i++) list[i] = index();
    return list;
  }
  //! Deallocates the memory for the given index
  static void delete_list(index *&index) {
    free(index); index = NULL;
  }

  //! Merges two indexes
  static void merge_lists(index *&destination, uint destination_size, index *source, uint source_size, uint offset) {
    if(destination_size >= source_size) {
      for(uint i = 0; i< source_size;i++) {
	if(source[i].get_length() > 0) {
	  destination[i].set_start_pos(offset + source[i].get_start_pos());
	  destination[i].set_length(source[i].get_length());
	}
      }
    } else fprintf(stderr, "!!\tAn error occured in merge(..) foudn in class 'index'\n");
  }
  
  /**
     @brief Enlarges the current buffer. Optional parameters may be set.
     @param <buffer> The buffer to be enlarged.
     @param <new_size> The new size of the buffer.
     @param <start_pos> The start position to set the old buffer (an optional param, normally set to '0').
     @param <copy_length> The length of the buffer to copy into the structure.
  */
  static void enlarge(index *&buffer, uint new_size, uint start_pos, uint copy_length) {
    index *buff_temp = (index*)malloc(sizeof(index)*new_size);
    if(!buff_temp) {
      fprintf(stderr, "!!\t Was not able to allucate %u=%.3E elements of class \"index\", each with a size of %u B.\n"
	      "-\t If this error is seen, please contact the developer at [oekseth@gmail.com].\n"
	      "Error at [%s]:%s:%d\n", new_size, (float)new_size, (uint)sizeof(index), __FUNCTION__, __FILE__, __LINE__);
      assert(buff_temp);
    }
    uint pos = start_pos;     if(copy_length > 1) pos+=copy_length-1;
    for(uint i = pos; i < new_size; i++) buff_temp[i] = index();     // Initiates the buffer:
    // Inserts the buffer in the front of the file (i.e. in
    // order for the indexes to maintain their correct position pointers:
    if(buff_temp != NULL) {  // Copys the old data into the new structure:
      memcpy(buff_temp +start_pos, buffer, sizeof(index)*copy_length);
      free(buffer);
      buffer = buff_temp;
    } else fprintf(stderr, "!!\tNot enough memory for storying %u elements of type index_t\n", new_size);
  }

  /**
     @return The total number of elements the lists represents in the file fo the span given.
  */
  static uint get_number_of_elements(index *list, uint start, uint end, uint size) {
    if(size > 0 && end <= size) { // has data set
      uint sum_outer = 0;
      for(uint protein_in = start; protein_in < end; protein_in++)
	sum_outer += list[protein_in].get_length();
      return sum_outer;
    }
    return 0;
  }

  /**
     @return true is data is set for the params.
  **/
  static bool has_data(index *list, uint index_in, uint size) {
    if(list != NULL) // has data stored
      if(index_in < size) {
	if(list[index_in].get_length() > 0) return true;
	else return false;
      }
    return false;
  }

  // End of methods to be considered removed.
  //! An alternative to a standard constructor.
  void init() {start=0, length=0;}

  //! The constructor.
 index() : start(0), length(0) {
};
  //! The constructor.
 index(uint _s, uint _l): start((uint)_s), length((uint)_l) {};
  /**
     @brief The main test function for this class.
     @remarks This method includes:
     - Formalized tests for setting of starting position and end position.
     - Valgrind used verifying memory usage.
     - Examples of how this class may be used.
     @author Ole Kristian Ekseth (oekseth)
     @date 04.01.2012 by oekseth.
  **/
  static void assert_class(const bool print_info);

}; 

/**
   @brief Maps a protein-id to a set of other proteins (i.e. identifies the protein pairs).
   @ingroup parsing_ops
 */
typedef class index index_t;

#endif
