#ifndef basket_parse_t
#define basket_parse_t
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
#include "p_rel.h"
#include "rel.h"  
#include <float.h>
/**
   @file 
   @class basket_parse
   @brief A wrapper class assuring that only data in memory are read.
   @ingroup filter_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 14.10.2011 by oekseth (initial)
   @date 14.10.2011 by oekseth (asserts)
**/
class basket_parse {
 private:
  p_rel_t *buffer;
  loint buffer_size;
  bool mem_allocated;
 public:
  //! @return the number of elements
  loint get_buffer_size(){return buffer_size;}

  //! @return true if buffer is set.
  bool has_data() {return (buffer != NULL);}
  //! Tests if index is inside the range of the list.
  bool has_data(uint index_offset, uint index) {
    return in_mem(index_offset, index);
  }

  /** 
      @brief Appends the total score sum for the data in this container.
      @remarks 
      # used for validation of correctness.
      # Increments the values (instead of the standard 0-initiating, we use the existing values in them, and append values on the alrady set).
  **/
  void get_properties_and_increment(uint &cnt, float &sum, uint &zeros) {
    if(buffer) {
      cnt += buffer_size;      
      for(uint i = 0; i < (uint)buffer_size; i++) {
	sum += buffer[i].get_distance();
	zeros += buffer[i].get_cnt_max();
      }
    }
  }
  //! Tests if index is inside the range of the list.
  bool in_mem(uint index_offset, uint index) {
    if(buffer) {
      const uint absolute_index = index_offset + index;
      if(absolute_index < buffer_size) {
	return true;
      } else {
	return false;
      }
    } else return false;
  }
  //! Returns the variable 'index out' variable from the struct at the index pos.
  uint get_index_out(uint index_offset, uint index) {
    const uint absolute_index = index_offset + index;
    if(in_mem(index_offset, index)) 	return buffer[absolute_index].get_ind_out();
    else {
      fprintf(stderr, "!!\tData not found for index(%u) with offset(%u) and buffer(%d). Error is found at line %d in file %s located in vunction %s. If this message is seen, please contact oekseth@gmail.com\n", index, index_offset, (buffer!=NULL), __LINE__, __FILE__, __FUNCTION__);
      assert(false);
      return UINT_MAX;      
    }
  }

  //! Returns the given overlap for either outer- or inner of the struct:
  overlap_t get_overlap(uint index_offset, uint index, bool is_inner) {
    const uint absolute_index = index_offset + index;
    assert(in_mem(index_offset, index)); // Do only validation in assert-mode
    if(is_inner) return buffer[absolute_index].get_overlap_in();
    else return buffer[absolute_index].get_overlap_out();
  }
  //! Returns the given overlap for the inner of the struct:
  overlap_t get_overlap_in(uint index_offset,uint index){return get_overlap(index_offset, index, true);}
  //! Returns the given overlap for the outer of the struct:
  overlap_t get_overlap_out(uint index_offset,uint index){return get_overlap(index_offset,index, false);}
  //! Returns the given sim-score
  float getSimScore(uint index_offset, uint index, float max_input_value) {
    const uint absolute_index = index_offset + index;
    if(in_mem(index_offset, index)) 	return buffer[absolute_index].getSimScore(max_input_value);
    else return FLT_MAX;
  }
  //! Frees the mem for 'this'.
  void free_mem() {
    if(mem_allocated) {
      p_rel::free_file_block(buffer, buffer_size);
    }
  }
  //! Frees the mem for the object.
  static void free_mem(basket_parse *&temp) {
    if(temp) {
      temp->free_mem(); delete temp; temp = NULL;
    }
  }

  //  static class basket_parse *init(uint buffer_size_, bool mem_allocated) {
  /**
     @return an object of type basket_parse using the input params as input to the constructor.
   */
  static class basket_parse *init(p_rel_t *buffer_, uint buffer_size_, bool mem_allocated) {
    return new basket_parse(buffer_, buffer_size_, mem_allocated);
  }
  //! Constructs the class setting variables to empty
 basket_parse() : buffer(NULL), buffer_size(0), mem_allocated(false) {};
  //! Constructs the class given the input params.
 basket_parse(p_rel *buffer_,
	      uint buffer_size_,
	      bool _mem_allocated) :
  buffer(buffer_),
    buffer_size(buffer_size_), 
    mem_allocated(_mem_allocated)
    {};

};
#endif
