#ifndef UINT_3x_OBJ_H
#define UINT_3x_OBJ_H
/**
   @file
   @brief Build a set of The id's for easy- and consisten index retrieval, used in eg, sending and receiving MPI data
   @remarks Used for MPI transfer.
   @ingroup blastfile_container
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

#include "log_builder.h"

/**
   @struct uint_3x_obj
   @brief Build a set of The id's for easy- and consisten index retrieval, used in eg, sending and receiving MPI data
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
struct uint_3x_obj {
  //! The id (or when used as object the length):
  uint taxon_id;
  //! The number of type T elements in the list to send
  uint total_length_of_rel;
  //! The number of type index_t elements in the list to send
  uint total_length_of_index;
  //! Prints the content:
  void print_object() {
    printf("-\ttaxon_id[%u] has relations(%u) and indexes(%u) at line %d in file %s\n", taxon_id, total_length_of_rel, total_length_of_index, __LINE__, __FILE__);
  }
  /**
     @return true if data set
     @remarks: If one of them is not set, it does not correspond to our expectation.
  **/
  bool is_set() {
    if(!total_length_of_index || !total_length_of_rel) return false;
    return true;
  }
  /**
     @brief The constructor for the sender.
     @remarks: If one of them is not set, it does not correspond to our expectation.
  **/
  uint_3x_obj(uint id, uint cnt_rel, uint cnt_index) :
    taxon_id(id), total_length_of_rel(cnt_rel), total_length_of_index(cnt_index)
  {
    if(!total_length_of_index || !total_length_of_rel) {
      assert(false);
    }
};
  //! The constructor for the receiver.
  uint_3x_obj() : taxon_id(0), total_length_of_rel(0), total_length_of_index(0) {};
}; // uint_3x_obj;

/**
   @brief Build a set of The id's for easy- and consisten index retrieval, used in eg, sending and receiving MPI data
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
typedef struct uint_3x_obj uint_3x_obj_t;
#endif
