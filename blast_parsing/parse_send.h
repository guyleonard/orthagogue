#ifndef parse_send_h
#define parse_send_h
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
#include "list_file_parse.h"
#include "protein_relation.h"

/**
   @file
   @struct parse_send
   @brief Holds the data between two pipes
   @ingroup transfer_container
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 25.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/
struct parse_send {
  //! The number of taxa this objects holds an collection for.
  int taxon_length;
  /**
     @brief Stores the overlap values for each protein pair.
     @remark In order to copy the overlap values and clean the memory used.
  **/
  overlap_t **arrOverlap; 
  //! The maximum found similarity value during this (the callers) operation.
  float max_sim_value; 
  //! The blast data stored.
  list_file_parse_t *parseData; 
  //! Information about the last protein found; may be used during the merging process.
  protein_relation prot;
  //! The global index position the data was found at.
  loint found_at_index_pos;
  //! Prints the content of the structs
  void print_parse_struct();
  //! @return true if data is set.
  bool has_data() {return (parseData!= NULL);}
  //! Prints the data stored in this object, i.e. the protein relations.
  void print();
  //! De-allocates the memory reserved for this object.
  void finalize();
  //! The constructor
  parse_send();
  //! The constructor
  parse_send(list_file_parse_t *_parseData, protein_relation _prot, overlap_t **arr, float _max_sim_value, int _taxon_length, loint found_at_index_pos);
};
/**
   @brief Holds the data between two pipes
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 25.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/
typedef parse_send parse_send_t;

#endif
