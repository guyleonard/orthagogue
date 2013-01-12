#ifndef parse_send_first_h
#define parse_send_first_h
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
#include "taxon_list_t.h"
/**
   @file
   @struct parse_send_first
   @brief Holds the data between class 'pipe_parse_parse' and class
   'pipe_pipe_merge' duirng the first (of two) parsing-operations.
   @ingroup transfer_container
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
**/
struct parse_send_first {
public:
  //! Contains the protein labels found in the blast-section processed by the class 'pipe_parse_parse'.
  taxon_list_t *t_list;
  //! Actual file position of the last line in the blast file.
  long long int last_block_pos;
  //! The block number set by the pipe reading the blocks prosessed in class 'pipe_parse_parse', implying class 'pipe_parse_file', and used to assure correctness in the merging of data in class 'pipe_parse_merge'
  uint block_number;
  //! By default set to false.
  bool data_resides_in_mem;
  //! @return The taxon-list produced in the previous pipe.
  taxon_list_t *getTaxonList() {return t_list;}
  //! Deallocates the memory
  void free_taxon_list_mem();
  //! Deallocates the pointer
  void finalize();
  //! The constructor
  parse_send_first(taxon_list_t *list, long long int _last, uint _block_cnt);
};

/**
   @brief Holds the data between class 'pipe_parse_parse' and class
   'pipe_pipe_merge' duirng the first (of two) parsing-operations.
   @ingroup transfer_container
   @author Ole Kristian Ekseth (oekseth)
   @date 16.09.2011 by oekseth (asserts)
**/
typedef struct parse_send_first  parse_send_first_t;
#endif
