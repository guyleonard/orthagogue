#ifndef taxon_pair_jj_h
#define taxon_pair_jj_h
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
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h> // some systems require it
#include <sys/stat.h>
#include <sys/termios.h> // for winsize
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/times.h>
#include "libs.h"
#include "taxa.h"
#include "list_file_parse.h"
#include "log_builder.h"
/**
   @file
   @class taxon_pair
   @brief Defines the next protein to work on.
   @ingroup filtering_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of this class as a libary)
   @date 31.12.2011 (cleanup)
**/
class taxon_pair {
 public:
  //! The taxon id of interest
  uint taxon; 
  //! The first protein of this collection.
  uint protein_start;
  //! The last protein of this collection.
  uint protein_end; 
 private:
  list_file_parse_t *listParseData; 
  list_file_struct_t *listStructData; 
  taxa_t *listTaxa;
  //! Generates a file for a set of this class:
  static void write_class_list(uint index, taxon_pair *pair, uint taxon_start, uint taxon_end,
				     uint taxon_length, uint n_threads, 
			       const bool only_inpa, list_file_parse_t *listParseData, list_file_struct_t *listStructData, taxa_t *listTaxa);

 public:
  //! @return the length of the proteins in the collection.
  uint get_length() {return (protein_end - protein_start);}
  //! @return the taxon for this object.
  uint get_taxon(){return taxon;}

  //! @return the first protein for this object.
  uint get_protein_start(){return protein_start;}

  //! @return the last protein for this object.
  uint get_protein_end(){return protein_end;}
  //! @return the list_file_struct_t for this object.
  list_file_struct_t *get_listStructData(){return listStructData;}
  //! @return the list_file_parse_t for this object.
  list_file_parse_t *get_listParseData(){return listParseData;}

  //! @return the list_file_parse_t for this object.
  taxa_t *get_listTaxa(){return listTaxa;}
  //! Sets the values for this object.
  void set_values(taxon_pair &obj);
  //! @returns the variables:
  void getVariables(uint &_taxon, uint &_protein_start, uint &_protein_end); 
//! Prints the data
  void print();   
  //! @returns true if is has more data
  bool has_data(); 
  //! @return an allocated region: A standard interface to ensure consistency
  static taxon_pair *initList(uint size);

  /**
     @return a taxon_pair list of data
  **/
  static taxon_pair *get_taxon_pairs(taxon_pair *listPair, uint &list_pair_pos);

  //! Closes the given taxon_pair object given.
  static void close(taxon_pair *&obj) { if(obj) {free(obj); obj = NULL;}}
  //! A standard interface to ensure consistency:
  static void free_taxonPairList(taxon_pair *list) {free(list);}  
  /**
     @brief Builds the set of blocks to be used during the parsing
     @param <taxon_start> The id of the first taxa to build the list for.
     @param <taxon_end> The id of the first taxa (and the following onwards) to not include in the list build.
     @param <taxon_length> The total number of taxa in the collection used as input.
     @param <n_threads>    The number of cpu's (threads) used for this operation.
     @param <only_inpa> If set, only calculates the size of its inparalogs.
     @param <listParseData> Holds data from the parsing operation.
     @param <listStructData> Holds data from the filtering operation.
     @param <listTaxa> Holds data about the taxa involved in this operation.
     @remarks 
     # Used by methods originating in class blast_filtering.
     # Depends on scheduling functions found in class list_file_parse.
     # Divid blocks using approximations of the best sizes of data.
     @author Ole Kristian Ekseth (oekseth).
  **/
  static taxon_pair *init_taxon_pair(uint taxon_start, uint taxon_end,
				     uint taxon_length, uint n_threads, 
				     const bool only_inpa,
				     list_file_parse_t *listParseData, list_file_struct_t *&listStructData, taxa_t *listTaxa
				     );
  //! The constructor
  taxon_pair(uint _taxon, uint _protein_start, uint _protein_end,
	     list_file_parse_t *&_listParseData, list_file_struct_t *&_listStructData, taxa_t *_listTaxa);
  //! The constructor
  taxon_pair();

  //! The main test function for this class
  static void assert_class(const bool print_info);
};

/**
   @brief Defines the next protein to work on.
   @ingroup filtering_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of this class as a libary)
   @date 31.12.2011 (cleanup)
**/
typedef class taxon_pair taxon_pair_t;

#endif
