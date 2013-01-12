/**
   @file
   @brief Holds a pair of proteins.
   @ingroup parsing_ops
**/
#ifndef protein_relation_h
#define protein_relation_h
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
   @class protein_relation
   @brief Holds a pair of proteins.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 08.01.2010 by oekseth (initial)
   @date 26.08.2011 by oekseth (asserts)
**/
class protein_relation {
  bool protein_relation_exsists;
 public:
  //! The inner taxon
  int taxon_in;
  //! The inner protein
  mem_loc protein_in;
  //! The outer taxon
  int taxon_out;
  //! The outer protein.
  mem_loc protein_out;  
  //! If set, pair will de discarded (i.e. not included) in the process of building the compressed list_file_parse object.
  void set_as_no_exsists(){protein_relation_exsists=false;}
  //! @returns true if data is set (implying data were found) for the pair this object represents.
  bool exsists(){return protein_relation_exsists;}
  //! Returns the inner protein set.
  mem_loc get_protein_in(){return protein_in;}
  //! Returns the outer protein set.
  mem_loc get_protein_out(){return protein_out;}
  //! Prints the fields set
  void print();
  //! Returns true if data is set.
  bool is_set();
  //! Copies the argument into this.
  void copy(protein_relation prot);
  //! @ Returns an a list of objects.
  static protein_relation* init(uint size) {return (protein_relation*)malloc(sizeof(protein_relation)*size);}
  //! Releases a list of object.
  static void free_list(protein_relation *&obj) {free(obj); obj= NULL;}
  //! The constructor.
  protein_relation();
  //! The constructor.
  protein_relation(int _taxon_in, mem_loc _protein_in, int _taxon_out, mem_loc _protein_out);

  //! The main test function for this class  
  static void assert_class(const bool print_info);

};

#endif
