#ifndef bp_container_h
#define bp_container_h 
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
#include "../configure.h"
#include "log_builder.h"
#include "list_norm.h"
/**
   @file
   @struct bp_container
   @brief The container holding result of blast parsing: wraps the output data from class blast_parsing.
   @ingroup transfer_container
   @author Ole Kristian Ekseth (oekseth)
   @date 28.12.2010 by Ole Kristian Ekseth (init)
**/
struct bp_container {
  //! Container storing the normative values
  list_norm_t *arrNorm; 
  //! Container storing info about each taxonfound.
  taxa_t *listTaxa; 
  //! The number of taxa in the collection.
  int taxon_length; 
  //! Holds the blast file in a modified form.
  list_file_parse_t *listParseData; 
  //! The maximum similairty score found in the blast file used
  float max_input_value; 

  //! the name (the whole path) of the input file
  char *FILE_INPUT_NAME; 
  //! Locates this objects varibles into the parameters defined.
  void get_variables(list_norm_t *&arr, taxa_t *&tax, int &t_length, list_file_parse_t *&parse, float &m_value) {
    assert(!arr);
    assert(!tax);
    assert(!parse);
    assert(!t_length);
    assert(!m_value);
    arr = arrNorm, tax = listTaxa, parse = listParseData, t_length = taxon_length, m_value = max_input_value;
    // TODO: Should the internal variables now be reset?
  }

  //! Returns true if data is set
  bool has_data() {
    if(listTaxa) {
      if(listParseData) {return true;}
    }    
    return false;
  }


  /**
     @brief Relases the memory located in tis object.
     @remarks To be with care, as other processes depends on them
  **/
  void free_memory() {
    if(listParseData) printf("data allocated for in blast_filtering.cxxlistParseData\n");
    taxa::close(listTaxa, taxon_length); 
    if(listParseData)  listParseData->free_memory(true);  
    list_norm::close(arrNorm); 
  }
  //! The constructor.
bp_container(list_norm_t *_arrNorm, taxa_t *_listTaxa, int _taxon_length, list_file_parse_t *_listParseData, float _max_input_value) :
  arrNorm(_arrNorm), listTaxa(_listTaxa), taxon_length(_taxon_length), listParseData(_listParseData), max_input_value(_max_input_value) {};
bp_container() : 
  arrNorm(NULL), listTaxa(NULL), taxon_length(0), listParseData(NULL), max_input_value(0), FILE_INPUT_NAME(NULL) {};
};

/**
   @brief The container holding result of blast parsing: wraps the output data from class blast_parsing.
   @author Ole Kristian Ekseth (oekseth)
   @date 28.12.2010 by Ole Kristian Ekseth (init)
**/
typedef struct bp_container bp_container_t;
#endif
