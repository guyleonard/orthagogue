#ifndef taxa_buffer_list_settings_h
#define taxa_buffer_list_settings_h
/** 
   @file
   @brief Holds the size information about the data containers used for mpi sending- and receiving.
   @author Ole Kristian Ekseth (oekseth)
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

#include "taxon_data.h"
#include "../configure.h"
#include "tbb_libs.h"
#include "log_builder.h"

#ifdef USE_MPI
/**
   @struct taxa_buffer_list_settings
   @brief Holds the size information about the data containers used for mpi sending- and receiving.
   @remarks 
   # Seperated for structure taxa_list_of_lists_settings in order to emobdy the seperation between what's sent seperatly, and what is sendt wrapped as an mpi container.
   # Only acitvated when amcro-variable USE_MPI is set.
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
   @date 10.07.2012 by oekseth (initial)
**/
struct taxa_buffer_list_settings {
  uint taxon_size; // number of taxa.
  uint taxonNamesArr_size; // Number of chars used holding the taxon names.
  uint total_cnt_proteins; // The total number of proteins in the collection.
  uint total_cnt_chars_used_representing_proteins; // The number of chars used holding the protein names.
  void print_variables(int myrank) {
    //    if(false && myrank) ;
    printf("[%d]\t taxon_size(%u), taxonNamesArr_size(%u), total_cnt_proteins(%u), total_cnt_chars_used_representing_proteins(%u) at line %d in file %s\n", myrank, taxon_size, taxonNamesArr_size, total_cnt_proteins, total_cnt_chars_used_representing_proteins, __LINE__, __FILE__);
  }
  taxa_buffer_list_settings() : taxon_size(0), taxonNamesArr_size(0), total_cnt_proteins(0), total_cnt_chars_used_representing_proteins(0) {};  
};

/**
   @struct taxa_buffer_list_settings
   @brief Holds the size information about the data containers used for mpi sending- and receiving.
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
   @date 10.07.2012 by oekseth (initial)
**/
typedef struct taxa_buffer_list_settings taxa_buffer_list_settings_t;

#endif

#endif
