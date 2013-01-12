#ifndef tsettings_input_h
#define tsettings_input_h
/**
   @file
   @brief Holds the properties for the blast file.
   @ingroup parsing_ops
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

#include "../blast_common/types.h"
/**
   @file
   @struct tsettings_input
   @brief Holds the properties for the blast file.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
**/
struct tsettings_input {
  //! The index in the blast file where the protein name is found
  int INDEX_IN_FILE_FOR_PROTEIN_NAME; 
  //! The index in the blast file where the taxon name is found
  int INDEX_IN_FILE_FOR_TAXON_NAME; 
  //! Use the score values in the column 12 of BLAST output instead of the e-values (default, column 11)
  bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE;
  //! Field separator in FASTA (blast) headers (dividing the taxon from the protein label).
  char SEPERATOR;
  //! Path to the BLAST output file in the tabular format
  char *FILE_NAME;
  //! The path to store the folders whom the binary data will be located at
  char *FILE_BINARY_LOCATION; 
  //! An initial guess on how many columns there are in the blast file with regard to the name-label.
uint DEFAULT_NUMBER_OF_COLUMNS_IN_NAME;
  //! If set, prints the data used as basis for the normalization procedure.
  bool DEBUG_NORM; 
  //! If set, prints the data used for normalization of the similarity scores.
  bool PRINT_NORMALIXATION_BASIS; 
  //! The numbers of cpu's to build the paralisation for: only approximate value is neccassary to be set.
  int CPU_TOT; 
  //! If set, uses all realtions in the input as basis for the normalization preocedure. Alternative, only those values filtered out, is to be used.
  bool USE_EVERYREL_AS_ARRNORM_BASIS; 
  //! The constructor
tsettings_input():
  INDEX_IN_FILE_FOR_PROTEIN_NAME(0), INDEX_IN_FILE_FOR_TAXON_NAME(1), USE_LAST_BLAST_CLOMUN_AS_DISTANCE(false), SEPERATOR('|'), FILE_NAME(NULL), DEFAULT_NUMBER_OF_COLUMNS_IN_NAME(5), DEBUG_NORM(false), PRINT_NORMALIXATION_BASIS(false), CPU_TOT(2), USE_EVERYREL_AS_ARRNORM_BASIS(false)
  {};
  //! The constructor
tsettings_input(int ind_prot, int ind_taxon, bool last_col, char sep, char *f_name, uint def_col_cnt):
  INDEX_IN_FILE_FOR_PROTEIN_NAME(ind_prot), INDEX_IN_FILE_FOR_TAXON_NAME(ind_taxon), USE_LAST_BLAST_CLOMUN_AS_DISTANCE(last_col), SEPERATOR(sep), FILE_NAME(f_name), DEFAULT_NUMBER_OF_COLUMNS_IN_NAME(def_col_cnt), DEBUG_NORM(false), PRINT_NORMALIXATION_BASIS(false), CPU_TOT(2), USE_EVERYREL_AS_ARRNORM_BASIS(false)
  {};
  //! The constructor
tsettings_input(int ind_prot, int ind_taxon, bool last_col, char sep, char *f_name, uint def_col_cnt, bool debug_norm, bool print_norm, int cpu, bool u_all):
  INDEX_IN_FILE_FOR_PROTEIN_NAME(ind_prot), INDEX_IN_FILE_FOR_TAXON_NAME(ind_taxon), USE_LAST_BLAST_CLOMUN_AS_DISTANCE(last_col), SEPERATOR(sep), FILE_NAME(f_name), DEFAULT_NUMBER_OF_COLUMNS_IN_NAME(def_col_cnt), DEBUG_NORM(debug_norm), PRINT_NORMALIXATION_BASIS(print_norm), CPU_TOT(cpu), USE_EVERYREL_AS_ARRNORM_BASIS(u_all)
  {};
};
/**
   @typedef tsettings_input_t
   @brief Holds the properties for the blast file.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
**/
typedef struct tsettings_input tsettings_input_t;
#endif
