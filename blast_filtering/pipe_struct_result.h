#ifndef pipe_struct_result_h
#define pipe_struct_result_h
/**
   @file
   @brief Sums up the result of the rows generated for the output files.
   @ingroup output_format
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
#include "../blast_common/enum_mcl.h"
/**
   @struct pipe_struct_result
   @brief Sums up the result of the rows generated for the output files.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
class pipe_struct_result {
 private:
  loint list[mcl_t_size+1];
  loint proteins_processed;
  loint taxon_pair_elements_processed;
  loint cnt_inparalogs;
  loint cnt_total_elements_in_list_file_struct;
 public:
  //! Get the list:
  void get_list(loint *&arg) {
    assert(!arg);
    arg = list ;
  }
  //! Builds the output.
  void print_result(FILE *f) {
    assert(f);
    char *mcl_names[mcl_t_size+1] = { //"..", ".."};
      "inpa",
	"inpa_number",
	"orth_inpa",
	"orth_inpa_number",
	"pair_orth",
	"pair_orth_number",
	"all",
	"all_number",
	"names_index", // holds the names of the proteins, and their corresponding index
	"none" // to be used as value when data is not to be specified
	};
    fprintf(f, "--\t The usage: ");
    for(uint i = 0; i < mcl_t_size; i++) {
      fprintf(f, "%s(%lld) ", mcl_names[i], list[i]);
    }
    if(taxon_pair_elements_processed) {fprintf(f, "taxon_pair_elements_processed(%llu)", taxon_pair_elements_processed);}
    if(proteins_processed) {fprintf(f, " proteins_processed(%llu)", proteins_processed);}
    if(cnt_inparalogs) {fprintf(f, " inpa(%llu)", cnt_inparalogs);}
    if(cnt_total_elements_in_list_file_struct) {fprintf(f, " total(%llu)", cnt_total_elements_in_list_file_struct);} 
    fprintf(f, "\n");
  }

  //! Set the total number of elements, as a reference value for comparing the outputs.
  void set_cnt_total_elements_in_list_file_struct(loint cnt) {cnt_total_elements_in_list_file_struct = cnt;}
  //! Append the inparalog counter.
  void append_inparalog_counter(loint cnt) {cnt_inparalogs += cnt;}
  //! Updates the number of proteins processed.
  void append_meta_data(uint taxon_pair_id, loint length) {
    if(taxon_pair_id > taxon_pair_elements_processed) taxon_pair_elements_processed = taxon_pair_id;
    proteins_processed+=length;
  }
  //! Increments the length for the index given.
  void append_length(uint index, loint length) {
    assert(index < (mcl_t_size+1));
    list[index] += length;
  }
  //! Intitiates the container holding the lengths to empty.
  pipe_struct_result() :
    proteins_processed(0), taxon_pair_elements_processed(0),
    cnt_inparalogs(0), cnt_total_elements_in_list_file_struct(0)
    {
      for(uint i = 0; i < (mcl_t_size+1); i++) list[i] = 0;
    }
};

/**
   @typedef pipe_struct_result_t
   @brief Sums up the result of the rows generated for the output files.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
typedef pipe_struct_result pipe_struct_result_t;

#endif
