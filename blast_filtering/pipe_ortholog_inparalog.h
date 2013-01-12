#ifndef pipe_ortholog_inparalog_h
#define pipe_ortholog_inparalog_h
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
#include "types.h"
#include "list_file_parse.h"
#include <stack>
#include "bucket_pipe_binary.h"
#include "taxon_pair.h"
#include "id_simil_list.h"
#include "basket_parse.h"
#include "algo_overlap.h"
/**
   @file
   @class pipe_ortholog_inparalog
   @brief Filters orthologs- and inparalogs in parallel.
   @ingroup pipe_filters
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 15.09.2011 by oekseth (asserts)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
   @date 25.12.2011 by oekseth (cleanup).
**/
class pipe_ortholog_inparalog: public tbb::filter {
 private:
  bool PRINT_NORMALIXATION_BASIS;
  bool DEBUG_NORM;
  loint debug_cnt_possible_orthologs;
  loint debug_cnt_possible_orthologs_inserted;
  loint debug_cnt_possible_inparalogs;
  loint debug_cnt_possible_inparalogs_inserted;
  taxa_t *listTaxa;
  bool MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL;
  char *FILE_BINARY_LOCATION;
  bool PRINT_OVERLAP_VALUES_ABOVE;
  bool use_improved_overlap_algo;
  bool DEBUG_PRINT_DISCARDED_PAIRS;
  id_simil_list listOrtho;  
  list_file_parse_t *listParseData; 
  short int AMINO_LIMIT;
  float max_input_value; // hold the highets input value read
  float MIN_SIMILARITY_LIMIT;
  log_builder_t *log;
  const uint taxon_length;
  const bool INPARALOG_OPERATION;
  const bool USE_EVERYREL_AS_ARRNORM_BASIS;
  const uint n_threads;
  uint list_pair_pos; // The position to get the 'pair block' from
  slock_t lock_set_ps_id;
  bool *in_use; // Holds the thread ids, where a field is set to true if its in use: othervise false
  list_norm_t **l_arrNorm;  
  //  /**! Allocates- and intiates memory for the 'norm_t'structure  */
  //  list_norm_t *init_arrNorm();
  /**! Returns the sim_score of the given index */
  float getSimScore(p_rel_t *buff, uint index);
  //    Returns the similarity score: if zero reciprocal vector is found, returns '0.0'
  float getOuterSimilarity(basket_parse *basket,//p_rel_t *buff_out_in,
			   uint taxon_out, uint taxon_in,
			   uint protein_out, uint protein_in,
			   overlap_t &o_v_right_in, overlap_t &o_v_right_out);
  //! Returns true if the value is above the given limit
  bool aboveInparalogLimit(uint taxon_id, uint protein_id, const float avg);

  //! Takes the average, and inserts it (in accordanse with the settings
  float insert_relation(uint my_id, uint taxon_in, uint taxon_out,
			uint protein_in, uint protein_out,
			float sim_in_out, float sim_out_in);
  ///! Find the relations for this protein towards other proteins belonging to an outer taxon
  void getProteinRelations(const uint my_id, basket_parse *basket_inn, basket_parse *basket_out, uint taxon_out, uint taxon_in, const uint real_protein_in); 
 public:
  //! Builds the set of blocks to be used during the parsing:
  struct taxon_pair *init_taxon_pair(uint taxon_start, uint taxon_length, uint __taxon_length, int n_threads);
  //! The object containing list_file_struct during the building process for each thread id.
  list_file_parse<rel> **l_fileStruct;
  //! De-allocates the data bounded by this class
  void free_data();
  //! De-allocates both the internal temporary objects of type list_file_struct for all of the threads, in addition to de-allocating the memory reserved for this object.
  void free_additional_blocks();
  //! The method of parallisation.
  void* operator()(void* item);
  //! The constructor.
  pipe_ortholog_inparalog(uint _taxon_length,
	      const bool _inparalog_operation, // set to true if its the inparalogs who shall be treated
	      const bool _use_everyrel_as_arrnorm_basis, const uint _n_threads, log_builder_t *_log,
	      id_simil_list &_listOrtho, list_file_parse_t *_listParseData, short int _AMINO_LIMIT,
	      float _max_input_value, float _MIN_SIMILARITY_LIMIT, bool _use_improved_overlap_algo,
	      bool _DEBUG_PRINT_DISCARDED_PAIRS, bool _PRINT_OVERLAP_VALUES_ABOVE,
	      bool _PRINT_NORMALIXATION_BASIS, bool _DEBUG_NORM, taxa_t *listTaxa,
	      bool _MODE_PAIRWISE_OUTPUT_ABC, bool MODE_PAIRWISE_OUTPUT_MCL, char *FILE_BINARY_LOCATION
	      );

#ifdef assert_code
  //! The test function for the private parts of this class.
  void assert_private_parts();
#endif
  //! The main test function for this class  
  static void assert_class(const bool print_info);

};



#endif
