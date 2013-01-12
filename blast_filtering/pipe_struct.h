#ifndef pipe_coortholog_result_h
#define pipe_coortholog_result_h
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
#include "log_builder.h"
#include "pipe_t.h"
#include <stack>
#include "taxa.h"
#include "list_file_parse.h"
#include "mcl_format.h"
#include "taxon_pair.h"
#include "bucket_norm.h"
#include "id_simil_list.h"
#include "pipe_struct_result.h"

#ifndef NDEBUG
/**
   Uses validating the elements inserted correspongning to our expectations:
   @remarks Intended to be handled sperately for each thread id.
**/
typedef struct meta_pipe_struct {
  //! The number of names ('root values') found:
  uint cnt_names;
  //! The number of inparalogs found:
  uint cnt_inpa;
  //! The number of orthologs found:
  uint cnt_ortho;
  //! The number of co-orthologs found:
  uint cnt_co_orth;
  //! The (empty) constructor:
  meta_pipe_struct() : cnt_names(0), cnt_inpa(0), cnt_ortho(0), cnt_co_orth(0) {};
} meta_pipe_struct_t;
#endif

/**
   @file
   @class pipe_struct
   @brief Either builds co-orthologs or builds the strings for the result file.
   @ingroup pipe_filters
   @return A 'bucket_norm object.
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
   @date 25.12.2011 by oekseth (cleanup).
**/
class pipe_struct: public tbb::filter {
 private:
#ifndef NDEBUG
  meta_pipe_struct_t *lst_elements_evaluated;
  meta_pipe_struct_t elements_expectation;
#ifdef USE_MPI
  int myrank;
#endif
#endif
  int number_of_nodes;
  //! For the abc file: if set, the data out pairwise in stead of as in a row
  bool MODE_PAIRWISE_OUTPUT_ABC; 
  //! for the mcl file: if set, the data out pairwise in stead of as in a row
  bool MODE_PAIRWISE_OUTPUT_MCL; 
  // The following variables decides what output is to be printed
  bool PRINT_IN_ABC_FORMAT;
  bool PRINT_IN_MCL_FORMAT;
  bool SORT_ABC_DATA; // if true, sorts teh abc files before outprint
  bool PRINT_NORMALIXATION_BASIS;
  bool DEBUG_PRINT_DISCARDED_PAIRS;
  bool DEBUG_NORM;
  int taxon_length;
  // --
  bool MODE_PAIRWISE_OUTPUT; // If set, prints the data out pairwise in stead of as in a row
  bool MODE_INTEGER_OUTPUT; // If set, prints the data out as integers in stead of as in a row
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT;
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT;
  id_simil_list listOrtho;  
  taxa_t *listTaxa;
  stack_rel *stackRel;// Holds the co orthologs for every protein in a stack; intialized in 'ortho_set' 
  // Data arrays to intialize and fill
  list_file_struct_t *listStructData; 
  list_file_parse_t *listParseData;
  float max_input_value;
  // --
  log_builder_t *log;
  uint n_threads;
  bool *in_use; // Holds the ps'es, where a field is set to true if uts in use, otherwise false
  //  const int g_protein_length;
  const int g_taxon_length;
  const uint SIZE_BUFFER; //the buffersize; when execided this limit, it is increased by the same amount
  uint list_pair_pos;
  const pipe_t PIPE_TYPE;
  const bool USE_EVERYREL_AS_ARRNORM_BASIS;
  float **arrAvgNorm; // holds the normative values
  uint _LENGTH_KEY_NAMES;
  uint _WIDTH_DISTANCE;
  taxon_pair *listPair; // Holds the blocsk to be used for the parsing
  slock_t lock_set_ps_id;
  list_norm_t **l_arrNorm;

  static const bool OPT_GET_FILE_LOCALLY = true;
  mcl_format_t **mclData;
  pipe_struct_result debug_dump_result;
  /**! Allocates- and intiates memory for the 'norm_t'structure  */
  list_norm_t *init_arrNorm();
  //! Produces a row for the protein given
  void produce_row(uint my_id, uint protein_in, uint taxon_in);
  //! Builds the array for the normalization values doing the averaging
  void build_normArr(const uint taxon_length, list_norm_t* arrNorm);
  /**! Inserts the similarity score into the stack*/
  void insert_into_stack_rel(uint my_id, uint taxon_in, uint taxon_out, uint world_in, uint world_out, float sim_score);

  ///! Insert the orthologs for the for the pair "protein_in <==> protein_out"
  void insertOrthologs(uint my_id, rel_t *buff_in_in, uint taxon_in, uint taxon_out, uint world_in, uint world_out, uint protein_in, uint protein_out);
  /**! Insert the inpa_inpa_relations   */
  void insertInpaInpaRelations(uint my_id, rel_t *buff_in_in, uint taxon_in, uint taxon_out/*, uint world_in*/, uint world_out, uint protein_in, uint protein_out);
 public:   
  //! Enables outprint of discarded pairs during the filtering procedure.
  void set_DEBUG_PRINT_DISCARDED_PAIRS(bool d) {DEBUG_PRINT_DISCARDED_PAIRS=d;}
  //! Opens a buffer and dumps it to memory if the space  for it is allocated
  void getFromBufferToStruct(const uint taxon_length, uint max_buffer_size, const bool only_inparalogs);
  //! Clears the memory allocated for this thread
  void finalize_memory(const uint taxon_length);
  //! The method of parallisation.
  /*orverride*/void* operator()(void* item);
  //! De-allocates the memory reserved for this object.
  void free_mem();
  /**
     @brief Constructor for the class.
     @remarks Operation decided when setting the 'PIPE_TYPE' variabel
     (PIPE_TYPE == DUMP):     Writes the data to a matrix, for further processing by another program
     (PIPE_TYPE == INPA_ORTH) Creates inparalogs based upon orthologs.
     (PIPE_TYPE == INPA_INPA) Creates inparalogs based upon orthologs' inparalogs.
  */
  pipe_struct(uint _nthread, uint taxon_start, uint taxon_length, pipe_t type, bool _USE_EVERYREL_AS_ARRNORM_BASIS, list_norm_t *normArr, log_builder_t *_log,
	      bool _MODE_PAIRWISE_OUTPUT, bool _MODE_INTEGER_OUTPUT, bool _DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT,
	      bool _DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT, id_simil_list &_listOrtho, taxa_t *_listTaxa, stack_rel *&_stackRel, list_file_struct_t *&_listStructData, list_file_parse_t *_listParseData, float _max_input_value,
	      bool _MODE_PAIRWISE_OUTPUT_ABC, bool _MODE_PAIRWISE_OUTPUT_MCL, bool _PRINT_IN_ABC_FORMAT,
	      bool _PRINT_IN_MCL_FORMAT, bool _SORT_ABC_DATA, bool _PRINT_NORMALIXATION_BASIS, bool _DEBUG_NORM, float **&arr_avgNorm
	      );


#ifdef assert_code
  //! The test function for the private parts of this class.
  void assert_private_parts();
#endif
  //! The main test function for this class  
  static void assert_class(const bool print_info);
};


#endif
