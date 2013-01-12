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
#ifndef pipe_parse_merge_h
#define pipe_parse_merge_h
#include "list_file_parse.h"
#include "parse_send.h"
#include "parse_send_first.h"
#include "taxon_list_t.h"
#include "parse_read_blocks.h"
#include "prot_list.h"
#include "taxa.h"
#include "list_file_parse.h"
#include "log_builder.h"
/**
   @file
   @class pipe_parse_merge
   @brief Processes in serial parsed elements of a blast file with task of mering these seperate objects into one consistent object.
   @ingroup pipe_parsing
   @remarks Inputs either a 'parse send_first' or a 'parse_send' structure, containing blocks of information to be merged.
   - The result is merged binary files of type 'file_parse'
   - Class includes methods validating reading of data by generating an empty test-set with a predefined set of newline-chars.
   @author: Ole Kristian Ekseth (oekseth)
   @date 29.12.2010 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/
class pipe_parse_merge: public tbb::filter {
  private:
  bool USE_BEST_BLAST_PAIR_SCORE; // If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
  taxa_t *listTaxa;
  char *FILE_BINARY_LOCATION;
  lint MAX_PARSE_BUFFER_SIZE;
  log_builder_t *log;
  prot_list_t *hashProtein; // The list of the hashed proteins in the collection
  int taxon_length;
  list_file_parse<p_rel> *parse_data_this;
  taxon_list_t *listTaxon;
  list<parse_send_t> lstParse; // holds the list of parsed elements
  list<parse_send_first_t> lstFirstParse; // holds the list of parsed elements
  queue<uint> blocks_end;
  int next_block;
  bool FIRST_READ;
  uint CPU_TOT; // Useful for logging.
 public:
  //! Update shte number of cpu's (for the logs).
  void set_cpu_number(uint c) {CPU_TOT = c;}
  //! Holds the set of the highest similarity scores found: Updated for each conainter reseived during the pipe, thereby representing the global maximum.
  float max_sim_score;
  //! If in debug-mode, this variable corresponds to the total sum of inserted ovelap values.
  uint debug_sum_inserted_overlap_values;
  //! A container holding the parsed blast file.
  parse_read_blocks *parseBlocks; // holds the blocks to read
  //! Holds the overlap values for each protein in the collection.
  overlap_t **arrOverlap;
  //! Defines if summing of the blast file scores for each similar pair is to be used.
  void set_best_blast_pair_score(bool b) {USE_BEST_BLAST_PAIR_SCORE = b;}
  //! @return The parsed blast file object.
  list_file_parse<p_rel> *getListParseData();
  //! @return The object containing the information about the taxa- and their proteins.
  taxon_list_t *getListTaxon();
  //!  Gets the overlap values into teh taxa structure, for later use. Frees the container given as input.
  void getInitTaxaArrOverlap(taxa_t *&listTaxa);
  /**! Called from the outside before the second pipeline is run  */
  void initSecondRead(taxa_t *&_listTaxa, int _taxon_length, prot_list_t *&_hashProtein);  
  //! Cleans the memory for an outer overlap array
  void free_overlap(overlap_t **arr);
  //! Frees the memory for arrOverlap and the 'taxon_list'
  void finalize_class();
  //! Frees the memory for the internal parse_read_blocks object.
  void finalize_parse_blocks();
  //! Finalize an input block when at the second read
  void finalize_block(file_parse<p_rel> **parse_block_c, struct protein_relation parse);
  //! Returns the maximal optimal size of the buffer to use
  mem_loc getMaxBufferSize(uint taxon_length);
  //! Writes the buffers to a file if the amount of data is above the assumed memory whom is required duirng the operation   
  void write_buffers_to_file(bool use_minimum_amount_of_memory);
  //! Merges the data when at the second read
  void merge_data(file_parse<p_rel> **parse_block_c, struct protein_relation rel, const bool debug);    
  //! Builds a consistent list of the overlap scores
  void merge_overlap(overlap_t **overl);  
  //! The method of parallisation.
  /*overrride*/void* operator()(void* item);
  //! The constructor
  pipe_parse_merge(uint disk_buffer_size, int _taxon_length, log_builder_t *_log, lint _MAX_PARSE_BUFFER_SIZE, taxa_t *_listTaxa, char *_binary_loc, uint _CPU_TOT);
    
#ifdef assert_code
  //! The test function for the private parts of this class.
  void assert_private_parts();
#endif
  //! The main test function for this class  
  static void assert_class(const bool print_info);
  
};
#endif
