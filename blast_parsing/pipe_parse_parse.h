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
#ifndef pipe_parse_parse_h
#define pipe_parse_parse_h
#include "list_file_parse.h"
#include "tsettings_input.h"
#include "blast_extractors.h"
//#include "module_blast_parsing.h"
#include "string_section.h"
#include "parse_send.h"
#include "parse_send_first.h"
#include "taxon_list_t.h"
#include "protein_vector.h" 
#include "parse_read_blocks.h"
#include "buffer_string_list.h"
#include "log_builder.h"
#include "../configure.h"
#include "list_norm.h"
/**
   @file
   @class pipe_parse_parse
   @brief Parses the block of chars given as input from a previous run.
   @ingroup pipe_parsing parsing_ops
   @remark Class wrapping a tbb container:
   - Reads a block of lines for each call to the operator.
   - To be called by the 'read-file-blocks' method or cia the tbbb operator using pipeline architecture.
   - Returns either a 'parse_send_first' or a 'parse_send' structure with the parsed data.
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 02.12.2011 by oekseth (added possiblity of storing the file on the first run in memory, and getting it from there on the second)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/
class pipe_parse_parse: public tbb::filter {
 private:
  // Below variables used int eh assertion of the code.
  loint reading_file_start_position;
  loint reading_file_length;

  // The following settings are for the parsing of columns in the blast file:
  tsettings_input_t blast_settings;
  uint DEFAULT_NUMBER_OF_COLUMNS_IN_NAME; 
  bool use_improved_overlap_algo;
  bool PRINT_NORMALIXATION_BASIS;
  bool DEBUG_NORM;
  char *FILE_BINARY_LOCATION;
  int INDEX_IN_FILE_FOR_PROTEIN_NAME; // The index in the blast file where the protein name is found
  int INDEX_IN_FILE_FOR_TAXON_NAME; // The index in the blast file where the taxon name is found
  bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE;
  bool USE_BEST_BLAST_PAIR_SCORE; // If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
  char *FILE_INPUT_NAME; // the name (the whole path) of the oinput file
  int CPU_TOT;
  taxa_t *listTaxa;
  char SEPERATOR;
  log_builder_t *log;
  slock_t taxa_upd;
  taxon_list** listProteins;
  bool send_all_data_to_pipe; // default=true: if set to false, writes the data in the second run to the files after completion
  // The below private variables are used controlling (verifying) the operation is followed, breaking the software if run using "undef NDEBUG"
  loint lines_in_file_found_first_read;
  loint lines_in_file_found_last_read;
  loint pairs_in_file_found;
  loint pairs_in_file_sent_to_merge;
  loint total_number_of_pairs_overlapping;
  loint total_number_of_chars_processed_in_first_read;
  loint total_number_of_chars_processed_in_second_read;
 public:
  //! @return the number of lines found in the last read
  uint get_lines_in_file_found_last_read(){return (uint)lines_in_file_found_last_read;}
  //! @return the number of pairs overlapping
  uint get_total_number_of_pairs_overlapping(){return (uint)total_number_of_pairs_overlapping;}
  //! Defines the type of parsing.
  bool FIRST_READ; // updated at 'initHashProtein'
  //! The number of different taxa found in the blast file.
  int taxon_length;
  //! Holds the total number of threads
  //  const int n_threads;
  //! If set updates the **norm structure.
  const bool USE_EVERYREL_AS_ARRNORM_BASIS;
  //! The list of the hashed proteins in the collection  .
  prot_list_t *hashProtein; 
 private:
  prot_list_t *hashTaxa;
 public:  
  //! Holds the thread-id's, where a field is set to true if uts in use, otherwise false
  bool *in_use; 
  /**
     @brief A container holding the parsed blast file, with an element for each thread id.
     @remarks Structure is: [for each of the treahds][a seperate reference]
  **/
  list_file_parse_t **parseData; 
  //! Used to know when a new protein in the collection is met in the second parsing run.
  protein_relation *first_protein;
  //! Holds the set of the highest similarity scores found, one for each thread id, thereby representing the local maximum.
  float *max_sim_value;
  //! Holds the data bout the previous protein in order the ease the finding of the indexes
  protein_vector_t *proteinVector; 
  //! Holds the local values of the **norm collection for each thread id.
  list_norm **local_arrNorm;
 private:
  parse_read_blocks_t *parseBlocks;
  buffer_string_list *stringBuffer;
 public:
  //! Inserted: If in debug-mode, this variable corresponds to the total sum of inserted ovelap values.
  uint debug_sum_inserted_overlap_values;

  //! Found in blast file: If in debug-mode, this variable corresponds to the total sum of inserted ovelap values.
  uint debug_sum_inserted_overlap_values_found;

  //! Returns the number of newlines (rows) found in the file:
  loint debug_get_cnt_newlines(char *buffer_one, char *logical_end);
  /**! Initiates the hash list, whom is used for retrieving the unique id's for
     the proteins, given as a string of chars, in the input file   */
  void initHashProtein(taxon_list_t *listProteins);
  /**! Initiates the global 'taxa_t listTaxa' structure   */
  taxa_t* intializeListTaxa();  
  //  /**! Allocates- and intiates memory for the 'norm_t'structure   */
  //  norm_t **init_arrNorm();
  // Returns the end of the column if the data contains all teh filds: retursn NULL if itdoes not contains all the data
  //  char *setNames(char *start_pos, char *&name, uint index_name, char *&taxon, uint index_taxon, char seperator, char *logical_end);
  // Gets each of the columns
  //  char *getIDColumn(const bool first_column, Parse &p, char *pos_column_start, char *logical_end, const char SEPERATOR);
  //! Gets the overlap data, inserts it into the strcture, and returns an updated position in the string
  char *getOverlapAndUpdate(Parse &p, char *&pos_column_start);
  //! Gets the overlap data, inserts it into the strcture, and returns an updated position in the string
  char *getDistanceUpdate(Parse &p, char *&pos_column_start, float &max_sim_val, char *logical_end, const bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE);

  /**
     @brief If NDEBUG is not set, asserts the parsing operation done.
     @return the number of elements combined for both input argument, and internally stored data.
  **/
  loint assert_parsing(list_file_parse_t *data);

 public:
    /**
     @brief Generate (writes) a log file describing the details of the memory signature for strings stored in memory during parsing
     @remarks Useful for analyzing- optimizing the memory allocation procedures.
  **/
  void log_generate_memory_allocation_overview(uint n_threads) {
    if(stringBuffer) stringBuffer->log_generate_memory_allocation_overview(n_threads);
  }
  //! Defines if summing of the blast file scores for each similar pair is to be used.
  void set_best_blast_pair_score(bool b) {USE_BEST_BLAST_PAIR_SCORE = b;}
  //! Initialises the data structures with the params given.
  void set_parse_blocks_second_read(parse_read_blocks_t *parseBlocks_);//, taxa_t *_listTaxa, int _taxon_length);
  //! @return the arrNorm strcture
  list_norm *getArrNorm(float max_input_value);
  /**! Initiates data used for the second read   */
  taxa_t* initSecondRead(taxon_list_t *listProteins, int &_taxon_length, int updated_cpu_cnt);    
  //! Frees the memory allocated specifially for this class
  list_file_parse_t *free_memory(list_file_parse_t *data); // argument set to NULL if zero merges shall be done

  /**
     @brief Inserts the protein pairs into data structures.
     @param <my_id> The id of the given thread.
     @param <buffer_one> The start position in memory to retrieve data from.
     @param <logical_end> The end position in memory to retrieve data from.
     @param <lines_in_file_found> Sets the number of lines found.
     @param <cnt_overlapping_pairs> The number of pairs overlapping: Used for validation in DEBUG-mode.
     @return The total number of proteins found
     @remarks The steps taken for each row in the input given (ie, the algorithm) is:
     -# Jumps to next line if (a) the line does not contain the two identifying labels or (b) the labels found does not match any found in the previous run through the file (ie,  it's not given as a left-right relation)
     -# Get the 4 indexes for the pair found (taxa and protein label).
     -# Updates the data structure: If the relation is new inserts it, else updates the values given.
     @author Ole Kristian Ekseth (oekseth)
  **/
  uint parse_blast_blocks_data(int my_id, char *buffer_one, char *logical_end, loint &lines_in_file_found, loint &cnt_overlapping_pairs);
  /**
     Find the protein identificators
     @param <my_id> The id of the given thread.
     @param <buffer_one> The start position in memory to retrieve data from.
     @param <logical_end> The end position in memory to retrieve data from.
     @param <lines_in_file_found> Sets the number of lines found.
  **/
  mem_loc parse_blast_blocks_ids(int my_id, char *buffer_one, char *logical_end, loint &lines_in_file_found);

  //! Allocates memory and inititates the data
  void intialize_parseData(int my_id, int block_cnt);
  //! The method of parallisation.
  /*overrride*/void* operator()(void* item);

  //! Constructor
  pipe_parse_parse(loint reading_file_start_position, loint reading_file_length, uint disk_buffer_size, int _taxon_length, const bool _USE_EVERYREL_AS_ARRNORM_BASIS,
		   log_builder_t *_log, char _SEPERATOR, char *_FILE_INPUT_NAME, int _CPU_TOT,
		   bool _USE_LAST_BLAST_CLOMUN_AS_DISTANCE, bool print_norm, bool debug_norm,
		   char *binary_loc, int index_protein, int index_taxon,
		   uint _DEFAULT_NUMBER_OF_COLUMNS_IN_NAME, bool _use_improved_overlap_algo);
  

#ifdef assert_code
  //! The test function for the private parts of this class.
  void assert_private_parts();
#endif
  //! The main test function for this class  
  static void assert_class(const bool print_info);
};
#endif
