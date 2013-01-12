#ifndef id_simil_list_h
#define id_simil_list_h
/**
   @file
   @struct id_simil_list
   @brief A linked list class
   @ingroup blastfile_container
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
#include "types.h"
#include "taxa.h"
#include "rel.h"
#include "log_builder.h"
#include "list_file_parse.h"
/**
   @struct id_simil_list
   @brief A linked list class
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by oekseth (init).
   @date 25.12.2011 by oekseth (cleanup).
**/
struct id_simil_list {
private:
  rel_t** arr;
  uint arr_size;
  taxa_t *listTaxa;
  // Below are config-settings for this scope of the processing.
  bool DEBUG_PRINT_DISCARDED_PAIRS;
  uint taxon_length;
  //  uint index_size;
  static const uint BASE_ARR_STEP = 10;
  public:
  //! @return the list, an sets the size:
  rel_t **get_internal_list(uint &_arr_size) {
    assert(!_arr_size);
    _arr_size = arr_size; return arr;
  }
  //! Test if condidtion is passed.
  static void test_condition_and_if_not_abort(bool condition, const int line, const char *function) {
    if(!condition) {
      fprintf(stderr, "!!\tWas not able allocating data due to memory constraints. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", line, __FILE__, function);
      exit(2);
    }
  }
  //! @return true if data is set.
  bool has_data() {return (arr && arr_size);}
  //! Enables printing of discarded pairs. 
  void set_DEBUG_PRINT_DISCARDED_PAIRS(bool d, uint l, taxa_t *li) {DEBUG_PRINT_DISCARDED_PAIRS = d; taxon_length=l; listTaxa = li;}
  /**
     @brief   Initiates a row of the '**arr'
     @remark  At is starting point not malloced before at this point (after the
     iteration of the orthologs, but before the iteration of the inparalogs).
  */
  void initGlobalRow(uint index_in, uint size);

  /**
     @brief   Initiates a row of the '**arr'
     @remark  At is starting point
     not malloced before at this point (after the iteration of the orthologs,
     but before the iteration of the inparalogs).
  */
  static void initLocalRow(uint index_in, uint size, rel_t **&arr);

  //! Returns the size of the global array at the specified index
  uint getGlobalSize(uint index_in);

  //! Returns the size of the local array (given as parameter) at the specified index
  static uint getLocalSize(uint index_in, rel_t **&arr);

  /**! Allocates, and intiates the global array to the length specified */
  void initGlobalArr(uint arr_size, uint index_size);

  /**! Allocates, and intiates the array, given as input parameter, to the length specified */
  static void initLocalArr(rel_t **&arr, uint arr_size, uint index_size);
  
  //! Frees the memory allocated for the global arr
  void freeGlobalArr();
  
  //! Frees the memory allocated for the global arr
  void freeLocalArr(rel_t **&arr, uint arr_size);

  //! Inserts an elements into the global array
  void insertGlobalElement(uint protein_in, uint world_protein_out, float sim_score);

/**
   @brief   Inserts an element into a row of the '**arrOrtho'
   @remark  Is first needed in the 'ortho_set_bin', therefore
   not malloced before at this point (after the iteration of the orthologs,
   but before the iteration of the inparalogs).
*/
  void insertLocalElement(rel_t **&arr, uint index_in, uint world_protein_out, float sim_score);

/**
   @param <protein_in> The protein-id to retrieve data for.
   @param <size> Sets the size of the array returned
   @return  Returns a list (array) of the proteins 'protein_in'
   has ortholog pairs with
*/
  rel_t* getGlobalRow(uint protein_in, uint &size);


  /**
     @param <arr> The list to use for the retrieval.
     @param <protein_in> The protein-id to retrieve data for.
     @param <size> Sets the size of the array returned
     @return  Returns a list (array) of the proteins 'protein_in' from the array given as input parameter has ortholog pairs with
  */
  static rel_t* getLocalRow(rel_t **&arr, uint protein_in, uint &size);

  //! Returns true if the index given as parameter has data set for ti
  bool hasLocalData(rel_t **&arr, uint index_in);

  //! Returns true if the index given as parameter has data set for ti
  bool hasGlobalData(uint index_in);

  //! Returns true if 'index_out' has 'index_in' as a relation
  bool hasRelation(rel_t **&arr, uint index_in, uint index_out);
  //! Returns true if 'index_out' has 'index_in' as a relation
  bool hasRelation(uint index_in, uint index_out) {
    if(arr) {
      if(index_out < arr_size) { 
	return hasRelation(arr, index_in, index_out);
      }
    }
    return false;
  }

  /**
     @brief Initalizes the data to zero
     @remark Does not de allocate the memory
  */
  void clearLocalData(rel_t *&arr);

  /**
     @brief Creates a bi partite matching
     @remarks Assumes every key (the array index and the index stored as an elemnt in the array follows the same key-code-rule
  */
  void execInternLocalReciProc(rel_t **&arr, uint arr_size);

  /**
     @brief Creates a bi partite matching
     @remarks Assumes every key (the array index and the index stored as an elemnt in the array follows the same key-code-rule
  */
  void execInternGlobalReciProc(uint arr_size) {
    execInternLocalReciProc(arr, arr_size);}
  
  //! Prints the array givens as input argument
  void printLocalArr(FILE *f, rel_t **arr, const bool use_names, taxa_t *listTaxa, uint taxon_length, uint arr_size/*,  bool to_file*/);
    
  //! Prints the global array:
  void printGlobalArr(FILE *f, const bool use_names, taxa_t *listTaxa, uint taxon_length) {
    printLocalArr(f, arr, use_names, listTaxa, taxon_length, arr_size);
  }

  //! De-allocates the memory for this object.
  void free_memory() {if(has_data()) freeGlobalArr();}
   /**
     @brief Prints the memory consumption for the object.
     @remarks For anaylsing the memory fingerprint.
  **/
  void print_getTotal_memoryConsumption_for_object(FILE *f) {
    assert(f);
    loint numbers_reserved = 0;
    if(arr) {
      for(uint i = 0; i < arr_size; i++) {
	numbers_reserved += getLocalSize(i, arr);
      }
    }
    const loint size_reserved = sizeof(id_simil_list) + sizeof(rel_t)*numbers_reserved;
    const float avg_cnt = (float)numbers_reserved/arr_size;
    //    const float fill_factor = (float)size_used / size_reserved;
    const float size_gb = (float)size_reserved/(1024*1024*1024);
    char end_char = ' ';
    // PRints the data, but in order to avoid clutter, newline is avoided if the goal is printing the data to the log file:
    if((f==stdout) || (f == stderr)) end_char = '\n';
    fprintf(f, "- Has %lld B of memory =~ %.5f GB => averaged_cnt(%.4f)%c", size_reserved, size_gb, avg_cnt, end_char); 
  }
  /**
     @return the number of pairs
     @remarks
     # The protein-pair distance is found at 'arr[protein_in][0].distance'
     # A protein-pair is inserted at index [protein_in][0].ind_out+1'
     # This gives the number of pairs: 'obj[0].ind_out+1'
  **/
  static uint get_number_of_pairs(rel_t *obj) {
    if(obj) return obj[0].ind_out; 
    else return 0;
  }
/**
   @return the sum of distances for pairs
   @remarks 
   # Used verifying that the content of two objects are equal
   # Practical when MPI sending/receving is to be verified.
**/
  static float get_total_sum_of_distances_for_pairs(rel_t **arr, uint arr_size);

/**
   @return the number of pairs
   @remarks Use a '+arr_size' offset to get the size, due to the first reference object.
**/
  static uint get_total_number_of_pairs(rel_t **arr, uint arr_size);
/**
   @return the number of pairs
   @remarks Use a '+arr_size' offset to get the size, due to the first reference object.
**/
  uint get_total_number_of_pairs() {return get_total_number_of_pairs(arr, arr_size);}

/**
   @return the number of pairs
   @remarks Use a '+arr_size' offset to get the size, due to the first reference object.
**/
  static uint get_total_number_of_pairs(rel_t **arr, uint arr_size, uint taxon_start, uint taxon_end, taxa *listTaxa, uint taxon_length, bool *list_of_nodes_taxa_responsilibties);
/**
   @return the number of pairs
   @remarks Use a '+arr_size' offset to get the size, due to the first reference object.
**/
  uint get_total_number_of_pairs(uint taxon_start, uint taxon_end, bool *list_of_nodes_taxa_responsilibties) {
    return id_simil_list::get_total_number_of_pairs(arr, arr_size, taxon_start, taxon_end, listTaxa, taxon_length, list_of_nodes_taxa_responsilibties);
  }


/**
   @return the sum of distances for pairs
   @remarks 
   # Used verifying that the content of two objects are equal
   # Practical when MPI sending/receving is to be verified.
**/
    float get_total_sum_of_distances_for_pairs() {return get_total_sum_of_distances_for_pairs(arr, arr_size);}
#ifdef USE_MPI
  /**
     @brief Performs intra-node building of reciprocal ortholog pairs
     @remarks Core code located in file "mpi_id_simil_list.h", making the seperation of the strictly MPI code verbouse.   
  **/
  void mpi_tx_rx_recip_tx_rx(int taxon_length);
  /**
     @brief Sends the co-orthologs to all the nodes.
     @remarks Core code located in file "mpi_id_simil_list.h", making the seperation of the strictly MPI code verbouse.   
  **/
  static void mpi_send_co_orthologs_accross_nodes(taxa *listTaxa, int taxon_length, stack_rel *stackRel, list_file_parse<rel> *listStructData);
#endif

  //! The constructor.
  id_simil_list() : arr(NULL), arr_size(0), listTaxa(NULL), DEBUG_PRINT_DISCARDED_PAIRS(false), taxon_length(0) {};
  //! The constructor.
  id_simil_list(uint _arr_size, taxa_t *_listTaxa) :
    arr(NULL), arr_size(_arr_size), listTaxa(_listTaxa), DEBUG_PRINT_DISCARDED_PAIRS(false), taxon_length(0)
  {
    if(arr_size)  initGlobalArr(arr_size, BASE_ARR_STEP);

  }
  //! The main test function for this class  
  static void assert_class(const bool print_info);
};


/**
   @brief A linked list class
   @ingroup blastfile_container
**/
typedef id_simil_list id_simil_list_t;
#endif
