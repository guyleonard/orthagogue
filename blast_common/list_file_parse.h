/**
   @file
   @brief Stores the blast file with information about overlap-values.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth).
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
#ifndef list_file_parse_h
#define list_file_parse_h
#include <string.h>
#include "log_builder.h"
#include "file_parse.h"
#include "prot_list.h"
#include "protein_relation.h"

/**
   @class list_file_parse
   @brief Stores the blast file with information about overlap-values.
   @tparam <T> A rel class for the relations_list object, ie, the implemented p_rel or rel classes
   @ingroup blastfile_container
   @remark: Serves as a wrapper matrix for the squares of the data, each squares
   representing a combination of two taxa, as shown in the documentation
   included in this code package.
   @author Ole Kristian Ekseth (oekseth).
   @date 12.01.2011 by oekseth (initial).
   @date 08.09.2011 by oekseth (asserts).
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of this class as a libary)
**/ 
template<class T> class list_file_parse {
 private:
  int taxon_length; // The number of taxon in the collection:
  taxa_t *listTaxa; // Stores a global pointer to this, but never changes it.
  bool MODE_PAIRWISE_OUTPUT_ABC;
  bool MODE_PAIRWISE_OUTPUT_MCL;
  char *FILE_BINARY_LOCATION;
  prot_list_t *hashProtein; // Stores a global pointer to this, but never changes it: Used to get the length of each of the proteins
  file_parse<T> ***list; // The actual data.
  const bool send_all_data_to_pipe; // default=true: if set to false, writes the data in the second run to the files after completion
  int my_id;
  bool USE_BEST_BLAST_PAIR_SCORE; // If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
#ifdef USE_MPI
  //! Identifies the taxa this object is responsible for, ie, the node this object belongs to; set from "mpi_transfer_list_file_parse" object.  
  bool *list_of_nodes_taxa_responsilibties;
#endif
 public:
  
  //! The maximum number of elements to reside in the buffer before they are written to a file.
  uint MAX_BUFFER_SIZE; 
  //! @return get MAX_BUFFER_SIZE
  uint get_MAX_BUFFER_SIZE() {return MAX_BUFFER_SIZE;}
#ifdef USE_MPI
  //! Sets the list identifing taxa this object is responsible for, ie, the node this object belongs to; set from  "mpi_transfer_list_file_parse" object.  
  void set_list_of_nodes_taxa_responsibilities(bool *&list, uint length) {
    assert(!list_of_nodes_taxa_responsilibties);
    assert(length == (uint)taxon_length);
    assert(list);
    list_of_nodes_taxa_responsilibties = list;
    list = NULL; 
  }
  //! Sets the list identifing taxa this object is responsible for, ie, the node this object belongs to; set from  "mpi_transfer_list_file_parse" object.  
  void copy_list_of_nodes_taxa_responsibilities(list_file_parse<p_rel> *obj) {
    //    uint length = 0;
    assert(!list_of_nodes_taxa_responsilibties);
    //    obj->take_and_remove_list_of_nodes_taxa_responsilibties(list_of_nodes_taxa_responsilibties, length);
    assert(taxon_length);
    list_of_nodes_taxa_responsilibties = new bool[taxon_length];
    assert(list_of_nodes_taxa_responsilibties);
    log_builder::test_memory_condition_and_if_not_abort(list_of_nodes_taxa_responsilibties!=NULL, __LINE__, __FILE__, __FUNCTION__);    
    memcpy(list_of_nodes_taxa_responsilibties, obj->get_list_of_nodes_taxa_responsilibties(), taxon_length*sizeof(bool));
  }
  //! @return the list identifing taxa this object is responsible for, ie, the node this object belongs to; set from "mpi_transfer_list_file_parse" object.  
  bool *get_list_of_nodes_taxa_responsilibties(){return list_of_nodes_taxa_responsilibties;}
  /**
     @brief Take the list from this object, and insert it into the parameters.
     @param<list> Sets the list identifing taxa this object is responsible for, ie, the node this object belongs to.
     @param<length> If correct performance, the variable is set to the number of taxa.
  **/
  void take_and_remove_list_of_nodes_taxa_responsilibties(bool *&list, uint &length){
    //! Expects the internal "list_of_nodes_taxa_responsilibties" to be set when this function is called.
    assert(!list);
    assert(!length);
    assert(list_of_nodes_taxa_responsilibties);
    assert(taxon_length);
    //! Sets the data:
    list = list_of_nodes_taxa_responsilibties;
    length = taxon_length;
    list_of_nodes_taxa_responsilibties = NULL;
  }
  //! @return true if the given inner taxon is "myranks" responsibility.
  bool taxon_is_myranks_responsibility(uint taxon_id) {
    //! Includes the follogin lines to show the expectations of use.
    assert(taxon_id < (uint)taxon_length); 
    assert(list_of_nodes_taxa_responsilibties); // SHolud be set.
    if(taxon_id < (uint)taxon_length) return list_of_nodes_taxa_responsilibties[taxon_id];
    else return false;
  }
#else
  //! Serves to hide complexity from the callers point of view.
  void set_list_of_nodes_taxa_responsibilities(list_file_parse<T> *obj) {
    if(false && obj) ;  // In order to hide the non-usage of the variable.
  }
  //! @return true whatsoever: Serves to hide complexity from the callers point of view.
  bool taxon_is_myranks_responsibility(uint taxon_id) {
    if(false && taxon_id) ; // In order to hide the non-usage of the variable.
    return true;
  }
#endif
  
  //! @return the length of the taxa-collection
  int get_taxon_length(){return (int)taxon_length;}
  //! @return the object at the given position; if not set returns NULL.
  file_parse <T> *get_element(uint taxon_in, uint taxon_out) {
    if(
       (taxon_in < (uint)taxon_length)
       && (taxon_out < (uint)taxon_length)
       && (list && list[taxon_in] && list[taxon_in][taxon_out])
       ) {
      return list[taxon_in][taxon_out];
    } else return NULL;
  }

  /**
     @brief Sets an element of type file_parse<T> without allocating memory:
     @param <taxon_in> The leftmost id of the taxon.
     @param <taxon_out> The rightmost id of the taxon.
     @param <element>   The element to insert.
     @warning Does not allocate memory, ie, element must not be deleted after inserten (from the calles point of view).
  **/
  void set_element_unsafe(uint taxon_in, uint taxon_out, file_parse<T> *&element) {
    assert(taxon_in < (uint)taxon_length);
    assert(taxon_out <(uint)taxon_length);
    assert(list != NULL);
    assert(list[taxon_in] != NULL);
    assert(list[taxon_in][taxon_out] == NULL); // For safeness should not be set.
    list[taxon_in][taxon_out] = element;
  }
  //! Test if condidtion is passed.
  static void test_condition_and_if_not_abort(bool condition, const int line,  const char *function) {
    if(!condition) {
      fprintf(stderr, "!!\tWas not able allocating data due to memory constraints. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", line, __FILE__, function);
      exit(2);
    }
  }

  //! @returns the listTaxa variable
  taxa_t *get_listTaxa(int &_length) {_length = taxon_length; return listTaxa;}
  //! @return the FILE_BINARY_LOCATION variable
  char *get_FILE_BINARY_LOCATION(){return FILE_BINARY_LOCATION;}
  //! @return the USE_BEST_BLAST_PAIR_SCORE variable
  bool get_USE_BEST_BLAST_PAIR_SCORE(){return USE_BEST_BLAST_PAIR_SCORE;}
  //! Sets the protein list of the class
  void setProtList(prot_list_t* prot) {hashProtein = prot;}
  //! @returns true if it has data about list
  bool has_data() {return (list != NULL);}
  /**! @return true if data is registered for the given protein in the specified matrix position  */
  bool has_data(uint taxon_in) {
    assert(list);
    assert(list[taxon_in]);
    assert(taxon_in < (uint)taxon_length);
    bool has_data_set = true;
    for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
      if(list != NULL && list[taxon_in]){
	if(!has_data(taxon_in, taxon_out)) {
	  has_data_set = false;
	}
      } else {
	fprintf(stderr, "!!\tHas no data set for [in=%u][out=%u...%d] (found at line %d in file %s in method %s). Please notify the developer at oekseth@gmail.com is this message is seen!\n", taxon_in, 0, taxon_length, __LINE__, __FILE__, __FUNCTION__);
	assert(false);
	has_data_set = false;
      }
    }
    return has_data_set;
  }

  //! @returns true if the structure has data at the position in the matrix
  bool has_data(uint taxon_in, uint taxon_out) {
    assert(taxon_in < (uint)taxon_length);
    assert(taxon_out < (uint)taxon_length);
    if(list && list[taxon_in] && list[taxon_in][taxon_out]) {
      return list[taxon_in][taxon_out]->has_data();
    }
    return false;
  }
  /**! @return true if data is registered for the given protein in the specified matrix position  */
  bool has_data(uint taxon_in, uint taxon_out, uint protein_in) {
    if(list && list[taxon_in] && list[taxon_in][taxon_out]){
      if(taxon_in < (uint)taxon_length && taxon_out < (uint)taxon_length) {
	return list[taxon_in][taxon_out]->has_data(protein_in);
      }
    } return false;
  }

 private:
  //! @return a visual the identifier for the tempalte fo this class
  const char *get_string_template(rel *obj) {
    if(false && obj) ; // Hides the argument.
    return "rel";
  }
  //! @return a visual the identifier for the tempalte fo this class
  const char *get_string_template(p_rel *obj) {
    if(false && obj) ; // Hides the argument.
    return "p_rel";
  }
 public:
  /**
     @brief Prepares the scheduling.
     @param <list_length> The total length of the input list.
     @param <thread_list> Holds the overview of the taxon-id to process; list is deleted in this function when the alst taxon is processed.
     @param <current_position> The internal index to access in the "thread_list" list.
     @param <taxon_start> The first taxon to work on in the iterval, given as [taxon_start, taxon_end -1]
     @param <taxon_end>   The last taxon to work on in the iterval, given as [taxon_start, taxon_end -1]
     @param <is_ortholog> If required ortholog data are in file(s), reads the data before processing, and removes them at the process' end (due to the fact that the "parsing elelments this involves are no longer of interest).
     @param <is_inpa> If required inparalog data are in file(s), reads the data before processing, and removes them at the process' end (due to the fact that the "parsing elelments this involves are no longer of interest).
     @param <is_co_orth>   If required co-ortholog data are in file(s), reads the data before processing, and removes only the inparalog data at the process' end (due to the need of them during generation of result file).     
     @return the interval to work on, and esnures that the given interval resides in memory.
     @remarks Used for both p_rel and rel objects, ie, in front of both parsing- and filtering.
  **/
  bool get_interval(uint &list_length, uint *&thread_list, uint &current_position, uint &taxon_start, uint &taxon_end, bool is_ortholog, bool is_inpa, bool is_co_orth) { 
    //! Internal tests verifying correct usage:
    if(is_ortholog) {
      assert(!is_inpa);
      assert(!is_co_orth);
    }
    if(is_inpa) {
      assert(!is_co_orth);
      assert(!is_ortholog);
    }
    if(is_co_orth) {
      assert(!is_inpa);
      assert(!is_ortholog);
    }
    
    taxon_start = 0, taxon_end = 0; // Initializes

    //! If software has data to read, the list is set:
    if(thread_list) {
      // Internal tests verifying correct usage:
      if(!current_position) {
	assert(!taxon_start);
      }
      if(current_position) { // Deallocates earlier allocations
	const uint old_taxon_start = thread_list[current_position-1];
	const uint old_taxon_end   = thread_list[current_position]; // The next position in the interval.
	for(uint i = old_taxon_start; i < old_taxon_end; i++) {
	  if(is_ortholog) {
	    file_parse<T>::close_non_inparalogs(list, true, old_taxon_start, old_taxon_end); 
	  } else if(is_inpa) {
	    file_parse<T>::close_For_inparalogs(list, true, old_taxon_start, old_taxon_end);
	  } else if(is_co_orth) {
	    file_parse<T>::close(list, true, old_taxon_start, old_taxon_end);
	  }
	}
      }	    
      if(current_position < list_length) { // Has more data to read.
	taxon_start = thread_list[current_position++];
	taxon_end   = thread_list[current_position]; // The next position in the interval.
	// Gets data into the buffer:
	if(is_ortholog) {
	  dump_buffer_if_above_limit(UINT_MAX, /*is_inparalog=*/false, taxon_start, taxon_end, false);
	} else if(is_inpa) {
	  dump_buffer_if_above_limit(UINT_MAX, /*is_inparalog=*/true, taxon_start, taxon_end, false);
	} else if(is_co_orth) {
	  // Need all the relations in memory, in order to verify the "edges".
	  include_all_relations(taxon_start, taxon_end);
	}
      } else {      
	// Deletes the allocation list if the end is reached:
	if(thread_list) {delete [] thread_list; thread_list = NULL; list_length = 0;}
      }
      return true;
    } else {
      if(is_ortholog) {
	file_parse<T>::close_non_inparalogs(list, true, taxon_start, taxon_end); 
      } else if(is_inpa) {
	file_parse<T>::close_For_inparalogs(list, true, taxon_start, taxon_end);
      } else if(is_co_orth) {
	//! Need the inparalogs when building the result file:
	file_parse<T>::close_non_inparalogs(list, true, taxon_start, taxon_end); 
	//file_parse<T>::close(list, true, taxon_start, taxon_end);
      }
      return false;
    }
  }
  /**
     @brief produces a sheduling list.
     @param <list_length> Holds the overview of the taxon-id to process; list is deleted in this function when the alst taxon is processed.
     @param <is_ortholog>   If one of the outer elements (required in the ortholog processing) is in a file, sets the "whole span" as a seperate "pipe job". 
     @param <is_inpa> If an element is found having data in a file, appends it seperate to the list.
     @param <is_co_orth>   If one of the outer elements (required in the co-ortholog processing) is in a file, sets the "whole span" as a seperate "pipe job". 
     @return a sheduling list to be used as input for "static bool get_interval(...) in this class.
     @remarks The procedure includes the considerations:
     # Build new blocks for "piping" when a taxon-pair is found to have data reciding in file(s).
     # Evaluates also reciprocal relations.
     # To be used in conjunction with method "init_taxon_pair(..)", found in class "taxon_pair", ie, all updates in this descrbied function, must be reflected in the taxon_pair function stated.
     @author: Ole Kristian Ekseth (oekseth).
  **/
  uint *get_thread_allocations_for_taxa(uint &list_length, bool is_ortholog, bool is_inpa, bool is_co_orth) {
    //! Internal tests verifying correct usage:
    if(is_ortholog) {
      assert(!is_inpa);
      assert(!is_co_orth);
    }
    if(is_inpa) {
      assert(!is_co_orth);
      assert(!is_ortholog);
    }
    if(is_co_orth) {
      assert(!is_inpa);
      assert(!is_ortholog);
    }
    if(list && taxon_length) {
      //! Declaration of the return variable:
      const uint taxa_list_allocations_size = taxon_length + 3;
      uint *taxa_list_allocations = new uint[taxa_list_allocations_size];
      log_builder::test_memory_condition_and_if_not_abort(taxa_list_allocations!=NULL, __LINE__, __FILE__, __FUNCTION__);
      memset(taxa_list_allocations, 0, sizeof(uint)*taxa_list_allocations_size);
      uint number_of_blocks = 0; // Hold the number of blocks processed.
      //! Starts at the first taxon:
      taxa_list_allocations[number_of_blocks] = 0; number_of_blocks = 1;
      //! Collects information about those taxa having data stored in files:
      for(int taxon_in = 0; taxon_in < taxon_length; taxon_in++) {
	bool taxon_has_data_in_file = false; // End the operation at the first occurence of data residing in memory.

	//! Only processes if the given "inner" taxon is "myranks" reposponisiblity; when MPI is not used, it's always ;)
	if(taxon_is_myranks_responsibility(taxon_in)) {
	  //! If below "taxon_has_data_in_file" is set, implies file content first must be extracted, ie, a new block is generated:
	  if(is_ortholog || is_co_orth) {
	    for(int taxon_out = 0; taxon_out < taxon_length && !taxon_has_data_in_file; taxon_out++) {
	      if(has_data(taxon_in, taxon_out) && list[taxon_in][taxon_out]->get_number_of_objects_in_file() ) {
		taxon_has_data_in_file = true;
	      }
	      if(list[taxon_out] && has_data(taxon_out, taxon_in) && list[taxon_out][taxon_in]->get_number_of_objects_in_file() ) {
		taxon_has_data_in_file = true;
	      }
	    }
	  } else if(is_inpa) {
	    if(has_data(taxon_in, taxon_in) && (list[taxon_in][taxon_in]->get_number_of_objects_in_file())) {
	      taxon_has_data_in_file = true;
	    }
	  } 
	  if(taxon_has_data_in_file) {
	    //! Generates a new block:
	    number_of_blocks++; // Increments in ordre to get a new block:
	    assert(number_of_blocks < taxa_list_allocations_size);
	    taxa_list_allocations[number_of_blocks] = taxon_in;
	  } else {
	    //! This inner taxon added as an continous block:
	    taxa_list_allocations[number_of_blocks]++;
	  }
	} else {
	  //! It's not our responsibility: In accordance with expectations in method "init_taxon_pair(..)",
	  //  located in class "taxon_pair", we add this as an "continous" block:
	  taxa_list_allocations[number_of_blocks]++;
	}
      }
      //! Sets "the end" mark, and fixes the length of the list:
      taxa_list_allocations[number_of_blocks] = taxon_length;
      taxa_list_allocations[number_of_blocks+2] = taxon_length;
      //! Updates the lwngth, to be used by the calling function:
      list_length = number_of_blocks; 
#ifndef NDEBUG
      //! Verifies that all the data are included into the container
      uint sum_taxa = 0;
      for(uint i = 0; i < number_of_blocks; i++) {
	const uint taxon_start = taxa_list_allocations[i];
	const uint taxon_end = taxa_list_allocations[i+1];
	sum_taxa += (taxon_end - taxon_start);	
      }
      //
      if(!(sum_taxa == (uint)taxon_length)) {
	log_builder::throw_warning(software_error, __LINE__, __FILE__, __FUNCTION__,"building of allocation list");
      //	printf("sum_taxa(%u) && list_length set to %u at line %d in file %s.\n", sum_taxa, list_length, __LINE__, __FILE__); 
	assert(sum_taxa == (uint)taxon_length);
      }
#endif
      
      return taxa_list_allocations;
    } else {
      list_length = 0;
      return NULL;
    }
  }
  //! @returns the length of the protein specified
  uint getProteinLength(uint in, uint out, int protein_in) {
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    assert(list);
    assert(list[in]);
    if(list[in][out]) return (uint)list[in][out]->getProteinLength(protein_in);
    return 0;
  }  

  /**
     @brief Gives the number of pairs in a subset of the blastp file.
     @param <taxon_in> The id of the leftmost taxon.
     @param <only_inpa> If set only calculates the number of inparalog pairs.
     @param <protein_start> The first protein index in the set.
     @param <protein_end> The last+1 protein index in the set.
     @returns the length of the proteins specified
  **/
  uint getProteinLength(uint taxon_in, bool only_inpa, uint protein_start, uint protein_end) {
    assert(taxon_in < (uint)taxon_length);
    assert(protein_end >= protein_start);
    if(list && list[taxon_in]) {
      uint sum = 0;
      for(uint i = protein_start; i < protein_end; i++) {      
	if(only_inpa) {
	  if(list[taxon_in][taxon_in]) sum += (uint)list[taxon_in][taxon_in]->getProteinLength(i);
	} else {
	  for(uint out = 0; out < (uint)taxon_length; out++) {
	    if(list[taxon_in][out]) sum += (uint)list[taxon_in][out]->getProteinLength(i);
	  }
	}
      }
      return sum;
    } else return 0;
  }  
  //! @returns the index start position of the protein specified
  uint getProteinIndexStart(uint in, uint out, int protein_in) {
    assert(list);
    assert(list[in]);
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    if(list[in][out]) return list[in][out]->get_index_start(protein_in);
    return 0;
  }  
  //! @returns the length of the file
  uint getLengthOfFile(uint in, uint out) {
    assert(list);
    assert(list[in]);
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    if(list[in][out])  return list[in][out]->getLengthOfFile();
    return 0;
  }
  //! @returns true if taxon-pair specified has data.
  bool buffer_is_not_null(uint taxon_in, uint taxon_out) {
    assert(list);
    assert(list[taxon_in]);
    assert(taxon_in < (uint)taxon_length);
    assert(taxon_out < (uint)taxon_length);
    if(list[taxon_in][taxon_out]) return list[taxon_in][taxon_out]->buffer_has_data_set();
    else return false;
  }

  //! Prints the keys and their lengths
  void print_listObjects() {
    printf("Prints the inded-list in list_file_parse.h:\n");
    if(list) {
      for(int in = 0; in< taxon_length; in++) {
	if(list[in]) {
	  for(int out = 0; out<taxon_length; out++) {
	    if(list[in][out]) {
	      if(list[in][out]->getTotalLengthOfData()) {
		printf("The key-list for [%d][%d] has (%d==%u) elements:\n",in,out, list[in][out]->getTotalLengthOfData(),  list[in][out]->getTotalLengthOfData_for_index());
		list[in][out]->print_list();
	      }
	    }
	  }
	}
      }
    }
  }
  //! Prints data about the class.
  void print_class_info() {
    FILE *f = stderr;
    fprintf(f, "--------------------------------------------------\n");
    fprintf(f, "Class Name: list_file_parse:\n");
    fprintf(f, "-\t Serves as a wrapper matrix for the squares of the data.\n");
    //  fprintf(f, "\t-\t taxon_length(%d) with taxons set to:\n", taxon_length); 
    taxa::print_class_info(f, listTaxa, taxon_length);  
  }

  //! Prints info about the elements size
  void print_info_size_matrix() {print_info_size_matrix(stdout);}
  //! Prints info about the elements size
  void print_info_size_matrix(FILE *f) {
    assert(list);
    fprintf(f, "taxon_length = %d\n", taxon_length);
    for(int in = 0; in< taxon_length; in++) {
      if(list[in]) {
	for(int out = 0; out<taxon_length; out++) {
	  if(list[in][out]) {
	    fprintf(f, "list[%d][%d] has %u elements reserved\n",in,out, list[in][out]->getLengthProteins());
	  }
	}
      } else fprintf(f, "data not set for row[%d]\n", in);
    }
  }
  //! Prints info about the pairs found in the matrix
  void print_info_pairs(FILE *f, const bool DEBUG_print_pairs_in_file_parse_log_file) {
    assert(list);
    assert(f);
    fprintf(f, "taxon_length = %d\n", taxon_length);
    for(int in = 0; in< taxon_length; in++) {
      if(list[in]) {
	for(int out = 0; out<taxon_length; out++) {
	  if(list[in][out]) {
	    list[in][out]->log_write_info(f, listTaxa, DEBUG_print_pairs_in_file_parse_log_file);//	  fprintf(f, "list[%d][%d] has %u protein-pairs.\n",in,out, list[in][out].get_total_pair_cnt());
	  }
	}
      } else fprintf(f, "data not set for row[%d]\n", in);
    }
    print_getTotal_memoryConsumption_for_object(f);
  }
  //! Prints data about the relations set in this object.
  void print_data(bool print_inner, bool print_outer) {
    // TODO: Think about alternatives nicely dumping data out.
    if(print_outer == print_inner && false) {;}
    fprintf(stderr, "!!\tFunction %s at line %d in class list_file_parse deactivated. Contact the developer if this method is of interest.\n", __FUNCTION__, __LINE__);
  }
/*   //! Prints data about the relations set in this object. */
/*   void print_data(uint taxon_in, uint taxon_out, bool print_inner, bool print_outer) { */
/*     if(print_outer == print_inner && false && taxon_out==taxon_in) {;} */
/*     // TODO Consider alternatives nicely dumping data out. */
/*     fprintf(stderr, "!!\tFunction %s at line %d in class list_file_parse deactivated. Contact the developer if this method is of interest.\n", __FUNCTION__, __LINE__); */
/*     //list[taxon_in][taxon_out].print(taxon_in, taxon_out, print_inner, print_outer, NULL, NULL, NULL, listTaxa); */
/*   } */
/*   //! Prints the data including the names */
/*   void print_names(bool print_inner, bool print_outer) { */
/*     // TODO: Think about alternatives nicely dumping data out. */
/*     if(print_outer == print_inner && false) {;} */
/*     fprintf(stderr, "!!\tFunction %s at line %d in class list_file_parse deactivated. Contact the developer if this method is of interest.\n", __FUNCTION__, __LINE__); */
/*   } */

//! Prints info about the relation:
  void print_relation_object_info(uint taxon_in, uint taxon_out) {
    assert(taxon_in < (uint)taxon_length);
    assert(taxon_out < (uint)taxon_length);
    if(list && list[taxon_in] && list[taxon_in][taxon_out]) {
      printf("# For taxon[%u][%u] the data set is (at line %d in file %s):\n", taxon_in, taxon_out, __LINE__, __FILE__);
      list[taxon_in][taxon_out]->print_relation_object_info(stdout);
    } else printf("(empty)\t For taxon[%u][%u] data not set; at line %d in file %s\n", taxon_in, taxon_out, __LINE__, __FILE__);
  }
  /**
     @brief Produces (writes) a log file holding the list of memory allocations for this object.
     @remarks 
     - Useful for verifying the parsing- and uptimizing the internal data structures.
     - In optimised mode does not produce the log file.
  **/
  void log_produce_memory_allocations(uint cpu_cnt, const bool DEBUG_print_pairs_in_file_parse_log_file) {
#ifndef NDEBUG
#ifdef LOG_WRITE_LIST_FILE_PARSE_CONTENT
    struct tms tmsstart; clock_t clock_log;
    if((clock_log = times(&tmsstart)) == -1) // starting values
      fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
    FILE *f_log = log_builder::get_log_file_pointer("list_file_parse", cpu_cnt, __FILE__, __LINE__);
    assert(f_log);

      fprintf(f_log, "For the %d taxa, the following data is set:\n", taxon_length);
      //      print_info_size_matrix(f_log);
      print_info_pairs(f_log, DEBUG_print_pairs_in_file_parse_log_file);

      // The size information of the object:
      const loint cnt_objects = taxon_length * taxon_length;
      const loint estimated_size = file_parse<T>::get_empty_object_size() *cnt_objects;
      const float size_gb = (float)estimated_size/(1024*1024*1024);
      fprintf(f_log, "\n- Detailes of memory consumption: (a) The file_parse<T>-object used as basis requuires each %lu B, (c) while a pointer to it requires %lu B, and (d) having in total %u taxa, for %lld objects, (e) each of size %lld B, (f) a total estimated size = %lld B =~ %0.4f GB\n", 
	      sizeof(file_parse<T>), sizeof(int),     taxon_length, cnt_objects, file_parse<T>::get_empty_object_size(), estimated_size, size_gb);

      // Generates the estimate on the time consumption:
      struct tms tmsend; clock_t end; long clktck = 0;
      if((end = times(&tmsend)) == -1) fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
      if((clktck =  sysconf(_SC_CLK_TCK)) < 0) fprintf(stderr, "!!\tAn error in sys-configuration at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
      const clock_t t_real = end - clock_log;
      const double log_time = t_real/(double)clktck;
      fprintf(f_log, "--\tUsed in total %10.8f seconds for this operation\n", log_time);
      fclose(f_log);
#else
    //    printf("not defined in list_file_parse.h\n"); 
#endif
#endif
  }
  //! Prints the data including the names with regard to the inparalogs
  void print_inpa(bool use_names) { print_data(true, use_names);}
  //! Simplified method for printing the data including the names of this object.
  void print_all(bool use_names) { print_data(false, use_names);  }
  //! @returns the list (for backward compability)
  file_parse<T> **getData(){return list;}
  //! @return true if the row at index argument given, is not NULL
  bool hasRow(uint in) {
    if(list!=NULL) {
      if(in < (uint)taxon_length) {
	if (list[in]!=NULL) return true;
      }
    }
    return false;
  }

  //! Compare the two classes:
  bool is_equal(list_file_parse *obj, const bool print_info) {
    bool is_equal = true;
    if(list) {
      if(!obj) {
	if(print_info) printf("!!\tArgument object of type list_file_parse not set\n");
	is_equal = false;
      } else {
	for(uint i = 0; i < (uint)taxon_length; i++) {
	  if(has_data(i)) {
	    if(!(obj->has_data(i))) {
	      if(print_info) printf("!!\tArgument object of type *file_parse, found at index[%u], not set\n", i);
	      is_equal = false;
	    } else {
	      for(uint out = 0; out < (uint)taxon_length; out++) {
		if(has_data(i)) {
		  if(!(obj->has_data(i, out))) {
		    if(print_info) printf("!!\tArgument object of type file_parse, found at index[%u][%u], not set\n", i, out);
		    is_equal = false;
		  } else {
		    if(!(list[i][out]->is_equal(*obj->get_element(i, out), print_info))) {
		      if(print_info) printf("!!\tArgument object of type file_parse, found at index[%u][%u], is not equal this\n", i, out);
		      is_equal = false;
		    }
		  }
		} else if(obj->has_data(i, out)) {
		  if(print_info) printf("!!\t The internal object (this) of type file_parse, found at index[%u][%u], not set\n", i, out);
		  is_equal = false;
		}
	      }
	    }      
	  } else if(obj->has_data(i)) {
	    if(print_info) printf("!!\tThe internal object (this) of type *file_parse, found at index[%u], not set\n", i);
	      is_equal = false;
	  }
	}
      }
    } else if(obj) {
      if(print_info) printf("!!\tThe internal object (this) of type list_file_parse not set\n");
      is_equal = false;
    }
    return is_equal;
  }
    

  //! Reduces the length of one index: Coresponds to discarding one element
  void decreaseLengthOfListElement(uint in, uint out, int protein_index) {
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    assert(list);
    assert(list[in]);
    if(hasRow(in) && list[in][out])  list[in][out]->decreaseLengthOfListElement(protein_index);
  }
  /**
     @brief Merges the argument with 'this'
     @param <data> The list_file_parse object to merge
     @remarks Deallocates the data given
  **/
  void merge(list_file_parse *&data) {
    protein_relation prot = protein_relation();
    uint cnt_elements_in_all_relation_lists = 0;
    merge_data(data, prot, cnt_elements_in_all_relation_lists);
  }
  
  /**
     @brief Merges the argument with 'this'
     @param <data> The list_file_parse object to merge
     @param <rel> Identifies poteintial merging zone.
     @remarks Deallocates the data given
  **/
  void merge_data(list_file_parse *&data, struct protein_relation rel, uint &cnt_elements_in_all_relation_lists) {
    assert(data);
    if(list != NULL) {
#ifndef NDEBUG
      const uint tot_elements = getTotalLengthOfData() + data->getTotalLengthOfData();
#endif
      for(int in = 0; in < taxon_length; in++) {
	if(data->hasRow(in)) {
	  assert(list[in]);	      
	  for(int out = 0; out<taxon_length; out++) {
	    if(data->has_data(in, out)) {
	      if(list[in][out]) {
		if(list[in][out]->has_data()) { 
		  if(rel.is_set() && list[in][out]->gives_overlap(in, out, rel.get_protein_in(), rel.get_protein_out())) {
		    // Adsjust the position by setting offset = '1':
		    data->decreaseLengthOfListElement(in, out, rel.protein_in); 	    
		    //		  list[in][out].merge_buffers(data->getListElement(in, out)); 
		  } //else
		  file_parse<T> *element = data->getListElementRef(in, out);
		  list[in][out]->merge_buffers(element, listTaxa[in], cnt_elements_in_all_relation_lists, &listTaxa[out]);
		} else {
		  list[in][out]->finalize(false); // Remove any memory allocations first.
		  data->steal_data_inverse(in, out, list[in][out]);
		  assert(!data->getListElement(in, out).has_data());
		}
	      } else { //no data is given for this, so grabs the data as-is
		data->steal_object_inverse(in, out, list[in][out]);   // Gets the object
	      }
	    }
	  }
	}
      }
#ifndef NDEBUG
      assert(tot_elements == getTotalLengthOfData()); // The 'root' should now consist of the sum of elements.
#endif
    } else fprintf(stderr, "!!\t\tinternal list is not set (list==NULL) at line %d in fucntion %s found in file %s. Please contact the  developer, oekseth@gmail.com.\n", __LINE__, __FUNCTION__, __FILE__);
    list_file_parse<T>::close(data, false); 
  }

  /**
     @brief Merges the argument with 'this'
     @param <data> The file_parse object to merge
     @param <taxon_in> The id of the innermost taxon.
     @param <taxon_out> The id of the outmost taxon.
     @remarks Deallocates the data given
  **/
  void merge_data(file_parse<T> *&data, uint taxon_in, uint taxon_out, uint &cnt_elements_in_all_relation_lists) {
    if(data == NULL) return; // No data set for the argument.
    if(list && data) {
#ifndef NDEBUG
      const uint this_length_old = getTotalLengthOfData();
      const uint data_length = data->getTotalLengthOfData();
      const uint tot_elements = this_length_old + data_length;
      const float data_sum_of_distances = data->get_sum_of_pair_distances_for_myrank_in_memory();
      const float sum_before_merging = get_sum_of_pair_distances_in_memory();
#endif
      assert(taxon_in < (uint)taxon_length);
      assert(taxon_out < (uint)taxon_length);
      assert(list[taxon_in]);
      if(!list[taxon_in][taxon_out]) {	
	list[taxon_in][taxon_out] = data;
	data = NULL;
	//	new file_parse<p_rel>((int)taxon_in, (int)taxon_out, (int)data->get_index_length(), 0, MAX_BUFFER_SIZE, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE); // UINT_MAX to avoid writing of data to the file
	//	list[taxon_in][taxon_out] = new file_parse<p_rel>(taxon_in, taxon_out, data->get_index_length(), 0, data->get_buffer_in_mem_pos(),  MAX_BUFFER_SIZE, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE); // UINT_MAX to avoid writing of data to the file
	if(!list[taxon_in][taxon_out]) { // If memory was not allocated:
	  fprintf(stderr, "!!\tMemory was not allocated for taxon-pair[%u, %u]. Aborts. Please refer to the manual or contact the developer at oekset@gmail.com. (This error is found in %s at line %d in method %s\n", taxon_in, taxon_out, __FILE__, __LINE__, __FUNCTION__);
	  exit(2);
	}
	assert(list[taxon_in][taxon_out]);
      } else {
	if(list[taxon_in][taxon_out]) {
	  if(list[taxon_in][taxon_out]->has_data()) { 		  
	    list[taxon_in][taxon_out]->merge_buffers(data, listTaxa[taxon_in],cnt_elements_in_all_relation_lists, &listTaxa[taxon_out]);
	  } else {
	    list[taxon_in][taxon_out]->finalize(false); // Remove any memory allocations first.
	    list[taxon_in][taxon_out] = data;
	    assert(data->buffer_has_data_set()); // Ensures that data is both set, and consistent (both list and size is set).
	    data = NULL; // Removes the pointer for the param given.
	    assert(buffer_is_not_null(taxon_in, taxon_out)); // Ensures that data is both set, and consistent (both list and size is set).
	  }
	}
      }
#ifndef NDEBUG
      //! Validate that the total number of elements before corresponds to the sum after:
      assert(tot_elements == getTotalLengthOfData()); // The 'root' should now consist of the sum of elements.

      //! Validates that the structure is correctly updated:
      const float expected_sum_of_distances_after_merge = get_sum_of_pair_distances_in_memory();
      const float inserted_sum = sum_before_merging + data_sum_of_distances;
      log_builder::compare_floats(expected_sum_of_distances_after_merge, inserted_sum, tot_elements, __LINE__, __FILE__, __FUNCTION__);
#endif
    } else if(data){
      fprintf(stderr, "!!\t\tinternal list is not set (list==NULL) at line %d in fucntion %s found in file %s. Please contact the  developer, oekseth@gmail.com.\n", __LINE__, __FUNCTION__, __FILE__);
    }
    //! Deallocates the object given:
    if(data) {file_parse<T>::close(data, false);}
  }


  //! Sets the given object at the given position
  void set_object(file_parse<T> *obj, uint in, uint out) {
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    assert(list);
    assert(list[in]);
    assert(!list[in][out]);
    list[in][out] = obj;
  }
  //! A wrapper for stealing data from the argument list.
  void steal_data(file_parse<T> *p, uint in, uint out) {
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    assert(list);
    assert(list[in]);
    if(list[in][out]) list[in][out]->steal_data(p);
  }
  //! A wrapper for stealing data from the argument list.
  void steal_data_inverse(uint in, uint out, file_parse<T> *p) {
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    assert(list);
    assert(list[in]);
    if(list[in][out]) list[in][out]->steal_data_inverse(p);
  }

  //! A wrapper for stealing a file_parse object from the argument list.
  void steal_object_inverse(uint in, uint out, file_parse<T> *&p) {
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    assert(list);
    assert(list[in]);
    p = list[in][out];
    list[in][out] = NULL; // Sets it to empty.
  }


  //! @returns the class at the index specified
  file_parse<T> *getListElementRef(uint in, uint out) {
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    assert(list);
    assert(list[in]);
    return list[in][out];
  }
  //! @returns the class at the index specified
  file_parse<T> getListElement(uint in, uint out) {
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    assert(list);
    assert(list[in]);
    return list[in][out];
  }
  /**
     A wrapper for stealing data from the argument class.
     @param <arg> The list_file_parse object to 'steal'.
     @remarks  
     - As taxa list has not changed at this stage in the processing, it is not included in the merging.
     - Sets the the lists residing in input to empy as those are taken.
  */
  void steal_data(list_file_parse arg) {
    MAX_BUFFER_SIZE = arg.MAX_BUFFER_SIZE;
    if(!listTaxa) listTaxa = arg.get_listTaxa(taxon_length);
    if(!FILE_BINARY_LOCATION) FILE_BINARY_LOCATION = arg.get_FILE_BINARY_LOCATION();
    for(int i = 0; i  < taxon_length; i++) {
      if(arg.hasRow(i)) {
	for(int out = 0; out<taxon_length; out++) {
	  if(list[i] && list[i][out]) {
	    list[i][out]->steal_data(arg.getListElementRef(i, out)); // Grabs the data, and intialises the values to the argument into zero (empy) values.
	  }
	}
      }
    }
  }
  /**
     @brief Purpose of the alogrithm below is to send data who is above the threshold size set, to the class 'pipe_merge.xx'
     @remarks The sending of data to the pipe builds on empirical evidence..
  **/
  list_file_parse *createPacket() {
    list_file_parse *data = new list_file_parse<T>(hashProtein, taxon_length, listTaxa, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE); 
    // Data for the optimal size:
    // - size=1000000: 11s for p, 4s for the r(emaining part)
    // - size=500000:  12s for p, 3s for the r(emaining part)
    // - size=300000:  13s for p, 2s for the r(emaining part)
    // - size=200000:  18s for1 p, 2s for the r(emaining part)
    // - size=200000:  14s for1 p, 2s for the r(emaining part)
    // - size=100000:  17s for p, 1s for the r(emaining part)
    
    //   const uint max_mem_size = UINT_MAX; 
    //    const uint max_mem_size = 100000; // TODO: Evaluate if this gives the best performance
    //    const uint max_mem_size = 10000; // TODO: Evaluate if this gives the best performance
    //    const uint max_mem_size = 1000; // TODO: Evaluate if this gives the best performance
    const uint max_mem_size = 0; // TODO: Evaluate if this gives the best performance
    bool data_added = false;
    for(int i = 0; i < taxon_length; i++) {
      //    for(int i = taxon_length-1; i >-1; i--) {
      for(int out = 0; i < taxon_length; i++) {
      //      for(int out = taxon_length-1; out >-1; out--) {
	assert(list[i]);
	if(list[i][out]) {
	  if(list[i][out]->getLengthBuffer() > max_mem_size) {
	    //	    const int length_prots = list[i][out]->getLengthProteins();
	    data->set_object(list[i][out], i, out); 
	    list[i][out] = NULL; // Ampties this.
	    //	    data->steal_data(list[i][out], i, out); 
	    // Then reallocates memory for the area
	    //file_parse<T>(i, out, length_prots, 0,  UINT_MAX, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE); // UINT_MAX to avoid writing of data to the file
	    data_added = true;
	  } //else fprintf(stderr, "in createPacket at list_file_parse; did not reserve data for id [%u][%u]\n", i, out);
	}
      }
    }
    if(true) return data;
    else { // TODO: make the below option (returning an NULL element without causing a crash later on) possible.
      if(data_added) return data;
      else {
	list_file_parse::close(data, false);
	return NULL;
      }
    }
  }

  /**
     @brief Produces the inparalogs, inserting them in the given container
     @return the number of inparalogs inserted
  ***/
  uint produceInparalogs(mcl_format_t *container, char *protein_name, uint taxon, uint world_index, uint local_index, float avg_div_factor) {
    if(has_data(taxon, taxon, local_index)) {
      return list[taxon][taxon]->produceInparalogs(listTaxa, MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, container, protein_name, taxon, world_index, local_index, /*listTaxa[taxon].rel_start, */avg_div_factor);
    } else return 0;
  }
  /**
     Inserts a new relation into the list
     @return true if protein was inserted as a new protein, and not updating an exsisting one.
     @remarks Returnvalue to be used in conjunction with the method named 'getTotalLengthOfData()' found in class list_file_parse.
  **/
  bool insert_new_rel(uint taxon_in, uint taxon_out, Parse p, protein_relation prot);
  /**! Inserts a relation into the structure.  */  
  void insert_relation(uint taxon_in, uint taxon_out, uint protein_in, uint protein_out, float sim_avg) {

    assert(list);
    assert(taxon_in < (uint)taxon_length);
    assert(taxon_out < (uint)taxon_length);
    assert(list[taxon_in]);
    if(list && list[taxon_in] && (taxon_out < (uint)taxon_length) && (taxon_in < (uint)taxon_length)) {
      if(list[taxon_in] != NULL){ 
	if(!list[taxon_in][taxon_out]) {
	  loint length = 0;
	  if(hashProtein) length = hashProtein[taxon_in].getLength();
	  list[taxon_in][taxon_out] = new file_parse<T>(taxon_in, taxon_out, length, 0,  MAX_BUFFER_SIZE, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE); // UINT_MAX to avoid writing of data to the file
	  if(!list[taxon_in][taxon_out]) { // If memory was not allocated:
	    fprintf(stderr, "!!\tMemory was not allocated for taxon-pair[%u, %u]. Aborts. Please refer to the manual or contact the developer at oekset@gmail.com. (This error is found in %s at line %d in method %s\n", taxon_in, taxon_out, __FILE__, __LINE__, __FUNCTION__);
	    exit(2);
	  }
	}
	list[taxon_in][taxon_out]->insert_new_rel(protein_in, protein_out, sim_avg, 0, 0);
      }
    }
  }




   /**
     @brief Prints the memory consumption for the object.
     @remarks For anaylsing the memory fingerprint.
  **/
  void print_getTotal_memoryConsumption_for_object(FILE *f) {
    assert(f);
    loint size_reserved = 0; 
    loint size_used = 0; 
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in]) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      size_reserved += list[taxon_in][taxon_out]->getTotal_memoryConsumption_for_object();
	      size_used     += list[taxon_in][taxon_out]->getTotal_used_memoryConsumption_for_object();
	    }
	  }
	}
      }         
    } 
    const float fill_factor = (float)size_used / size_reserved;
    const float size_gb = (float)size_reserved/(1024*1024*1024);
    char end_char = ' ';
    // Prints the data, but in order to avoid clutter, newline is avoided if the goal is printing the data to the log file:
    if((f==stdout) || (f == stderr)) end_char = '\n';
    fprintf(f, "- Has %lld B of memory =~ %.5f GB => fill_factor(%.4f)%c", size_reserved, size_gb, fill_factor, end_char);
  }

  //! Writes the buffer to stdout, using integres as representation.
  void print_buffer_as_integer_file(uint taxon_in, uint taxon_out) {
    if(list && taxon_in < (uint)(taxon_length) && taxon_out < (uint)taxon_length && list[taxon_in] && list[taxon_in][taxon_out]) {
      list[taxon_in][taxon_out]->print_buffer_as_integer_file(taxon_in, taxon_out);
    }
  }

  //! Writes the buffer to the pointer, using integres as representation.
  void write_buffer_as_integer_file(FILE *f) {
    assert(f);
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in]) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      if(list && taxon_in < (uint)(taxon_length) && taxon_out < (uint)taxon_length && list[taxon_in] && list[taxon_in][taxon_out]) {
		list[taxon_in][taxon_out]->print_buffer_as_integer_file(taxon_in, taxon_out, f);
	      }
	    }
	  }
	}
      }
    }
  }

  /**
     @brief Sums the memory consumption for the object.
     @remarks For anaylsing the memory fingerprint.
     @return the total length of data in Bytes.
  **/
  loint getTotal_memoryConsumption_for_object() {
    loint length = sizeof(list_file_parse<T>);
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in]) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      length += list[taxon_in][taxon_out]->getTotal_memoryConsumption_for_object();
	    }
	  }
	}
      }         
    } 
    return length;
  }

  /**
     @brief Sums the of references in the listIndex.
     @remarks For verifying the copying (merging) process.
     @return the total length of data.
  **/
  uint getTotalLengthOfData_for_index() {
    uint length = 0;
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in]) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      length += list[taxon_in][taxon_out]->getTotalLengthOfData_for_index();
	    }
	  }
	}
      }         
    } 
    return length;
  }  

  /**
     @brief Sums the of references in the listIndex.
     @remarks For verifying the copying (merging) process.
     @return the total number of index_t objects.
  **/
  uint get_index_length() {
    uint length = 0;
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in]) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      length += list[taxon_in][taxon_out]->get_index_length();
	    }
	  }
	}
      }         
    } 
    return length;
  }  

  /**
     @brief Sums the relations.
     @remarks For verifying the copying (merging) process.
     @return the total length of data.
  **/
  uint getTotalLengthOfData() {
    uint length = 0;
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in]) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      length += list[taxon_in][taxon_out]->getTotalLengthOfData();
	    }
	  }
	}
      }         
    } 
    return length;
  }

  /**
     @brief Sums the relations for the orthologs
     @remarks For verifying the copying (merging) process.
     @return the total length of data.
  **/
  uint getTotalLengthOfData_for_orthologs() {
    uint length = 0;
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in]) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      if(taxon_in != taxon_out) {
		length += list[taxon_in][taxon_out]->getTotalLengthOfData();
	      }
	    }
	  }
	}
      }         
    } 
    return length;
  }

  /**
     @brief Sums the relations for the inparalogs
     @remarks For verifying the copying (merging) process.
     @return the total length of data.
  **/
  uint getTotalLengthOfData_for_inparalogs() {
    uint length = 0;
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in][taxon_in]) {
	  length += list[taxon_in][taxon_in]->getTotalLengthOfData();
	}
      }         
    } 
    return length;
  }


  /**
     @brief Sums the relations, given the relation is part of myranks work.
     @remarks For verifying the copying (merging) process.
     @return the total length of data.
  **/
  uint getTotalLengthOfData_for_myrank() {
    uint length = 0;
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in] && taxon_is_myranks_responsibility(taxon_in)) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      length += list[taxon_in][taxon_out]->getTotalLengthOfData();
	    }
	  }
	}
      }         
    } 
    return length;
  }


  /**
     @brief Sums the relations, given the relation is part of myranks work.
     @remarks For verifying the copying (merging) process.
     @return the total length of data.
  **/
  uint getTotalLengthOfData_for_inparalogs_myrank(uint taxon_start, uint taxon_length) {
    uint length = 0;
    if(list) {
      for(uint taxon_in = taxon_start; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in] && taxon_is_myranks_responsibility(taxon_in)) {
	  if(list[taxon_in][taxon_in]) {
	    length += list[taxon_in][taxon_in]->getTotalLengthOfData();
	  }
	}
      }         
    } 
    return length;
  }

  /**
     @brief Sums the relations, given the relation is part of myranks work.
     @remarks For verifying the copying (merging) process.
     @return the total length of data.
     @remarks If macro variable 'USE_MPI' is defined, only return the sum of values
     for those considered relevant.
  **/
  float get_sum_of_pair_distances_for_myrank_in_memory() {
    float length = 0;
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in] && taxon_is_myranks_responsibility(taxon_in)) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      /*
#ifdef USE_MPI
	      if(list_of_nodes_taxa_responsilibties) {
		if(list_of_nodes_taxa_responsilibties[taxon_in] || list_of_nodes_taxa_responsilibties[taxon_out]) {
		  length += list[taxon_in][taxon_out]->get_sum_of_pair_distances_for_myrank_in_memory();
		}
	      } else length += list[taxon_in][taxon_out]->get_sum_of_pair_distances_for_myrank_in_memory();
#else
*/
	      length += list[taxon_in][taxon_out]->get_sum_of_pair_distances_for_myrank_in_memory();
	      //#endif
	    }
	  }
	}
      }         
    } 
    return length;
  }

  /**
     @brief Sums the relations, given the relation is part of myranks work.
     @remarks For verifying the copying (merging) process.
     @return the total length of data.
     @remarks If macro variable 'USE_MPI' is defined, only return the sum of values
     for those considered relevant.
  **/
  float get_sum_of_pair_distances_in_memory() {
    float length = 0;
    if(list) {
      for(uint taxon_in = 0; taxon_in < (uint)taxon_length; taxon_in++) {
	if(list[taxon_in]) {
	  for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	    if(list[taxon_in][taxon_out]) {
	      length += list[taxon_in][taxon_out]->get_sum_of_pair_distances_for_myrank_in_memory();
	    }
	  }
	}
      }         
    } 
    return length;
  }

  /**
     @brief Sums the relations, given the relation is part of myranks work.
     @remarks For verifying the copying (merging) process.
     @return the total length of data.
  **/
  float get_sum_of_pair_distances(uint taxon_in, uint taxon_out) {
    float length = 0;
    if(list) {
      if(list[taxon_in]) { // && taxon_is_myranks_responsibility(taxon_in)) {
	if(list[taxon_in][taxon_out]) {
	  length += list[taxon_in][taxon_out]->get_sum_of_pair_distances_for_myrank_in_memory();
	}
      }         
    } 
    return length;
  }

  /**
     @brief  Sums the relations
     @param <taxon_in> The taxon id to sum the columns for.
     @return The total number of pairs for a given inner taxon
  **/
  uint get_total_number_of_pairs_for_taxon(uint taxon_in) {
    assert(taxon_in < (uint)taxon_length);
    uint length = 0;
    if((taxon_in < (uint)taxon_length) && list) {      
      if(list[taxon_in]) {
	for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	  if(list[taxon_in][taxon_out]) {
	    length += list[taxon_in][taxon_out]->getTotalLengthOfData();
	  }
	}
      }
    } 
    return length;
  }
  /**
     @brief Sums the relations.
     @param <taxon_start> The id of the first taxa to build the list for.
     @param <taxon_end> The id of the first taxa (and the following onwards) to not include in the list build.
     @param <is_inparalog> If set, only calculates the size of its inparalogs.
     @param <biggest_collection_size> The Biggest collection, used as basis setting the worst-case list size in the building of sheduling set.
     @remarks Useful when building lists for thread sheduling.
     @return the total length of data for the given range
  **/
  uint getTotalLengthOfData(uint taxon_start, uint taxon_end, const bool is_inparalog, uint &biggest_collection_size) {
    assert((int)taxon_start < taxon_length);
    assert((int)taxon_end <= taxon_length);
    assert(list);
    assert(list[taxon_start]);
    biggest_collection_size = 0;
    uint length = 0;
    for(uint taxon_in = taxon_start; taxon_in < (uint)taxon_length; taxon_in++) {
      uint temp_length = 0;
      if(is_inparalog) {if(list[taxon_in][taxon_in]) temp_length += list[taxon_in][taxon_in]->getTotalLengthOfData();}
      else {
	for(uint taxon_out = 0; taxon_out < (uint)taxon_length; taxon_out++) {
	  if(taxon_in != taxon_out) {if(list[taxon_in][taxon_out]) temp_length += list[taxon_in][taxon_out]->getTotalLengthOfData();}
	}
      }
      length += temp_length;
      biggest_collection_size = max(temp_length, biggest_collection_size);
    }
    return length;
  }
  /**! @return the maximal optimal size of the buffer to use */
  mem_loc getMaxBufferSize() {
    mem_loc size = 0;
    for(uint taxon_id = 0; taxon_id < (uint)taxon_length; taxon_id++) {
      assert(list[taxon_id]);
      if(list[taxon_id][taxon_id]) {
	const uint this_size = 10+list[taxon_id][taxon_id]->getTotalLengthOfData();
	if(this_size > size) {
	  loint temp_pos, temp_size;
	  T *temp = T::init_list(this_size*taxon_length, USE_BEST_BLAST_PAIR_SCORE, temp_pos); 
	  if(temp != NULL) {
	    size=this_size;
	    T::close(temp, temp_pos, temp_size);
	  }	
	}
      }
    }
    return size;
  }

  /**! Returns the size of data for the given protein  */
  uint getSizeofBuffer(uint taxon_in, uint taxon_out, uint protein_start, uint protein_end)  {
    assert(taxon_in < (uint)taxon_length);
    assert(taxon_out < (uint)taxon_length);
    if(has_data(taxon_in, taxon_out))
      return list[taxon_in][taxon_out]->getTotalRelations(protein_start, protein_end);
    return 0;
  }

  //! Gets the buffer of type file_parse<T> 
  file_parse<T> ***getBuffer() {
    return list;
  }

//! Gets the buffer of type file_parse<T> 
  file_parse<T> ***get_buffer() {
    return list;
  }

  //! Gets the buffer specified
  rel_t *getBuffer(uint in, uint out, uint protein_in) {
    assert(in < (uint)taxon_length);
    assert(out < (uint)taxon_length);
    if(has_data(in, out, protein_in)) {
      return list[in][out]->getBuffer(protein_in);
    }
    return NULL;
  }

  //! Returns the number of elements in the buffer, given the arguments
  uint getLength(uint taxon_one, uint taxon_two, uint protein_in)  {
    // Below asserts included in order to throw an error if it was a "non-deliberate" use of the taxon-length-sizes.
    assert(taxon_one < (uint)taxon_length);
    assert(taxon_two < (uint)taxon_length);
    if(list && list[taxon_one] && list[taxon_one][taxon_two] && (taxon_one < (uint)taxon_length) && (taxon_two < (uint)taxon_length)) {
      return list[taxon_one][taxon_two]->getProteinLength(protein_in);
    }
    return 0;
  }
  /**! Gets the buffer: if not already opendes, if retrives it from the a file*/
  basket_parse *getBufferFromFile(uint taxon_one, uint taxon_two) {
    if(list && (taxon_one< (uint)taxon_length) && (taxon_two< (uint)taxon_length) && list[taxon_one] && list[taxon_one][taxon_two] && (taxon_one < (uint)taxon_length) && (taxon_two < (uint)taxon_length)) 
      return list[taxon_one][taxon_two]->getBuffer();
    else return NULL;
  }
  /**! Makes the data consitent, and accessible  */
  void dumpBufferOrFile(const bool is_inparalog, uint taxon_start, uint taxon_end) {
    dump_buffer_if_above_limit(UINT_MAX, is_inparalog, taxon_start, taxon_end);
  }
  /**! Makes the data consitent, and accessible  */
  void getBufferIntoMemory(const bool is_inparalog, uint taxon_start, uint taxon_end) {
    dump_buffer_if_above_limit(UINT_MAX, is_inparalog, taxon_start, taxon_end, false);
  }
  //! Gets all the relations into the buffer, given the space defined by input parameters.
  void include_all_relations(uint taxon_start, uint taxon_end) {
    dump_buffer_if_above_limit(UINT_MAX, false, taxon_start, taxon_end, true);
  }
  /**
     @brief Write the buffer to memory if it's above threshold
     @param <max_buffer_size> The threshold size the buffer must be above in order to be written to memory.
     @param <only_inparalogs> If set only the inparalogs are written, else the others.
     @param <taxon_start> The first of the inner pairs to process
     @param <taxon_end> The last+1 of the inner pairs to process, ie, "taxon_length" is to be used if the remaining pairs are to be written.
     @param <get_all_relations> Like setting "only_inparalogs" parameter to boht true and false.
     @remarks
     # Called at before operation(s), like the inpa-inpa operation.
  **/
  void dump_buffer_if_above_limit(const uint max_buffer_size, const bool only_inparalogs, uint taxon_start, uint taxon_end, bool get_all_relations) {
    // Variables to be set correctly before call
    assert(taxon_start < (uint)taxon_length);
    assert(taxon_end <= (uint)taxon_length);
    // Criterias to be tested at run-time:
    if(list && list[taxon_start]) {
      if(only_inparalogs) {
	for(uint i = taxon_start; i < (uint)taxon_end; i++) {
	  if(list[i] && list[i][i]) list[i][i]->dumpBufferOrFile(max_buffer_size); // if possible, gets the data into the buffer for the whole inparalog
	}
      } else {
	for(uint i = taxon_start; i < taxon_end; i++) {
	//	for(uint i = 0; i < (uint)taxon_length; i++) {
	  if(list[i]) {
	    for(uint out = taxon_start; out < (uint)taxon_end; out++) {
	    //	    for(uint out = 0; out < (uint)taxon_length; out++) {
	      if(get_all_relations || (i != out))  {if(list[i][out]) list[i][out]->dumpBufferOrFile(max_buffer_size);}
	    }
	  }
	}
      }
    }
  }


  
  /**
     @brief Opens a buffer and dumps it to memory if the space  for it is allocated
     @remarks
     # Called at before operation(s), like the inpa-inpa operation.
  **/
  void getFromBufferToStruct(uint max_buffer_size,const bool only_inparalogs)  {
    dump_buffer_if_above_limit(max_buffer_size, only_inparalogs, 0, taxon_length, false);
/*     for(uint i = 0; i < (uint)taxon_length; i++) { */
/*       if(!only_inparalogs) { */
/* 	if(list[i] != NULL) { */
/* 	  for(uint out = 0; out < (uint)taxon_length; out++) { */
/* 	    list[i][out].dumpBufferOrFile(max_buffer_size); */
/* 	  } */
/* 	} else fprintf(stderr, "!!\tStructure not filled in list_file_struct_depricated at line %d in function %s\n",__LINE__, __FUNCTION__); */
/*       } else list[i][i].dumpBufferOrFile(max_buffer_size); // only the inparalogs */
/*     } */
  }
  /**! Dump the either the inparalogs- or the 'other relations' to the file specified, removing the data from the list.*/
  void dump_to_file(const bool is_inparalog, uint taxon_start, uint taxon_end) {
    // Variables to be set correctly before call
    assert(taxon_start < (uint)taxon_length);
    assert(taxon_end <= (uint)taxon_length);
    // Criterias to be tested at run-time:
    if(list && list[taxon_start]) {
      if(is_inparalog) {
	for(uint i = taxon_start; i < (uint)taxon_end; i++) {
	  if(list[i] && list[i][i]) list[i][i]->dump_buffer_to_file();
	}
      } else {
	for(uint i = taxon_start; i < (uint)taxon_end; i++) {
	  if(list[i]) {
	    for(uint out = taxon_start; out < (uint)taxon_end; out++) {
	      if(i != out) {if(list[i][out]) list[i][out]->dump_buffer_to_file();}
	    }
	  }
	}
      }
    }
  }
  /**
     @brief Writes the buffers to a file if the amount of data is above the assumed memory whom is required duirng the operation
     @param <use_minimum_amount_of_memory> is set to true, does not use the optimal amount of memory iot enhance speedup
  */
  void write_buffers_to_file(bool use_minimum_amount_of_memory) {
    mem_loc max_buffer_size = 0; // TODO: This value is unsed; should it be used?
    if(!use_minimum_amount_of_memory) max_buffer_size = getMaxBufferSize();       // Dump the buffer, if it is to large, to a file: else the data from the file into the buffer
    for(int taxon_in = 0;taxon_in<taxon_length;taxon_in++) {
      if(list[taxon_in]) {
	for(int taxon_out = 0;taxon_out<taxon_length;taxon_out++) {
	  if(list[taxon_in][taxon_out]) list[taxon_in][taxon_out]->dumpBufferOrFile(0); // dumps everything
	}
      }
    } 
  }
  /**
     @brief Sums the number of pairs in the collection.
     @remarks Some pipes need a pool of proteins. The task of this function is returning the last
     protein at the end of each buffer, implying that it would be the start of the next bufer as well.
     @param <start_pos> The protein id to start the count of pairs from.
     @param <taxon_in>  The inner taxon to count the pairs from.
     @param <only_inpa> If set to true, do not count the pairs from other taxa than the given taxon_in value. Else iterates over all the taxa in collection.
     @param <max_cnt_proteins> Included to ensure that all threads get approx equal amount of data.
     @return the protein counts.
     @author Ole Kristian Ekseth.
  */
  uint getProteinStartOfNextBuffer(uint start_pos, uint taxon_in, const bool only_inpa, uint max_cnt_proteins) {
    assert((int)taxon_in < taxon_length);
    uint cnt_prots = 0;
    if(list && list[taxon_in]) {
      if(listTaxa) {
	for(uint prot_id = start_pos; prot_id < (uint)listTaxa[taxon_in].total_cnt; prot_id++) {
	  if(cnt_prots >= max_cnt_proteins) {
	    return prot_id;
	  }
	  else {	
	    if(only_inpa) {
	      if(list[taxon_in][taxon_in]) cnt_prots += list[taxon_in][taxon_in]->getProteinLength(prot_id);
	    } else { // calculates for the outer relations
	      for(uint taxon_out = 0; taxon_out <  (uint)taxon_length; taxon_out++) {
		if(taxon_out != taxon_in) {
		  if(list[taxon_in][taxon_out] && list[taxon_in][taxon_out]->has_data(prot_id)) {
		    cnt_prots += list[taxon_in][taxon_out]->getProteinLength(prot_id);
		  }
		}
	      }
	    }
	  }
	}
      } else {
	fprintf(stderr, "!!\tglobal list listTaxa in class 'list_file_parse' not set at line %d in file %s. Contact oekseth@gmail.com.\n", __LINE__, __FILE__);
	assert(false);
      }
    } else {
      fprintf(stderr, "!!\trelations for  taxon %u in class 'list_file_parse' not set at line %d in file %s. Contact oekseth@gmail.com.\n", taxon_in, __LINE__, __FILE__);
      assert(false);
    }
    return listTaxa[taxon_in].total_cnt -1; // if at this point, the bucket was not completely filled with data
  }

  //! Deallocates the hash. @remarks Ends the life of the hash labeled proteins.
  void free_hashProteins() {
    if(hashProtein) {
      for(uint i =0; i< (uint)taxon_length;i++)
	hashProtein[i].delete_buffer();
      free(hashProtein); hashProtein=NULL;
    }
  }
 private:
  //! @return the available memory in Bytes on the system kernel.
  static loint get_available_memory_on_system() {
    loint amount_of_free_memory_b = 0;
    const char *cmd = "free -lb | grep \"Mem\""; 
    FILE* pipe = popen(cmd, "r");
    if(pipe) {
      char buffer[128];
      log_builder::test_memory_condition_and_if_not_abort(buffer!=NULL, __LINE__, __FILE__, __FUNCTION__);
      memset(buffer, '\0', 128);
       while(!feof(pipe)) { 
 	if(fgets(buffer, 128, pipe) != NULL) { 
	  printf("string gave %s\n", buffer);
	  
	  // Mem:           125         42         83          0          0          0 
 	  if(!sscanf(buffer, "%*s %*d %*d %lld   %*d %*d  %*d", &amount_of_free_memory_b)) { 
	    log_builder::throw_warning(software_dependencies, __LINE__, __FILE__, __FUNCTION__, "Did not find a value in the fird column of numbers using the tool 'free -b| grep \"Mem\"' on your system, therefore unable to predict the amount of free memory available. As a susbstitution, it's assumed that 16GB of data is avialable.");
	    //! Note: Below updated after tip from Mark DeBeen, which received a message named "integer overflow in expression" when setting "amount_of_free_memory_b = 16*1024*1024*2014;"
	    const loint max_on_system = std::numeric_limits<int>::max();
	    assert(max_on_system);
	    if((max_on_system / 1024) > (16*1024*1024)) {	      
	      amount_of_free_memory_b = 16*1024*1024*2014;
	    } else {amount_of_free_memory_b = max_on_system;}
 	  }  
	} 
       }
       //      delete buffer, buffer = NULL;
      pclose(pipe);
    } else {
      log_builder::throw_warning(software_dependencies, __LINE__, __FILE__, __FUNCTION__, "Did not find the tool 'free -b| grep \"Mem\"' on your system, therefore unable to predict the amount of free memory available. As a susbstitution, it's assumed that 16GB of data is avialable.");
      amount_of_free_memory_b = 16*1024*1024*2014;
    }
    return (amount_of_free_memory_b*0.5); // Halfes the sum in order to enable more use of it.
  }

  loint get_estimated_memory_usage_for_proteins(int taxon_length, float estimated_degree_filling_taxa, float estimated_degree_filling_proteins, loint cnt_proteins) {  
    return (loint)(estimated_degree_filling_taxa *taxon_length *cnt_proteins*sizeof(int)) // The length of the reference tables; each column has all of its proteins represented
      + (estimated_degree_filling_taxa*taxon_length *taxon_length * cnt_proteins
	 * cnt_proteins*estimated_degree_filling_proteins * sizeof(T)) // The containers holding the paris itself
      ;
  }
 public:

  /**
     @return true if place found fitting the expected data into the structure, else returns false.
  **/
  bool set_max_buffer_size_based_on_file_properties(taxa *obj, int taxon_length) {
    if(obj && taxon_length && false) { // TODO: remove the false-clause!
      const loint amount_of_free_memory_b = list_file_parse<T>::get_available_memory_on_system(); 
      log_builder::test_memory_condition_and_if_not_abort(listTaxa!=NULL, __LINE__, __FILE__, __FUNCTION__);
      const loint cnt_proteins = listTaxa[taxon_length-1].rel_end;
      const float estimated_degree_filling_proteins = 0.001; // TODO: Should be considered to be set from the terminal, specified by the user.
      // TODO: Use the info about the file size to calculate/estimate the value below!
      const float estimated_degree_filling_taxa = 0.5; // TODO: Should be considered to be set from the terminal, specified by the user.
      
      const loint size_est_in_b_taxon = sizeof(list_file_parse<T>)
	+ (estimated_degree_filling_taxa*taxon_length *taxon_length*sizeof(file_parse<T>)) // The number of the internal data strctures
	;
      const loint size_est_in_b_depend = get_estimated_memory_usage_for_proteins(taxon_length, estimated_degree_filling_taxa, estimated_degree_filling_proteins, cnt_proteins);
      
      
      const loint size_est_in_b = size_est_in_b_taxon + size_est_in_b_depend;
      const float size_est_in_gb = (float)size_est_in_b / (1024*1024*1024);
      
      if(amount_of_free_memory_b > size_est_in_b) {
	MAX_BUFFER_SIZE = UINT_MAX; // keeps everything in memory
	return true;
      } else if(amount_of_free_memory_b > size_est_in_b_taxon) {
	// Need here to adjust the number of proteins in the collection
	const loint free_space = amount_of_free_memory_b - size_est_in_b_taxon;
	assert(free_space < size_est_in_b_depend);
	const loint difference_to_be_overcomed = size_est_in_b_depend - free_space;
	loint cnt_proteins_temp = 0, cnt_taxa_temp = 0;
	
	int possible_numbers_of_protein_in_each = 1000;
	int updated_diff=get_estimated_memory_usage_for_proteins(taxon_length, estimated_degree_filling_taxa, estimated_degree_filling_proteins, possible_numbers_of_protein_in_each);
	uint cnt_iter = 0;
	// TODO: Not set correctly!
	do {
	  taxa::get_cnt_taxons_above_limit(obj, taxon_length, possible_numbers_of_protein_in_each, cnt_proteins_temp, cnt_taxa_temp);
	  int updated_diff=get_estimated_memory_usage_for_proteins(taxon_length, estimated_degree_filling_taxa, estimated_degree_filling_proteins, possible_numbers_of_protein_in_each);
	} while(updated_diff < difference_to_be_overcomed && cnt_iter++ < 1000);
	MAX_BUFFER_SIZE = possible_numbers_of_protein_in_each; // keeps everything in memory
	return true;
      } else {
	char temp[100]; memset(temp, '\0', 100);
	sprintf(temp, "Unable storing the complete taxa in the matrix, with the count=%d", taxon_length);
	log_builder::throw_warning(memory_constraint, __LINE__, __FILE__, __FUNCTION__, temp);
	return false; // Has not ability storing all the taxa in memory
      }
    } else return true;

  }
  /**
     @brief Inits a structure of this:
     @param <size> The number of objects to make for each class of this class type (class list_file_parse).
     @param <hashProtein> A reference to the hash-map linking protein labels to integers (i.e. only to be read).
     @param <_taxon_length> The number of different taxa in the collection.
     @param <_listTaxa> The properties of each taxon.
     @param <_FILE_BINARY_LOCATION> The location to store the binary (temporary) files.
     @param <USE_BEST_BLAST_PAIR_SCORE> If set to true, instead of summing the input values for the same protein pair, only uses the best pair.
     @author Ole Kristian Ekseth.
  **/
  static list_file_parse<T> ** init(uint size, prot_list *hashProtein, int _taxon_length, taxa_t* _listTaxa, char *_FILE_BINARY_LOCATION, bool USE_BEST_BLAST_PAIR_SCORE) {
    assert(_taxon_length >= 0);
    assert(size > 0);
    list_file_parse **parseData = NULL; // new list_file_parse<T>*[size]; //*)malloc(sizeof(list_file_parse*)*size); // Holds the refernces to the binary files created under the parsing
    try {parseData = new list_file_parse<T>*[size];;} 
    catch (std::exception& ba) {
      if(!log_builder::catch_memory_exeception(size, __FUNCTION__, __FILE__, __LINE__)) {
	fprintf(stderr, "!!\t An interesting error was discovered: %s."
		"The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		"Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
      }
    }
    for(uint my_id = 0;my_id<size;my_id++)
      parseData[my_id] = new list_file_parse(my_id, hashProtein, _taxon_length, _listTaxa, _FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);
    return parseData;
  }
  /**
     @brief Inits a structure of this:
     @param <_my_id> The id of this container. (To be used when mering the containers from multiple threads.)
     @param <_hash> A reference to the hash-map linking protein labels to integers (i.e. only to be read).
     @param <_taxon_length> The number of different taxa in the collection.
     @param <_listTaxa> The properties of each taxon.
     @param <_FILE_BINARY_LOCATION> The location to store the binary (temporary) files.
     @param <_USE_BEST_BLAST_PAIR_SCORE> If set to true, instead of summing the input values for the same protein pair, only uses the best pair.
     @author Ole Kristian Ekseth.
  **/
  list_file_parse(int _my_id, prot_list_t *_hash, int _taxon_length, taxa_t* _listTaxa, char *_FILE_BINARY_LOCATION, bool _USE_BEST_BLAST_PAIR_SCORE)   :
   taxon_length(_taxon_length), listTaxa(_listTaxa), FILE_BINARY_LOCATION(_FILE_BINARY_LOCATION),

    hashProtein(_hash), send_all_data_to_pipe(true),  my_id(_my_id), USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE), 
#ifdef USE_MPI
    list_of_nodes_taxa_responsilibties(NULL),
#endif
    MAX_BUFFER_SIZE(UINT_MAX)     {
    list = file_parse<T>::init(taxon_length);
  } 
  /**     Deallocates the memory  **/
  void free_mem(bool delete_file) {free_memory(delete_file);}
  //! Frees the reserved memory;
  void free_memory(const bool delete_file)  {
    file_parse<T>::close(list, taxon_length, delete_file);
#ifdef USE_MPI
    if(list_of_nodes_taxa_responsilibties) {delete [] list_of_nodes_taxa_responsilibties; list_of_nodes_taxa_responsilibties = NULL;}    
#endif
  }
  //! Frees the memory reserved for the list of file_parse object, but not the matrix holding 'potentials'.
  void free_file_parse_memory(const bool delete_file, uint index_start, uint index_end_pluss_one) {
    assert(index_start < (uint)taxon_length);
    assert(index_end_pluss_one <= (uint)taxon_length);
    file_parse<T>::close(list, delete_file, index_start, index_end_pluss_one);
  }
  //! Frees the reserved memory for the specified taxon-pair.
  void free_memory_for_taxon_pair(uint taxon_in, uint taxon_out)  {
    assert(taxon_in < (uint)taxon_length);
    assert(taxon_out < (uint)taxon_length);
    if(list && list[taxon_in] && list[taxon_in][taxon_out]) {
      list[taxon_in][taxon_out]->free_mem(true);
    }
  }

  //! Frees the reserved memory ontly for the non-inparalogs
  void free_memory_non_inparalogs(uint taxon_start, uint taxon_end)  {
    file_parse<T>::close_non_inparalogs(list, true, taxon_start, taxon_end);    
  }
  //! Frees the reserved memory ontly for the non-inparalogs
  void free_memory_For_inparalogs(uint taxon_start, uint taxon_end)  {
    file_parse<T>::close_For_inparalogs(list, true, taxon_start, taxon_end);    
  }
  /**
     @brief Closes the object, de-allocating the memory reserved.
     @param <obj> The object to close.
     @param <delete_file> If set, deletes intermediary files if they exists.
   **/
  static void close(list_file_parse *&obj, const bool delete_file) {
    if(obj) {obj->free_memory(delete_file); delete obj; obj = NULL;}
  }
  /**
     @brief Closes the list of objects used through parallel processing.
     @param <obj> The list of objects to close.
     @param <length> The number of objects the list contains.
     @param <delete_file> If set, deletes intermediary files if they exists.
   **/
  static void close(list_file_parse **&obj, uint length, const bool delete_file) {
    if(obj) {
      for(uint i =0; i < length; i++) {
	close(obj[i], delete_file); 
      }
      delete [] obj; obj = NULL;
    }
  }
  /**
     @brief De-allocates only the pointer to the class wrapper, and not the content itself
     @remarks Used after parallel processing, when one object consisted of multiple threads, and each thread allocated- and deallocated such objects.
  **/
  static void close_only_pointer(list_file_parse **&obj) {
    if(obj) {delete [] obj; obj = NULL;}
  }
  //! Constructor
  list_file_parse(prot_list_t *_hash, int _taxon_length, taxa_t* _listTaxa, char *_FILE_BINARY_LOCATION, bool _USE_BEST_BLAST_PAIR_SCORE) 
   : taxon_length(_taxon_length), listTaxa(_listTaxa), MODE_PAIRWISE_OUTPUT_ABC(false), MODE_PAIRWISE_OUTPUT_MCL(false),
    FILE_BINARY_LOCATION(_FILE_BINARY_LOCATION),
    hashProtein(_hash),
    send_all_data_to_pipe(true), my_id(0),
    USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE),
#ifdef USE_MPI
    list_of_nodes_taxa_responsilibties(NULL),
#endif
    MAX_BUFFER_SIZE(UINT_MAX)
    {
      list = file_parse<T>::init(taxon_length, _USE_BEST_BLAST_PAIR_SCORE); //init_list();
    }

  //! Constructor
 list_file_parse()   : taxon_length(0), listTaxa(NULL),
   MODE_PAIRWISE_OUTPUT_ABC(false), MODE_PAIRWISE_OUTPUT_MCL(false),
   FILE_BINARY_LOCATION(NULL), hashProtein(NULL), list(NULL), send_all_data_to_pipe(true), my_id(-1), USE_BEST_BLAST_PAIR_SCORE(false), 
#ifdef USE_MPI
    list_of_nodes_taxa_responsilibties(NULL),
#endif
   MAX_BUFFER_SIZE(UINT_MAX)
   {
  }


  //! Constructor
   list_file_parse(int _taxon_length, taxa_t* _listTaxa, bool _USE_BEST_BLAST_PAIR_SCORE)   :
     taxon_length(_taxon_length), listTaxa(_listTaxa),
     MODE_PAIRWISE_OUTPUT_ABC(false), MODE_PAIRWISE_OUTPUT_MCL(false),
     FILE_BINARY_LOCATION(NULL), hashProtein(NULL), list(NULL),
     send_all_data_to_pipe(true), my_id(-1), USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE) {
  }

  //! Constructor
 list_file_parse(list_file_parse &arg)   : hashProtein(NULL), list(NULL), send_all_data_to_pipe(true),
   MODE_PAIRWISE_OUTPUT_ABC(false), MODE_PAIRWISE_OUTPUT_MCL(false),
   my_id(-1) 
#ifdef USE_MPI
   , list_of_nodes_taxa_responsilibties(NULL)
#endif
   { 
    steal_data(arg);

  } 

  //! Constructor, called from the merger, i.e. the serial part of the code
   list_file_parse(prot_list_t *_hash, long int max_buff_size, int _taxon_length, taxa_t* _listTaxa, char *_FILE_BINARY_LOCATION, bool _USE_BEST_BLAST_PAIR_SCORE)   :
     taxon_length(_taxon_length), listTaxa(_listTaxa),
     MODE_PAIRWISE_OUTPUT_ABC(false), MODE_PAIRWISE_OUTPUT_MCL(false),
     FILE_BINARY_LOCATION(_FILE_BINARY_LOCATION),
     hashProtein(_hash), send_all_data_to_pipe(true),
     my_id(0), USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE),
#ifdef USE_MPI
    list_of_nodes_taxa_responsilibties(NULL),
#endif
     MAX_BUFFER_SIZE(max_buff_size) 
    {
     list = file_parse<T>::init(taxon_length);
  } 

  /**! Allocates data  */
  static list_file_parse<T> **allocate_single_list(uint size, uint size_each, int _taxon_length, taxa_t *_listTaxa, 
						 bool _MODE_PAIRWISE_OUTPUT_ABC,  // for the abc file: if set, the data out pairwise in stead of as in a row
						 bool _MODE_PAIRWISE_OUTPUT_MCL,  // for the mcl file: if set, the data out pairwise in stead of as in a row
						 char *_FILE_BINARY_LOCATION)  {
    assert(_taxon_length >= 0);
    assert(size > 0);
    list_file_parse **parseData = new list_file_parse<T>*[size]; //*)malloc(sizeof(list_file_parse*)*size); // Holds the refernces to the binary files created under the parsing
    test_condition_and_if_not_abort(parseData != NULL, __LINE__, __FUNCTION__);
    for(uint my_id = 0;my_id<size;my_id++) {
      parseData[my_id] = list_file_parse<T>::init_class(_taxon_length, _listTaxa, _MODE_PAIRWISE_OUTPUT_ABC, _MODE_PAIRWISE_OUTPUT_MCL, _FILE_BINARY_LOCATION);
    }
    //  printf("in list_file_struct_depricated allocated (%u, %u) at line %d\n", n_threads, size_each, __LINE__);
    return parseData;
  }


  //! Initiates the class.
  static list_file_parse<T> *init_class(int _taxon_length, 
				      taxa_t *_listTaxa, 
				      bool _MODE_PAIRWISE_OUTPUT_ABC,  // for the abc file: if set, the data out pairwise in stead of as in a row
				      bool _MODE_PAIRWISE_OUTPUT_MCL,  // for the mcl file: if set, the data out pairwise in stead of as in a row
				      char *_FILE_BINARY_LOCATION) {
    return new list_file_parse<T>(UINT_MAX, _taxon_length, _listTaxa, _MODE_PAIRWISE_OUTPUT_ABC, _MODE_PAIRWISE_OUTPUT_MCL, _FILE_BINARY_LOCATION);
  }

  //! The constructor:
  list_file_parse(uint max_buff_size, int _taxon_length, 
		   taxa_t *_listTaxa, 
		   bool _MODE_PAIRWISE_OUTPUT_ABC,  // for the abc file: if set, the data out pairwise in stead of as in a row
		   bool _MODE_PAIRWISE_OUTPUT_MCL,  // for the mcl file: if set, the data out pairwise in stead of as in a row
		   char *_FILE_BINARY_LOCATION)   : 
     taxon_length(_taxon_length), listTaxa(_listTaxa),
    MODE_PAIRWISE_OUTPUT_ABC(_MODE_PAIRWISE_OUTPUT_ABC), MODE_PAIRWISE_OUTPUT_MCL(_MODE_PAIRWISE_OUTPUT_MCL),
     FILE_BINARY_LOCATION(_FILE_BINARY_LOCATION),
     hashProtein(NULL), send_all_data_to_pipe(true),  my_id(0),
     USE_BEST_BLAST_PAIR_SCORE(false), 
#ifdef USE_MPI
    list_of_nodes_taxa_responsilibties(NULL),
#endif
MAX_BUFFER_SIZE(max_buff_size)
      {
	list = file_parse<rel>::init(taxon_length, false);
      }

#ifdef USE_MPI
    /**
       @brief Sends only selective data accross nodes.
       @todo Requires that at least the same number of taxa must be found in the file as the number of nodes used for this processing.
    **/
    void mpi_send_only_selective_data(MPI_Comm comm, int myrank, int number_of_nodes, int t_taxon_length);
    /**
       @brief Sends- and receives the list of p_rel objects accross nodes.
    **/
    void mpi_make_data_consistent_accross_nodes(MPI_Comm comm, int myrank, int number_of_nodes, uint index_start, uint index_end_pluss_one);

    /**
       @brief Sends- and receives the list of rel objects accross nodes.
       @remarks Each node only sends- and receives those taxon-pairs regarded as interesting in the building of co-orthologs.
**/
    void mpi_make_data_consistent_accross_nodes(int taxon_length);
#endif


  //! Builds a list for testing.
  void build_default_list(uint length_prots)  {
#ifdef assert_code
    const uint block_number  = 0;
    const uint mem_reserved = UINT_MAX;
    if(taxon_length > 0) {
      list = new file_parse<T>**[taxon_length];//parseData_c[my_id];
      for(uint i =0; i<(uint)taxon_length; i++) {
	list[i] = new file_parse<T>*[taxon_length]; 
	for(uint out =0; out<(uint)taxon_length; out++) {
	  list[i][out] = new file_parse<T>(i, out, length_prots, block_number, mem_reserved, FILE_BINARY_LOCATION, false);
	}
      }
    } else fprintf(stderr, "!!\tTaxon length not set: Discard initializing the list in 'list_file_parse' class\n");
#endif
  }
  //! The main test function for this class  
  static void assert_class(const bool print_info) {
    const static char *class_name = "list_file_parse";
    if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
#ifndef NDEBUG
    //  const uint old_taxon_length = taxon_length;
    //  taxa *old_listTaxa = listTaxa;
    int taxon_length = 3;   const int length_in = 3, length_out = 3;
    const uint length_prots = 5;
    taxa_t *listTaxa = taxa::init_test_array(taxon_length, length_prots);
    list_file_parse *data = new list_file_parse<T>(taxon_length, listTaxa, false); 
    data->build_default_list(length_prots);// builds some sets of file_parse-objects:
    for(uint in =0; in < (uint)taxon_length; in++) {
      for(uint out =0; out < (uint)taxon_length; out++) {
	for(int prot_in = 0; prot_in< length_in; prot_in++) {
	  for(int prot_out = 0; prot_out< length_out; prot_out++) {
	    Parse p = Parse(); p.set_distance(2), p.set_overlap_in(4);
	    protein_relation pr = protein_relation();
	    pr.protein_in = prot_in, pr.protein_out = prot_out, pr.taxon_in = in, pr.taxon_out = out;
	    bool ret = data->insert_new_rel(in, out, p, pr);
	    if(false) assert(ret == !((in == out) && (prot_in == prot_out)));
	    if(prot_in != prot_out) {
	      assert(data->has_data(in, out));
	      assert(data->has_data(in, out, prot_in));
	    }
	    p.free_memory();
	  }
	}
      }     
    }
    data->free_mem(true);
    taxa::delete_taxa(taxon_length, listTaxa);
    delete data;
#endif
#endif
    if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
    //  if(false)   print_constants();

  }

  };

/**
   @brief Stores the blast file with information about overlap-values.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth).
**/
typedef class list_file_parse<p_rel> list_file_parse_t;

/**
   @brief Stores the blast file with information about overlap-values.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth).
**/
typedef class list_file_parse<rel> list_file_struct_t;


#endif
 
