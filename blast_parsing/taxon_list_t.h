#ifndef taxon_list_h 
#define taxon_list_h
/** 
   @file
   @brief Builds a linked list data structure.
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
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h> // some systems require it
#include <sys/stat.h>
#include <sys/termios.h> // for winsize
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/times.h>
#include "taxon_data.h"
#include "../configure.h"
#include "tbb_libs.h"
#include "log_builder.h"
#ifdef USE_MPI
#include "taxa_buffer_list_settings.h"
#endif
/**
   @class taxon_list
   @brief Builds a linked list data structure.
   @ingroup parsing_container
   @remark The task of this class is to handle protein names, bulding linked list,
   sorted on the diffe rent taxa. An example is given below:
   list[human] = (protein_1) -> (protein_2) -> ...->(protein_n)
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)  
**/
class taxon_list {
 private:
  //! @param <size> the size of the name-pointer to be (1) allocated, (2) Intialized to '\0'
  char* init_name(uint size);
  //! Inserts a new taxon into the list
  uint insert_taxon(char *taxon);
  /**
     @brief Enlarges the taxon list if not big enough
     @date 14.01.2011 by oekseth
  */
  void enlarge_taxon_list();
  /**@return the listTaxon given. */
  taxon_data_t *listTaxon;
  //! The name of the taxon first read: Stored with regard to merging of the list.
  char *taxon_name_first_read; 
  loint taxon_name_first_read_size;
  //! The name of the taxon last read: Stored with regard to merging of the list.
  char *taxon_name_last_read;
  loint taxon_name_last_read_size;
  //! The name of the protein first read: Stored with regard to merging of the list.
  char *protein_name_first_read;
  loint protein_name_first_read_size;
  //! The name of the protein last read: Stored with regard to merging of the list.
  char *protein_name_last_read;
  loint protein_name_last_read_size;
  //! The number of taxa used.
  int taxon_used;
  //! The number of taxa reserved
  int taxon_reserved;
 public:


  //! @return the name of the last protein inserted
  char *get_protein_name_last_read() {return protein_name_last_read;}
  //! @return the length of the last protein inserted.
  loint get_protein_name_last_read_length() {
    return protein_name_last_read_size;
/*     if(protein_name_last_read) return strlen(protein_name_last_read); */
/*     else return 0; */
  }
  //! @return the name of the last taxon inserted
  char *get_taxon_name_last_read() {return taxon_name_last_read;}
  //! @return the length of the last taxon inserted.
  loint get_taxon_name_last_read_length() {
    return taxon_name_first_read_size;
/*     if(taxon_name_last_read) return strlen(taxon_name_last_read); */
/*     else return 0; */
  }
  //! @return true if equal
  bool is_equal_taxon_name_first(char *name) {
    return compare_strings(taxon_name_first_read, name);
  }
  //! @return true if equal
  bool is_equal_protein_name_first(char *name) {
    return compare_strings(protein_name_first_read, name);
  }
  //! Copies protein name from 'src' to 'dest', updating the 'dest_size' variable.
  void copy_name(char *&dest, loint &size_dest, char *src, const loint size_src);
  //! Copies taxon name from 'src' to 'dest'.
  void copy_taxon_name(char *&dest, char *src);
  /**
     @brief Inserts the name in the taxon_data object list *buffTaxon
     @param <taxon> The label of the taxon.
     @param <name> The label of the protein.
     @return True if a new taxon is inserted.
  **/
  bool insert_protein(char *taxon, char *name);
  
  //! If argument set to true, prints the proteins as well
  void printf_buffer(const bool print_prots);

  /**
     @brief Produces (writes) a log file holding the list of memory allocations.
     @remarks Useful for optimizing the memory allocation procedures for building of lists for different taxa.
  **/
  void log_produce_memory_allocations();

  /**
     @return the number of taxa in teh collection
     @date 01.12.2010
     @deprecated Replaced by getTaxaLength
  **/
  uint getLength();

  //! @return the number of taxa in teh collection
  uint getTaxaLength() {return getLength();}
  //! @return the number of proteins for the given taxon
  uint getProteinLength(uint i);

  //! @return the given buffer from 'lisTaxon'
  char **getBuffer(uint i);
  /**
     @brief  @return the given name from listTaxon   
     @date 02.12.2010
     @deprecated Replaced by getTaxonName
  **/
  char *getName(uint taxon_index);
  //! @return the given name from listTaxon   
  char *getTaxonName(uint taxon_index) {return getName(taxon_index);}
  //! Deletes the pointers/memory references in this structure
  void delete_buffer();
  /**
     @brief Merges two arrays.
     @remark Must be in the same order
  */
  void merge(taxon_list *arg);

  /**
     @brief Merges an element into the list    
     @param <arg> The input to merge.
     @remarks Assumes there is no overlap between proteins due to disjoint regions.
  */
  void merge(taxon_data &arg);

  /**
     @brief Takes a lsit of 'taxon_data', and for those it has a zero relation with, merges
     @remark Called from the reciprocal method 'merge': 
  */
  void recipZeroYEsMerge(taxon_data *&arg, uint outer_taxon_used);
  /**
     @return true if the taxon given as argument is the same as the first
     taxon in both the argument and 'this', and in addition that
     the proteins are equal, and therefore the last shall be skipped.
  */
  bool is_overlap(taxon_list arg, char *taxon);


  //! @return true if this lsit has the given taxon set as argument
  bool hasTaxon(char *taxon_name, taxon_data &buff);
  /**
     @param <taxon_name> The name of the taxon to search for.
     @param <taxon_name_index> The index to be updated in this method.
     @return True if if its found for the given param.     
     @remark Used by getProteinLength(char *taxon)
  **/
  bool getTaxonIndex(char *taxon_name, int &taxon_name_index);
  /**
     @param <taxon> The name of the taxon to search for.
     @param <taxon_name_index> The index to be updated in this method.
     @return True if if its found for the given param.     
     @remark Assumes that the taxon exsits: If not, returns the 0-index(0)
  **/
  bool getProteinLength_n(char *taxon, int &taxon_name_index);  

  /**
     @brief Speeds up the exectuion for files having a bunch of taxa with close to empty collections.
     @param <threshold> The number of proteins a taxon must have more of in order to be evaluated.
     @param <taxon_length> The variable to update.
  **/
  void remove_taxa_with_proteins_below_threshold(uint threshold, int &taxon_length);

#ifdef USE_MPI
 private:
  //! Builds an mpi object, including it into the two structures given as parameters.
  void build_mpi_1d_lists(taxa_buffer_list_settings_t &obj_lists_length, char *&taxonNamesArr, uint *&taxon_prots_used, char *&proteinLabels, uint *&bufferProteinsLength, int myrank, int number_of_nodes);
 public:
  /**
     @brief Sends- and receives the list of taxa- and protein labels accross nodes.
  **/
  void mpi_make_data_consistent_accross_nodes(MPI_Comm comm, int myrank, int number_of_nodes, const uint LIMIT_MINIMUM_NUMBER_OF_PROTEINS_FOR_EACH_TAXA, int &taxon_length_to_be_updated);

#endif

  /**
     @brief De-allocates only the pointer to the class wrapper, and not the content itself
     @remarks Used after parallel processing.
  **/
  static void close_only_pointer(taxon_list *&obj) {
    if(obj) {delete [] obj; obj = NULL;}
  }
  //! Deallocates the memory given the param.
  static void close(taxon_list *&obj);
  //! Frees the memory
  void free_mem();
  //! Frees the memory
  void free_memory() {free_mem();}
  //! @return an object of this class.
  static taxon_list *init();
  //! The constructor
  taxon_list();

 private:
  //! Verifies the length set of the taxa and the prots
  void assert_getLength();
  //!  @return True if the strings given are equal.
  bool compare_strings(char *name_1, char *name_2);
  //!  @return True if the strings given are equal.
  bool compare_strings(char *name_1, char *name_2, uint length_string);
  /**
     @remarks 
   - Used to verify that the buffers are equal
   - Uses pointers, meaning extra vlidation testing they are set is included.
   @return true if the input buffers are equal on all of their elements.
  **/
  bool compare_buffers(char **buff_1, char **buff_2, uint length_buffer, uint length_string);
  //! Tests that the insertion of the proteins works correctly
  void assert_insert_protein();
 public:
  //! The main test function for this class
  static void assert_class(const bool show_positive_results);
};

/**
   @brief Builds a linked list data structure.
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
   @date 16.09.2011 by oekseth (asserts)  
**/
typedef class taxon_list taxon_list_t;

#endif
