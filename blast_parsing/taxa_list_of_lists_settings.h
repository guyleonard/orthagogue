#ifndef taxa_list_of_lists_settings_h
#define taxa_list_of_lists_settings_h
/** 
   @file
   @brief Holds the data containers used for mpi sending- and receiving.
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
#include "taxon_list_t.h"
#include "taxon_data.h"
#ifdef USE_MPI
#include "taxa_buffer_list_settings.h"
/**
   @struct taxa_list_of_lists_settings
   @brief Holds the data containers used for mpi sending- and receiving.
   @remarks 
   # Seperated for structure taxa_taxa_buffer_list_settings in order to emobdy the seperation between what's sent seperately, and what is sendt wrapped as an mpi container.
   # Only acitvated when amcro-variable USE_MPI is set.
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
   @date 10.07.2012 by oekseth (initial)
**/
struct taxa_list_of_lists_settings {
  char *taxonNamesArr; // For storing all the taxon-names in the same 1-d-list, accepting variable-sizes of the char-strings
  uint *taxon_prots_used; // The number of proteins for each taxa.
  char *proteinLabels;// The buffer of protein labels
  uint *bufferProteinsLength; // The list of protein label-sizes.

  //! Deallocates the memory:
  void free_memory() {
    if(taxonNamesArr) {delete [] taxonNamesArr, taxonNamesArr = NULL;}
    if(taxon_prots_used) {delete [] taxon_prots_used; taxon_prots_used = NULL;}
    if(proteinLabels) {delete [] proteinLabels; proteinLabels = NULL;}
    if(bufferProteinsLength) {delete [] bufferProteinsLength; bufferProteinsLength = NULL;}
  }

  //! Sends the data to the root.
  void send_data(MPI_Comm comm, int myrank, int tag, taxa_buffer_list_settings_t obj) {
    const int root_id = 0; send_data(comm, myrank, tag, obj, root_id);
  }

  //! Sends the data to given node.
  void send_data(MPI_Comm comm, int myrank, int tag, taxa_buffer_list_settings_t obj, int receive_id) {
    int ret_val = 0;
    ret_val = MPI_Send(taxonNamesArr, obj.taxonNamesArr_size, MPI_CHAR, receive_id, tag, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 

    ret_val = MPI_Send(taxon_prots_used, obj.taxon_size, MPI_UNSIGNED, receive_id, tag, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 

    ret_val = MPI_Send(proteinLabels, obj.total_cnt_chars_used_representing_proteins, MPI_CHAR, receive_id, tag, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 

    ret_val = MPI_Send(bufferProteinsLength, obj.total_cnt_proteins, MPI_UNSIGNED, receive_id, tag, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 
  }

  //! Sends the data using broad_cast method, ie, to every node
  void broadcast_send_data(MPI_Comm comm, int myrank, int number_of_nodes, int tag, taxa_buffer_list_settings_t obj) {
    int ret_val = 0;
    assert(myrank == 0); // The assumption when the function was first used. {
    const int root_id = myrank;
    ret_val = MPI_Bcast(taxonNamesArr, obj.taxonNamesArr_size, MPI_CHAR, root_id, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 
    ret_val = MPI_Bcast(taxon_prots_used, obj.taxon_size, MPI_UNSIGNED, root_id, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 

    ret_val = MPI_Bcast(proteinLabels, obj.total_cnt_chars_used_representing_proteins, MPI_CHAR, root_id, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 

    ret_val = MPI_Bcast(bufferProteinsLength, obj.total_cnt_proteins, MPI_UNSIGNED, root_id, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 
    
  }

  //! Merges the buffers into the taxon_list object:
  void merge(taxon_list_t *this_sender, taxa_buffer_list_settings_t obj) {
    assert(this_sender);
    uint index_startpos = 0, buffer_startpos = 0, taxon_name_index = 0;
    for(uint i = 0; i < obj.taxon_size; i++) {
      //! Validate the transformation of this new-built lists into the original type (ie, a set of taxon_data objects).
      taxon_data taxon_obj = taxon_data(taxon_prots_used[i]);
      //! Sets the taxon label:
      taxon_obj.setTaxonName(taxonNamesArr + taxon_name_index);

      //! Sets the set of protein labels:
      char *listOfProteins = proteinLabels    + buffer_startpos;
      uint *listOfLengths  = bufferProteinsLength + index_startpos;
      const uint protein_length = taxon_prots_used[i];
      buffer_startpos += taxon_obj.set_labels_using_1d_object(protein_length, listOfProteins, listOfLengths);

      //! Merges the object into the xsisting one:
      this_sender->merge(taxon_obj);

      //! Update the variables for the next run:
      taxon_name_index += strlen(taxonNamesArr + taxon_name_index) +1; // Taxon-name-variables
      index_startpos += taxon_prots_used[i]; // Protein-label-variables.
    }    
  }
  //! Receives and merges the data int the given object
  void receive_and_merge_data(taxon_list_t *this_sender, MPI_Comm comm, int myrank, int sender_id, int tag, taxa_buffer_list_settings_t obj) {
    //! Receives the data:
    int ret_val = 0;
    if(myrank == 0) {
      MPI_Status status_object_info;
      ret_val = MPI_Recv(taxonNamesArr, obj.taxonNamesArr_size, MPI_CHAR, sender_id,  tag, comm, &status_object_info);
      if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 
      
      ret_val = MPI_Recv(taxon_prots_used, obj.taxon_size, MPI_UNSIGNED, sender_id,  tag, comm, &status_object_info);
      if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 
      
      ret_val = MPI_Recv(proteinLabels, obj.total_cnt_chars_used_representing_proteins, MPI_CHAR, sender_id,  tag, comm, &status_object_info);
      if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 
      
      ret_val = MPI_Recv(bufferProteinsLength, obj.total_cnt_proteins, MPI_UNSIGNED, sender_id,  tag, comm, &status_object_info);
      if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__);
    } else {
      // TODO: If principles changes to more advaced sheduling in the futre, a change here is mandatory.
      const int root_id = sender_id;
      ret_val = MPI_Bcast(taxonNamesArr, obj.taxonNamesArr_size, MPI_CHAR, root_id, comm);
      if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 
      ret_val = MPI_Bcast(taxon_prots_used, obj.taxon_size, MPI_UNSIGNED, root_id, comm);
      if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 

      ret_val = MPI_Bcast(proteinLabels, obj.total_cnt_chars_used_representing_proteins, MPI_CHAR, root_id, comm);
      if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 

      ret_val = MPI_Bcast(bufferProteinsLength, obj.total_cnt_proteins, MPI_UNSIGNED, root_id, comm);
      if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__); 
    }

    //! Merges the data:
    merge(this_sender, obj);
  }
  //! Constructor setting everything to NULL
  taxa_list_of_lists_settings() : taxonNamesArr(NULL), taxon_prots_used(NULL), proteinLabels(NULL), bufferProteinsLength(NULL) {};

  //! Constructor with objects
  taxa_list_of_lists_settings(char *_taxonNamesArr, uint *_taxon_prots_used, char *_proteinLabels, uint *_bufferProteinsLength) : taxonNamesArr(_taxonNamesArr), taxon_prots_used(_taxon_prots_used), proteinLabels(_proteinLabels), bufferProteinsLength(_bufferProteinsLength) {};

  //! Constructor with lengths for allocating objects, initialising them, with intention of getting data through mpi receive calls.
  taxa_list_of_lists_settings(taxa_buffer_list_settings_t obj) : taxonNamesArr(NULL), taxon_prots_used(NULL), proteinLabels(NULL), bufferProteinsLength(NULL)
  {
    //! Allocates in the same order as the varialbes are listed:
    if(obj.taxonNamesArr_size) {
      taxonNamesArr = new char[obj.taxonNamesArr_size+1];
      memset(taxonNamesArr, '\0', obj.taxonNamesArr_size+1);
    }
    if(obj.taxon_size) {
      taxon_prots_used = new uint[obj.taxon_size];
      memset(taxon_prots_used, 0, obj.taxon_size);
    }
    if(obj.total_cnt_proteins) {
      bufferProteinsLength = new uint[obj.total_cnt_proteins];
      memset(bufferProteinsLength, 0, obj.total_cnt_proteins);
    }
    if(obj.total_cnt_chars_used_representing_proteins) {
      proteinLabels = new char[obj.total_cnt_chars_used_representing_proteins+1];
      memset(proteinLabels, '\0', obj.total_cnt_chars_used_representing_proteins+1);
    }
  };
};

/**
   @struct taxa_list_of_lists_settings
   @brief Holds the data containers used for mpi sending- and receiving.
   @ingroup parsing_container
   @author Ole Kristian Ekseth (oekseth)
   @date 10.07.2012 by oekseth (initial)
**/
typedef struct taxa_list_of_lists_settings taxa_list_of_lists_settings_t;

#endif

#endif
