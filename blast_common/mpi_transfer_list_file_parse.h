#ifndef mpi_transfer_list_file_parse_h
#define mpi_transfer_list_file_parse_h
#ifdef USE_MPI
/**
   @file
   @brief Holds the enitre list_file_parse object
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
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

#include "log_builder.h"
#include "mpi_transfer_p_rel.h"
#include "mpi_transfer_rel.h"
#include "list_file_parse.h"

#include "uint_3x_obj.h"
/**
   @class mpi_transfer_list_file_parse
   @brief Holds the enitre list_file_parse object
   @tparam <T> A rel class for the relations_list object, ie, the implemented p_rel or rel classes
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
template<class T> class mpi_transfer_list_file_parse {    
 private:
  int myrank;
  int number_of_nodes;
  list_file_parse<T> *this_caller;
  uint taxon_length;
  taxa *listTaxa;
  MPI_Datatype mpi_type_uint_3x; 
  MPI_Datatype mpi_type_uint_2x;
  MPI_Datatype mpi_type_T;
#ifndef NDEBUG
  //! Below variable to validate that receiving of data is correct, and thereby conclude that any difference with our expectations is due to mering of data): If in debug, sets the length of distances to transfer to each node, later used validating that the received value correspons to the actual.
  float *debug_sum_distances_total_to_send;
  float *debug_sum_distances_total_received;
#endif

  //! The different types of info objects:
  enum send_types {object_info, object_info_in_list, index_list, T_list, debug_sum_distance};
  //! Specifiers for the type of sending operation to use:
  enum collective_type {mpi_bcast, mpi_send_blocking, mpi_Isend_nonblocking};
  //! Gets the 2d-coordinate given the id.
  void get_2d_coordinate(uint id, uint &in, uint &out, uint taxon_length) {
    if(id) {
      assert(taxon_length);
      out = id % taxon_length;
      //    in = 0; printf("--\t id(%u) --> in(%u) and out(%u) at line %d in file %s\n", id, in, out, __LINE__, __FILE__);
      assert(id >= out);
      in = (id - out) / taxon_length;
    } else {in = 0, out = 0;}
    //  printf("--\t id(%u) --> in(%u) and out(%u) (taxon_length=%u) at line %d in file %s\n", id, in, out, taxon_length, __LINE__, __FILE__);
  }

  //! @return the 1d-coordinate given the id.
  uint get_1d_coordinate(uint in, uint out, uint taxon_length) {  
    const uint id = ((in*taxon_length)+out);
#ifndef NDEBUG
    //! Asserts that the element-id corresponds to the correct inverse transformation:
    uint in_temp, out_temp;	  
    get_2d_coordinate(id, in_temp, out_temp, taxon_length);
    assert(in == in_temp);
    assert(out == out_temp);
#endif
    return id;
  }


  //! @return the 1d node coordinate for a given taxon.
  uint get_1d_node_taxon_coordinate(uint node_rank, uint taxon_id, uint taxon_length) {
    return (node_rank * taxon_length) + taxon_id; 
  }
 private:
  /**
     @brief Build a datatype consisting of 2 uint-elements.
     @param <mpi_type> The datatype to set.
  **/
  void mpi_build_datatype_2_uint(MPI_Datatype &mpi_type_uint_3x) {
    static const int cnt_types = 2;// The number of types to represent.
    int block_lengths[cnt_types] = {1,1};
    MPI_Datatype old_types[cnt_types] = {MPI_UNSIGNED, MPI_UNSIGNED};
    //! The location of each element;
    // TODO: Is the below "mpi initialisation" correct?
    MPI_Aint indices[cnt_types] = {0, sizeof(uint)};
    for(int i = 1; i < cnt_types;i++) {indices[i] +=indices[i-1];}
    //! Builds the mpi representation of the structure:
    MPI_Type_struct(cnt_types, block_lengths, indices, old_types, &mpi_type_uint_3x);
    MPI_Type_commit(&mpi_type_uint_3x);
  }


  /**
     @brief Build a datatype consisting of 3 uint-elements.
     @param <mpi_type> The datatype to set.
  **/
  void mpi_build_datatype_3_uint(MPI_Datatype &mpi_type_uint_3x) {
    static const int cnt_types = 3;// The number of types to represent.
    int block_lengths[cnt_types] = {1,1, 1};
    MPI_Datatype old_types[cnt_types] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED};
    //! The location of each element;
    // TODO: Is the below "mpi initialisation" correct?
    MPI_Aint indices[cnt_types] = {0, sizeof(uint), sizeof(uint)};
    for(int i = 1; i < cnt_types;i++) {indices[i] +=indices[i-1];}
    //! Builds the mpi representation of the structure:
    MPI_Type_struct(cnt_types, block_lengths, indices, old_types, &mpi_type_uint_3x);
    MPI_Type_commit(&mpi_type_uint_3x);
  }

 private:
  /**
     @brief Does preprocessing, to know the length of the lists, avoiding to resize:
     @remarks 
     # Gives some extra costs, but hope the reduced memory consumption may benefit it.
     # Note that "total_length_of_index_data" is not a subsitute for the value found in
     "listTaxa", as the varibalbe may contain less proteins than the maximum allowable.
  **/
  void get_metadata_information(uint index_start, uint index_end_pluss_one, uint &cnt_taxa_pairs, uint &cnt_total_length_of_rel_data, uint &cnt_total_length_of_index_data, list_file_parse<T> *this_caller) {
    cnt_taxa_pairs = 0, cnt_total_length_of_rel_data = 0, cnt_total_length_of_index_data = 0;
    file_parse<T> ***list = this_caller->get_buffer();
    if(list) {
      for(uint i = index_start; i < index_end_pluss_one; i++) {    
	if(list[i]) {
	  for(uint out = index_start; out < index_end_pluss_one; out++) {
	    if(list[i][out] && list[i][out]->has_data()) {
	      cnt_taxa_pairs++;
	      cnt_total_length_of_index_data += list[i][out]->get_index_length();
	      cnt_total_length_of_rel_data   += list[i][out]->getTotalLengthOfData();
	    }
	  }
	}
      }
      //! Verifies that our count is correct:
      assert(cnt_total_length_of_rel_data == this_caller->getTotalLengthOfData());
      assert(cnt_total_length_of_index_data == this_caller->get_index_length());
    }
  }


  //! @return true if data is something the receiving node is interested in:
  bool receive_node_is_interested_in_these_data(uint id_receive_node, uint in, uint out, uint taxon_length, uint *data_to_send) {
    //! Calculates the data the receiver is interested in
    assert(taxon_length);
    assert(in < taxon_length);
    assert(out < taxon_length);
    assert(data_to_send);
    const uint index_id_receiver_innermost = get_1d_node_taxon_coordinate(id_receive_node, in, taxon_length); 
    const uint index_id_receiver_outermost = get_1d_node_taxon_coordinate(id_receive_node, out, taxon_length);
    if((data_to_send[index_id_receiver_innermost] != 0) || (data_to_send[index_id_receiver_outermost] != 0)) { // Data the receiver is interested in:
      return true;
    } else return false;
  }
  //! @return information about 'this' object:
  uint *get_metadata_information(struct uint_3x_obj &obj_info, uint &cnt_received_rel_not_reciprocal) {
    return get_metadata_information(myrank, NULL, taxon_length, number_of_nodes, obj_info, this_caller, cnt_received_rel_not_reciprocal);
  }

  /**
     @brief Does preprocessing, to know the length of the lists, avoiding to resize:
     @param <id_receiver> The id (rank) of the receiver.
     @param <lst_lengths> The list holding the lengths (making it faster doing the calculations, ie, the data the receiver is interested in
     @param <obj_info>    The object to insert the metadata into
     @param <this_caller> The mother object to use.
     param <cnt_received_rel_not_reciprocal> Used for the receving nod eto validate that the numer of paris received is correct.
     @return The list holding the imaginary next position in file, used to avoid sending data received from others
     @remarks 
     # Tested using similar to the function (the one with same naame, but different set of parameters).
     # If 'list_of_data_to_send' is not set, assumes all data is to be sent:
  **/
  uint *get_metadata_information(uint id_receiver, uint *list_of_data_to_send, uint taxon_length, int number_of_nodes, struct uint_3x_obj &obj_info, list_file_parse<T> *this_caller, uint &cnt_received_rel_not_reciprocal) {
    // Asserts, intitilises, and gets easy access to data:
    cnt_received_rel_not_reciprocal = 0;
    assert(taxon_length);
    //    assert(list_of_data_to_send);
    uint *lastPosInOwnData = new uint[taxon_length*taxon_length];
    log_builder::test_memory_condition_and_if_not_abort(lastPosInOwnData!=NULL, __LINE__, __FILE__, __FUNCTION__);
    memset(lastPosInOwnData, 0, sizeof(uint)*taxon_length*taxon_length);
    uint cnt_taxa_pairs = 0, cnt_total_length_of_rel_data = 0, cnt_total_length_of_index_data = 0;
    file_parse<T> ***list = this_caller->get_buffer(); // For easy access.
    if(list) {
      for(uint i = 0; i < taxon_length; i++) {    
	if(list[i]) {
	  for(uint out = 0; out < taxon_length; out++) {
	    if(list[i][out] && list[i][out]->has_data()) {	      
#ifndef NDEBUG
	      const uint index_id_receiver_innermost = get_1d_node_taxon_coordinate(id_receiver, i, taxon_length); 
	      //! Value appended (for the caller) in debug mode for the receiver, ensuring that the new vount is correct.
	      if(!list_of_data_to_send || (list_of_data_to_send[index_id_receiver_innermost] != 0)) {
		cnt_received_rel_not_reciprocal += list[i][out]->getTotalLengthOfData();
	      }
#endif
	      //! Calculates the data the receiver is interested in
	      if(!list_of_data_to_send || receive_node_is_interested_in_these_data(id_receiver,i,out, taxon_length, list_of_data_to_send)) {
		cnt_taxa_pairs++;
		/**
		   As the list of index is enlarged (when receiving),
		   an easy solution solving this be to use the number of
		   proteins for the given taxa as basis. 
		   (
		   An alternative call would have been cnt_total_length_of_index_data += list[i][out]->get_index_length(); 
		   )
		**/
		cnt_total_length_of_index_data += list[i][out]->get_index_reserved();
		cnt_total_length_of_rel_data   += list[i][out]->getTotalLengthOfData();
		//! Sets the position; avoid sending received data back- and forth.
		lastPosInOwnData[get_1d_coordinate(i, out, taxon_length)] = list[i][out]->get_index_T_end_of_this_object();
	      } else lastPosInOwnData[get_1d_coordinate(i, out, taxon_length)] = 0;

	    } else lastPosInOwnData[get_1d_coordinate(i, out, taxon_length)] = 0;
	  }
	}
      }
    }
    //! Initiates local variables, making the code similar to the code sed sending the objects:
    obj_info.taxon_id = cnt_taxa_pairs; // 'taxon_id' represents there the total number of taxon-pairs.
    obj_info.total_length_of_index = cnt_total_length_of_index_data;
    obj_info.total_length_of_rel = cnt_total_length_of_rel_data;
    return lastPosInOwnData;
  }

  //! If macro variable NDEBUG is not set, then function verifies that they are equal
  void verify_is_equal_and_set(index_t *arr_index_current, index_t *temp_index, uint index_length, uint taxon_in, uint taxon_out, file_parse<T> ***list) {
#ifndef NDEBUG
    assert(arr_index_current);
    assert(temp_index);
    assert(list);
    assert(list[taxon_in]);
    assert(list[taxon_in][taxon_out]);
    assert(list[taxon_in][taxon_out]->get_file_content());
    for(uint i = 0; i < index_length; i++) {
      //! Compares the buffer-fragment with the 1d-list:
      assert(arr_index_current[i].get_start_pos() == temp_index[i].get_start_pos());
      assert(arr_index_current[i].get_length() == temp_index[i].get_length());
	      
      //! Compares the real object with the buffer:
      assert(arr_index_current[i].get_start_pos() == list[taxon_in][taxon_out]->get_index_start(i));
      assert(arr_index_current[i].get_length() == list[taxon_in][taxon_out]->get_index_length(i));
    }
#endif
  }
  /**
     @brief Build the datatype in questions.
     @param <temp> In order to underline (specify) that it's the "mpi_transfer_p_rel" that we are interested in calling.
  **/
  void build_datatype_T(p_rel *temp, MPI_Datatype &mpi_type_T) {
    if(false && !temp) ; // In order to hide the fact that we do not use the variable for any other purpose than ensuring that the function we called is the correct one.
    mpi_transfer_p_rel::mpi_build_datatype(mpi_type_T);
  }
  /**
     @brief Build the datatype in questions.
     @param <temp> In order to underline (specify) that it's the "mpi_transfer_p_rel" that we are interested in calling.
  **/
  void build_datatype_T(rel *temp, MPI_Datatype &mpi_type_T) {
    if(false && !temp) ; // In order to hide the fact that we do not use the variable for any other purpose than ensuring that the function we called is the correct one.
    mpi_transfer_rel::mpi_build_datatype(mpi_type_T);
  }


  //! Sets the metadata for the given taxon:
  void set_1d_metadata(file_parse<T> ***list,uint &a_special_overlapping_case, uint *arr_lastPosInOwnData, uint taxon_in, uint taxon_out, uint taxon_length, struct uint_3x_obj &meta_taxa_1d, int id_receive_node) {
    assert(!a_special_overlapping_case);
    a_special_overlapping_case = 0; 
    //! Adjusts: If set to null, is regarded as an "orindary element to be sent as a whole"
    if(arr_lastPosInOwnData && (arr_lastPosInOwnData[get_1d_coordinate(taxon_in, taxon_out, taxon_length)] < list[taxon_in][taxon_out]->get_index_T_end_of_this_object())) { 
      a_special_overlapping_case = arr_lastPosInOwnData[get_1d_coordinate(taxon_in, taxon_out, taxon_length)]; 
    }
    //! Inserts the metadata:
    meta_taxa_1d.total_length_of_index = list[taxon_in][taxon_out]->get_index_length();
    
    assert(0== list[taxon_in][taxon_out]->getLengthOfFile()); // Should not have data in a file when sending
    meta_taxa_1d.total_length_of_rel   = list[taxon_in][taxon_out]->getTotalLengthOfData();
#ifndef NDEBUG
    //! Below variable to validate that receiving of data is correct, and thereby conclude that any difference with our expectations is due to mering of data).
    assert(debug_sum_distances_total_to_send);
    if(id_receive_node != (int)UINT_MAX) {
      debug_sum_distances_total_to_send[id_receive_node] += list[taxon_in][taxon_out]->get_sum_of_pair_distances_for_myrank_in_memory();
    }
#endif    

    //! If back-and-forth-operations are earlier performed, only sends data having the original elements residing.
    if(a_special_overlapping_case>0) {
      meta_taxa_1d.total_length_of_rel = a_special_overlapping_case;
    }    
  }

  //! Inserts the index, ie, data about the buffer:
  void set_1d_index_list(file_parse<T> ***list, uint taxon_in, uint taxon_out, struct uint_3x_obj meta_taxa_1d, uint &cnt_index_inserted, index_t *&arr_index_current, const uint cnt_total_length_of_index_data, const uint a_special_overlapping_case) {
    index_t *temp_index = NULL;
    temp_index = list[taxon_in][taxon_out]->get_index_list();
    assert(temp_index);
    if(temp_index && meta_taxa_1d.total_length_of_index) {
      //! Copies the data, and updates appropriate information:
      assert((cnt_index_inserted + meta_taxa_1d.total_length_of_index) <= cnt_total_length_of_index_data);
      memcpy(arr_index_current, temp_index, sizeof(index_t)*meta_taxa_1d.total_length_of_index);
      verify_is_equal_and_set(arr_index_current, temp_index, meta_taxa_1d.total_length_of_index, taxon_in, taxon_out, list); // Verifies the insertion:
      //! If back-and-forth-operations are earlier performed, only sends data having the original elements residing.
      if(a_special_overlapping_case>0) { // Those values above "threshold" are cleares
	for(uint t = 0; t < meta_taxa_1d.total_length_of_index; t++) {
	  //! If the start position in the buffer is greater tahn before the merging is started, it's a new relation; therefore we do not send it.
	  if(arr_index_current[t].get_start_pos() >= a_special_overlapping_case) {
	    arr_index_current[t] = index_t(); // makes empty.
	  }
	}
      }
      arr_index_current  += meta_taxa_1d.total_length_of_index;  // Updates the pointer.
      cnt_index_inserted += meta_taxa_1d.total_length_of_index; // The size information is updated.
      assert(cnt_index_inserted <= cnt_total_length_of_index_data);            // Validates that our estimate is correct.
    }
  }

  //! Inserts the buffer
  void set_1d_rel_buffer(file_parse<T> ***list, uint taxon_in, uint taxon_out, struct uint_3x_obj meta_taxa_1d, uint &cnt_rel_inserted, T *&arr_rel_current, const uint cnt_total_length_of_rel_data) {
    assert(list[taxon_in][taxon_out]->get_number_of_objects_in_file() == 0); // assumes all data are in memory, before running this operation.
    T *temp_rel = list[taxon_in][taxon_out]->get_file_content();
    assert(temp_rel); // Higher sensitivity for debugging, as we are interested in understanding more.
    if(temp_rel && meta_taxa_1d.total_length_of_rel) {
      //! Copies the data, and updates appropriate information:
      assert((cnt_rel_inserted+meta_taxa_1d.total_length_of_rel) <= cnt_total_length_of_rel_data);            // Validates that our estimate is correct.
      memcpy(arr_rel_current, temp_rel, sizeof(T)*meta_taxa_1d.total_length_of_rel);
      arr_rel_current += meta_taxa_1d.total_length_of_rel;  // Updates the pointer.
      cnt_rel_inserted += meta_taxa_1d.total_length_of_rel; // The size information is updated.
      assert(cnt_rel_inserted <= cnt_total_length_of_rel_data);            // Validates that our estimate is correct.
    }
  }


  //! Transforms the 2d-objects into 1d:
  void transform_2d_into_1d(list_file_parse<T> *this_caller, taxa *listTaxa,
			    file_parse<T> ***list, uint index_start, uint index_end_pluss_one, uint cnt_taxa_pairs, enum collective_type type_method_sending,
			    const uint taxon_length, uint *arr_lastPosInOwnData, uint &elements_used, struct uint_3x_obj *&meta_taxa_1d,
			    int id_receive_node, uint *data_to_send,
			    uint &cnt_rel_inserted, T *&arr_rel_start, uint &cnt_total_length_of_rel_data,
			    uint &cnt_index_inserted, index_t *&arr_index_start, uint &cnt_total_length_of_index_data) {
    //! Allocates, and initiates the lists of metadata
    assert(!meta_taxa_1d);
    meta_taxa_1d = new struct uint_3x_obj[cnt_taxa_pairs];
    log_builder::test_memory_condition_and_if_not_abort(meta_taxa_1d!=NULL, __LINE__, __FILE__, __FUNCTION__);

    //! Allocates, and initiates the index-list:
    assert(!arr_index_start);
    arr_index_start = new index_t[cnt_total_length_of_index_data];
    log_builder::test_memory_condition_and_if_not_abort(arr_index_start!=NULL, __LINE__, __FILE__, __FUNCTION__);
    //    index_t *arr_index_current = arr_index_start;
    //       uint cnt_index_inserted = 0;

    //! Allocates, and initiates the list for the data buffer:
    assert(!arr_rel_start);
    arr_rel_start = new T[cnt_total_length_of_rel_data];
    log_builder::test_memory_condition_and_if_not_abort(arr_rel_start!=NULL, __LINE__, __FILE__, __FUNCTION__);

    index_t *arr_index_current = arr_index_start;
    T *arr_rel_current = arr_rel_start;
    if(list) {
      for(uint i = index_start; i < index_end_pluss_one; i++) {    
	if(list[i]) {
	  for(uint out = index_start; out < index_end_pluss_one; out++) {
	    bool pair_is_to_be_sent = true;
	    if(data_to_send && type_method_sending == mpi_Isend_nonblocking) { // Sends only the specific pairs (of interest):
	      //! Calculates the data the receiver is interested in
	      pair_is_to_be_sent = receive_node_is_interested_in_these_data(id_receive_node, i, out, taxon_length, data_to_send);
	    }

	    if(arr_lastPosInOwnData && (0== arr_lastPosInOwnData[get_1d_coordinate(i, out, taxon_length)])) {pair_is_to_be_sent = false;}
	    int myrank = 0;  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
	    if(pair_is_to_be_sent && list[i][out] && list[i][out]->has_data()) {
	      //! Sets the metadata for the given taxon:
	      uint a_special_overlapping_case = 0;  set_1d_metadata(list, a_special_overlapping_case, arr_lastPosInOwnData, i, out, taxon_length, meta_taxa_1d[elements_used], id_receive_node);

	      //! Inserts the index, ie, data about the buffer:
	      set_1d_index_list(list, i, out, meta_taxa_1d[elements_used], cnt_index_inserted, arr_index_current, cnt_total_length_of_index_data, a_special_overlapping_case);

	      //! Inserts the buffer:
	      set_1d_rel_buffer(list, i, out, meta_taxa_1d[elements_used], cnt_rel_inserted, arr_rel_current, cnt_total_length_of_rel_data);
      
	      //! Sets the taxon-id (due to the spraseness of the data set,
	      //  The id must be included for identification
	      meta_taxa_1d[elements_used].taxon_id = get_1d_coordinate(i, out, taxon_length);
	      elements_used++;
	      assert(elements_used <= cnt_taxa_pairs);
	    }
	  }
	}
      }
#ifndef NDEBUG       
      //! Count the number of relations found in the index-length-list:
      uint temp_cnt_total_length_of_index_data = 0;
      file_parse<T> ***list = this_caller->get_buffer(); // For easy access.
      if(list) {
	for(uint i = 0; i < (uint)taxon_length; i++) {    
	  if(list[i]) {
	    for(uint out = 0; out < (uint)taxon_length; out++) {
	      if(list[i][out] && list[i][out]->has_data()) {
		//! Calculates the data the receiver is interested in
		bool pair_is_to_be_sent = true;
		if(data_to_send && type_method_sending == mpi_Isend_nonblocking) { // Sends only the specific pairs (of interest):
		  pair_is_to_be_sent = receive_node_is_interested_in_these_data(id_receive_node, i, out, taxon_length, data_to_send);
		}		
		if(arr_lastPosInOwnData && (0== arr_lastPosInOwnData[get_1d_coordinate(i, out, taxon_length)])) {pair_is_to_be_sent = false;}

		if(pair_is_to_be_sent) {
		  temp_cnt_total_length_of_index_data += list[i][out]->get_index_length();
		}
	      }
	    }
	  }
	}
      }           
      assert(cnt_index_inserted == temp_cnt_total_length_of_index_data);
      assert(cnt_rel_inserted   == cnt_total_length_of_rel_data);
    
      //
      //! Asserts the merging of the list:

      //! Resets the values:
      arr_index_current = arr_index_start;
      cnt_index_inserted = 0;
      //! Allocates, and initiates the list for the data buffer:
      arr_rel_current = arr_rel_start;
      cnt_rel_inserted = 0;
      //	uint elements_used_debug = 0;
      //! Generates a global object to insert data into:
      char *FILE_BINARY_LOCATION = this_caller->get_FILE_BINARY_LOCATION();
      const bool USE_BEST_BLAST_PAIR_SCORE = this_caller->get_USE_BEST_BLAST_PAIR_SCORE();
      list_file_parse<T> *temp_global          = new list_file_parse<T>(NULL, taxon_length, listTaxa, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);
      list_file_parse<T> *temp_global_onthefly = new list_file_parse<T>(NULL, taxon_length, listTaxa, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);
      assert(temp_global);
      assert(temp_global_onthefly);
      //! Iterates through the 1d-list, covering the different aspects with a variety of tests:
      for(uint id = 0; id<cnt_taxa_pairs; id++) {
	const uint element_id = meta_taxa_1d[id].taxon_id;
	uint taxon_in, taxon_out; get_2d_coordinate(element_id, taxon_in, taxon_out, taxon_length);
	assert(taxon_in < (uint)taxon_length);
	assert(taxon_out < (uint)taxon_length);

	//! If coordinates are correct, the object should be set:
	assert(list[taxon_in][taxon_out]);
    
	//! Asserts the 1d-index-list
	assert(list[taxon_in][taxon_out]->get_file_content());
	for(uint i = 0; i < meta_taxa_1d[id].total_length_of_index; i++) {
	  assert(arr_index_current[i].get_start_pos() == list[taxon_in][taxon_out]->get_index_start(i));
	  assert(arr_index_current[i].get_length() == list[taxon_in][taxon_out]->get_index_length(i));
	}

	//! Asserts the 1d-buffer
	T *temp_rel = list[taxon_in][taxon_out]->get_file_content();	
	assert(temp_rel);
	for(uint i = 0; i < meta_taxa_1d[id].total_length_of_rel; i++) {
	  assert(arr_rel_current[i].equals(temp_rel[i]));
	}
    
	/**
	   Assert a transformation back to the original object:
	   1. Build a new list_parse object.
	   2. Insert the list into it, allocating memory.
	   3. Compare with the internal object, ensuring it's equal.
	   4. Insert it into a global test object
	   5. Inserts it into the global objct, merging "on the fly"
	**/
	file_parse<T> *temp = new file_parse<T>();
	assert(temp);
	file_parse<T> *temp_onthefly = new file_parse<T>();
	assert(temp_onthefly);
	//!       2. Insert the list into it, allocating memory.
	const uint MAX_BUFFER_SIZE = this_caller->get_MAX_BUFFER_SIZE();
	temp->insert_object_using_1d_lists(taxon_in, taxon_out, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE, arr_rel_current, meta_taxa_1d[id].total_length_of_rel, arr_index_current, meta_taxa_1d[id].total_length_of_index, MAX_BUFFER_SIZE); 
	temp_onthefly->insert_object_using_1d_lists(taxon_in, taxon_out, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE, arr_rel_current, meta_taxa_1d[id].total_length_of_rel, arr_index_current, meta_taxa_1d[id].total_length_of_index, MAX_BUFFER_SIZE); 
	//!  3. Compare with the internal object, ensuring it's equal.
	assert(temp->is_equal(*list[taxon_in][taxon_out], true));
	assert(temp_onthefly->is_equal(*list[taxon_in][taxon_out], true));

	//!  4. Insert it into a global test object
	temp_global->set_element_unsafe(taxon_in, taxon_out, temp);
    
	//!  5. Inserts it into the global objct, merging "on the fly"      
	temp_global_onthefly->merge_data(temp_onthefly, taxon_in, taxon_out);
	assert(!temp_onthefly);

	//! Updates the pointers, making them consisten for the next run:
	arr_index_current  += meta_taxa_1d[id].total_length_of_index;  // Updates the pointer.
	arr_rel_current    += meta_taxa_1d[id].total_length_of_rel;  // Updates the pointer.
      }
      if(type_method_sending != mpi_Isend_nonblocking) {
	//! Verifies that the global test object is equal to this, ie, that our 1d-object correspons to the real object:
	assert(temp_global->is_equal(this_caller, true));
	assert(temp_global_onthefly->is_equal(this_caller, true));
      }
      //! Deallocates the data:
      list_file_parse<T>::close(temp_global, false);
      list_file_parse<T>::close(temp_global_onthefly, false);
#endif
    }
  }

  /**
     @brief Transfers the mpi data residing in a list_file_parse object.
     @remarks A deadlock is caused if this function is not composed with a "reciprocal" receive fuctntion-call.
  **/
  void mpi_transfer_data(struct uint_3x_obj obj_info,MPI_Comm comm,int myrank, int id_receive_node, int number_of_nodes, uint index_start, uint index_end_pluss_one, list_file_parse<T> *this_caller, enum collective_type type_method_sending) {
    uint *temp = NULL;  uint *temp_2 = NULL;
    //! The container- and datatype for the object holding information about the lengths to be transferred:
    struct uint_3x_obj *meta_taxa_1d = NULL;
    //! The container and datatype for the index_t object:
    index_t *arr_index_start = NULL;
    //! The container and datatype for the T object:
    T *arr_rel_start = NULL;
    MPI_Request meta_send_req, index_send_req, T_send_req; // These variables are provided, but not used in this context.
    mpi_transfer_data(obj_info, comm, myrank, id_receive_node, number_of_nodes, index_start, index_end_pluss_one, this_caller, type_method_sending, temp, temp_2, meta_taxa_1d, arr_index_start, arr_rel_start, 
		      meta_send_req, index_send_req, T_send_req
);
  }

  /**
     @brief Transfers the mpi data residing in a list_file_parse object.
     @remarks A deadlock is caused if this function is not composed with a "reciprocal" receive fuctntion-call.
  **/
  void mpi_transfer_data(struct uint_3x_obj obj_info,MPI_Comm comm,int myrank, int id_receive_node, int number_of_nodes, uint index_start, uint index_end_pluss_one, list_file_parse<T> *this_caller, enum collective_type type_method_sending, uint *data_to_send, uint *arr_lastPosInOwnData,
			 struct uint_3x_obj *&meta_taxa_1d, index_t *&arr_index_start, T *&arr_rel_start,
			 MPI_Request &meta_send_req, MPI_Request &index_send_req, MPI_Request &T_send_req
			 ) {
    assert(this_caller);
    assert(number_of_nodes);
    assert(index_end_pluss_one);
    //! Get some variables from the caller for easy access:
    file_parse<T> ***list = this_caller->get_buffer();

    //! Sets the length of the lists, avoiding to resize:
     uint cnt_taxa_pairs = obj_info.taxon_id, cnt_total_length_of_rel_data = obj_info.total_length_of_rel, cnt_total_length_of_index_data = obj_info.total_length_of_index;

    //! Sends the length of the lists to the node(s) redeving data:
    if(type_method_sending != mpi_Isend_nonblocking) {
      const int root_id = id_receive_node;
      //! Sends the metadata, either to the root, of from the root to all the other processes:
      if(myrank != 0) MPI_Send(&obj_info, 1, mpi_type_uint_3x, root_id, object_info, comm);
      else           MPI_Bcast(&obj_info, 1, mpi_type_uint_3x, root_id, comm);
    } // else an other sending procedure is used.

    
    if(true) {
      //! Transforms the 2d-objects into 1d:
      uint cnt_index_inserted = 0, cnt_rel_inserted = 0, elements_used = 0;
      assert(!meta_taxa_1d);
      assert(!arr_index_start);
      assert(!arr_rel_start);
      transform_2d_into_1d(this_caller, listTaxa, list, index_start, index_end_pluss_one, cnt_taxa_pairs, type_method_sending,
			   taxon_length, arr_lastPosInOwnData, elements_used, meta_taxa_1d,
			   id_receive_node, data_to_send,
			   // Data for the T type:
			   cnt_rel_inserted, arr_rel_start, cnt_total_length_of_rel_data,
			 // Data for the index_t type:
			   cnt_index_inserted, arr_index_start, cnt_total_length_of_index_data);
      assert(elements_used == cnt_taxa_pairs); // The number of inserted elements should correspond to the reserved length.
    }

    //
    //! The mpi sending procedure, summing the values using the same pointer for the task:

    //! If non-blocking send is used, store the reqeust info in them:

#ifndef NDEBUG
    //! Sends metadata to be verified by the receiver:
    if(type_method_sending == mpi_Isend_nonblocking) {
      assert(debug_sum_distances_total_to_send);
      assert(debug_sum_distances_total_to_send[id_receive_node] != 0);
      // As this is the smallest element to send, we assume it has arrived when the longer data has arrived:
      MPI_Request temp_req;
      MPI_Isend(&debug_sum_distances_total_to_send[id_receive_node], 1, MPI_FLOAT, id_receive_node, debug_sum_distance, comm, &temp_req);
    }
#endif

    //! Sends the length of the subsections the "main list" consists of, to the node(s) redeving data:
    if(type_method_sending != mpi_Isend_nonblocking) {
      const int root_id = id_receive_node;
      if(myrank != 0) MPI_Send(meta_taxa_1d, cnt_taxa_pairs, mpi_type_uint_3x, root_id, object_info_in_list, comm);
      else           MPI_Bcast(meta_taxa_1d, cnt_taxa_pairs, mpi_type_uint_3x, root_id, comm);
    } else MPI_Isend(meta_taxa_1d, cnt_taxa_pairs, mpi_type_uint_3x, id_receive_node, object_info_in_list, comm, &meta_send_req);

    //! Sends the index_list:
    if(type_method_sending != mpi_Isend_nonblocking) {
      const int root_id = id_receive_node;
      if(myrank != 0) MPI_Send(arr_index_start, cnt_total_length_of_index_data, mpi_type_uint_2x, root_id, index_list, comm);
      else           MPI_Bcast(arr_index_start, cnt_total_length_of_index_data, mpi_type_uint_2x, root_id, comm);
    } else MPI_Isend(arr_index_start, cnt_total_length_of_index_data, mpi_type_uint_2x, id_receive_node, index_list, comm, &index_send_req);

    //! Sends the data of type T
    if(type_method_sending != mpi_Isend_nonblocking) {
      const int root_id = id_receive_node;
      if(myrank != 0) MPI_Send(arr_rel_start, cnt_total_length_of_rel_data, mpi_type_T, root_id, T_list, comm);
      else           MPI_Bcast(arr_rel_start, cnt_total_length_of_rel_data, mpi_type_T, root_id, comm);
    } else MPI_Isend(arr_rel_start, cnt_total_length_of_rel_data, mpi_type_T, id_receive_node, T_list, comm, &T_send_req);

    //! Deallocates the objects at the end:
    if(type_method_sending != mpi_Isend_nonblocking) {
      if(meta_taxa_1d) {delete [] meta_taxa_1d, meta_taxa_1d = NULL;}
      if(arr_index_start) {delete [] arr_index_start, arr_index_start = NULL;}
      if(arr_rel_start) {delete [] arr_rel_start, arr_rel_start = NULL;}    
    }
  }

  //! Deallocates the remaning data containers used for sendning:
  void deallocate_data_containers_used_for_sending(int myrank, int number_of_nodes,
						   struct uint_3x_obj **nodelist_meta_taxa_1d, index_t **nodelist_arr_index_start,
						   T **nodelist_arr_rel_start, MPI_Request *nodelist_meta_send_req,
						   MPI_Request *nodelist_index_send_req, MPI_Request *nodelist_T_send_req) {
    assert(number_of_nodes);
    assert(nodelist_meta_taxa_1d);
    assert(nodelist_arr_index_start);
    assert(nodelist_arr_rel_start);
    assert(nodelist_meta_send_req);
    assert(nodelist_index_send_req);
    assert(nodelist_T_send_req);
    for(int id_receiver = 0; id_receiver < number_of_nodes; id_receiver++) {
      if(id_receiver !=myrank) {
	if(nodelist_arr_rel_start[id_receiver] != NULL) {
	  //! The data is set, implying this node sendt some data, ie, deallocating must be performed
	  //! - Deallocates for the meta-data when data is sent:	  
	  MPI_Wait(&nodelist_meta_send_req[id_receiver], MPI_STATUS_IGNORE);
	  if(nodelist_meta_taxa_1d[id_receiver]) {delete [] nodelist_meta_taxa_1d[id_receiver], nodelist_meta_taxa_1d[id_receiver] = NULL;}
	  //! - Deallocates for the index-data when data is sent:	  
	  MPI_Wait(&nodelist_index_send_req[id_receiver], MPI_STATUS_IGNORE);
	  if(nodelist_arr_index_start[id_receiver]) {delete [] nodelist_arr_index_start[id_receiver], nodelist_arr_index_start[id_receiver] = NULL;}
	  //! - Deallocates for the rel-data when data is sent:	  
	  MPI_Wait(&nodelist_T_send_req[id_receiver], MPI_STATUS_IGNORE);
	  if(nodelist_arr_rel_start[id_receiver]) {delete [] nodelist_arr_rel_start[id_receiver], nodelist_arr_rel_start[id_receiver] = NULL;}	  
#ifndef NDEBUG
	  //! Confirms that we are at the end.
	  //! Note: Also works as a show-case for the 'MPI_Test(..)' functionality:
	  int flag = 0; MPI_Test(&nodelist_meta_send_req[id_receiver], &flag, MPI_STATUS_IGNORE);
	  assert(flag == true);
	  flag = 0; MPI_Test(&nodelist_index_send_req[id_receiver], &flag, MPI_STATUS_IGNORE);
	  assert(flag == true);
	  flag = 0; MPI_Test(&nodelist_T_send_req[id_receiver], &flag, MPI_STATUS_IGNORE);
	  assert(flag == true);
#endif
	}	
      }
    }  
  }

  /**
     @brief Does not wait: Only deallocates the data containers if data is received.
     @remarks
     # Purpose is to reduce the memory consumption without introducing lag (waits)
     # The order of the deallocating goes from the smallest chunk sent, towards the biggest (The *T object).
     # Only test if data is set for the biggest container
  **/
  void non_final_deallocate_data_containers_used_for_sending(int myrank, int number_of_nodes,
						   struct uint_3x_obj **nodelist_meta_taxa_1d, index_t **nodelist_arr_index_start,
						   T **nodelist_arr_rel_start, MPI_Request *nodelist_meta_send_req,
						   MPI_Request *nodelist_index_send_req, MPI_Request *nodelist_T_send_req) {
    assert(number_of_nodes);
    assert(nodelist_meta_taxa_1d);
    assert(nodelist_arr_index_start);
    assert(nodelist_arr_rel_start);
    assert(nodelist_meta_send_req);
    assert(nodelist_index_send_req);
    assert(nodelist_T_send_req);
    for(int id_receiver = 0; id_receiver < number_of_nodes; id_receiver++) {
      if(id_receiver !=myrank) {
	if(nodelist_arr_rel_start[id_receiver] != NULL) {
	  int flag = 0; MPI_Test(&nodelist_meta_send_req[id_receiver], &flag, MPI_STATUS_IGNORE);
	  if(flag) {
	    //! - Deallocates for the meta-data when data is sent:	  
	    if(nodelist_meta_taxa_1d[id_receiver]) {delete [] nodelist_meta_taxa_1d[id_receiver], nodelist_meta_taxa_1d[id_receiver] = NULL;}
	    flag = 0; MPI_Test(&nodelist_index_send_req[id_receiver], &flag, MPI_STATUS_IGNORE);
	    if(flag) {
	      //! - Deallocates for the index-data when data is sent:	  
	      if(nodelist_arr_index_start[id_receiver]) {delete [] nodelist_arr_index_start[id_receiver], nodelist_arr_index_start[id_receiver] = NULL;}	      
	      flag = 0; MPI_Test(&nodelist_T_send_req[id_receiver], &flag, MPI_STATUS_IGNORE);
	      if(flag) {
		//! - Deallocates for the rel-data when data is sent:	  
		if(nodelist_arr_rel_start[id_receiver]) {delete [] nodelist_arr_rel_start[id_receiver], nodelist_arr_rel_start[id_receiver] = NULL;}
	      }
	    }
	  }
	}	
      }
    }  
  }
  /**
     @brief Merges the data into the list_file_parse object
     @remarks
     # Expands the lis
     # Note: The functions used are verified thoroughly in the "transfer" function, found in their class
  **/
  void merge_data(struct uint_3x_obj *meta_taxa_1d, index_t *arr_index_start, T *arr_rel_start, struct uint_3x_obj obj_info) {
    //! Assumes that when this function is called, the following variables are set.
    assert(meta_taxa_1d);
    assert(arr_index_start);
    assert(arr_rel_start);
    assert(obj_info.is_set());
    //! Initiates local variables, making the code similar to the code send sending the objects:
    const int taxon_length = this_caller->get_taxon_length();
    char *FILE_BINARY_LOCATION = this_caller->get_FILE_BINARY_LOCATION();
    const bool USE_BEST_BLAST_PAIR_SCORE = this_caller->get_USE_BEST_BLAST_PAIR_SCORE();
    const uint cnt_taxa_pairs = obj_info.taxon_id; // 'taxon_id' represents there the total number of taxon-pairs.
    const uint cnt_total_length_of_index_data = obj_info.total_length_of_index;
    const uint cnt_total_length_of_rel_data   = obj_info.total_length_of_rel;
    //! Resets the values, and initiates:
    index_t *arr_index_current = arr_index_start;
    T *arr_rel_current = arr_rel_start;
#ifndef NDEBUG
    float debug_sum_total = 0;
    const float sum_before_merging = this_caller->get_sum_of_pair_distances_in_memory();
#endif
    //! Merges each of the taxon pairs given:
    for(uint id = 0; id<cnt_taxa_pairs; id++) {
      //! Inititiate values
      const uint element_id = meta_taxa_1d[id].taxon_id;
      uint taxon_in, taxon_out; get_2d_coordinate(element_id, taxon_in, taxon_out, taxon_length);
      assert(taxon_in < (uint)taxon_length);
      assert(taxon_out < (uint)taxon_length);	    
      //! Build a new list_parse object.
      file_parse<T> *temp = new file_parse<T>();	    
      assert(temp);

      //!  Insert the list into it, allocating memory.
      const uint MAX_BUFFER_SIZE = this_caller->get_MAX_BUFFER_SIZE();
      //printf("before[%u][%u] (with buffer=%u) at line %d in file %s\n", taxon_in, taxon_out, meta_taxa_1d[id].total_length_of_rel, __LINE__, __FILE__); // TODO: remove this line
      //index::print_list(arr_index_current, meta_taxa_1d[id].total_length_of_index); // TODO: remove this line
      //T::print_buffer(arr_rel_current, meta_taxa_1d[id].total_length_of_rel); // TODO: remove this line      
      temp->insert_object_using_1d_lists(taxon_in, taxon_out, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE, arr_rel_current, meta_taxa_1d[id].total_length_of_rel, arr_index_current, meta_taxa_1d[id].total_length_of_index, MAX_BUFFER_SIZE); 

#ifndef NDEBUG
      //! Validates the inserted sum of distances corresponds to our calculations:
      const float expected_sum_of_distances = temp->get_sum_of_pair_distances_for_myrank_in_memory();
      float sum = 0;
      for(uint i = 0; i < meta_taxa_1d[id].total_length_of_rel; i++) {
	sum += arr_rel_current[i].get_distance();
      }
      const float sum_local_tmp = relations_list<T>::get_sum_of_pair_distances_in_memory(arr_rel_current, meta_taxa_1d[id].total_length_of_rel);
      assert(sum == sum_local_tmp);
      assert(sum == expected_sum_of_distances);
      debug_sum_total += sum;
#endif
      assert(temp->buffer_has_data_set()); // Ensures that data is both set, and consistent (both list and size is set).
      this_caller->merge_data(temp, taxon_in, taxon_out);
      assert(!temp); // Ensures it is freed.
      assert(this_caller->buffer_is_not_null(taxon_in, taxon_out)); // Ensures that data is both set, and consistent (both list and size is set).
#ifndef NDEBUG
      //! Validates that the structure is correctly updated:
      const float expected_sum_of_distances_after_merge = this_caller->get_sum_of_pair_distances_in_memory();
      const float inserted_sum = sum_before_merging + debug_sum_total;
      log_builder::compare_floats(expected_sum_of_distances_after_merge, inserted_sum, meta_taxa_1d[id].total_length_of_rel, __LINE__, __FILE__, __FUNCTION__);
#endif

      //! Updates the pointers, making them consisten for the next run:
      arr_index_current  += meta_taxa_1d[id].total_length_of_index;  // Updates the pointer.
      assert(arr_index_current <= (arr_index_start+cnt_total_length_of_index_data));
      arr_rel_current    += meta_taxa_1d[id].total_length_of_rel;  // Updates the pointer.
      assert(arr_rel_current <= (arr_rel_start+cnt_total_length_of_rel_data));
    }
  }
 /**
     @brief Performs the receive operation, both for the collective operation(myrank=0) and the "spreading" operation (myrank!=0).
  **/
  void mpi_receive_data_and_merge_it_into_this(MPI_Comm comm, int myrank, int number_of_nodes, uint index_start, uint index_end_pluss_one, list_file_parse<T> *this_caller) {
    const int root_id = 0; // Expessive identification of the root:
    if(myrank != 0) {
      //! Receives the metadata, either to the root, of from the root to all the other processes:
      struct uint_3x_obj obj_info = uint_3x_obj_t();
      MPI_Bcast(&obj_info, 1, mpi_type_uint_3x, root_id, comm);
      //! Initiates local variables, making the code similar to the code sed sending the objects:
      const uint cnt_taxa_pairs = obj_info.taxon_id; // 'taxon_id' represents there the total number of taxon-pairs.         
      const uint cnt_total_length_of_index_data = obj_info.total_length_of_index;
      const uint cnt_total_length_of_rel_data   = obj_info.total_length_of_rel;

      
      //! Allocates, and initiates the list for the lists (buffers):
      struct uint_3x_obj *meta_taxa_1d = new struct uint_3x_obj[cnt_taxa_pairs];
      index_t *arr_index_start = new index_t[cnt_total_length_of_index_data];
      T *arr_rel_start = new T[cnt_total_length_of_rel_data];
      assert(meta_taxa_1d);
      assert(arr_index_start);
      assert(arr_rel_start);
      //! Receives the length of the subsections the "main list" consists of, to the node(s) redeving data:      
      MPI_Bcast(meta_taxa_1d, cnt_taxa_pairs, mpi_type_uint_3x, root_id, comm);
      //! Receives the index_list:
      MPI_Bcast(arr_index_start, cnt_total_length_of_index_data, mpi_type_uint_2x, root_id, comm);
      //! Receives the data of type T
      MPI_Bcast(arr_rel_start, cnt_total_length_of_rel_data, mpi_type_T, root_id, comm);
      
      //
      //! Expands the list, by first deleting the internal data in "this_caller" object!
      // TODO: This operation might not be the best, as some relations are first deleted, and then re-inserted, therefore
      //       it might be worthwile considering an alternative approach.
      this_caller->free_file_parse_memory(true, index_start, index_end_pluss_one); 
      merge_data(meta_taxa_1d, arr_index_start, arr_rel_start, obj_info);

      //! Deallocates the objects at the end:
      if(meta_taxa_1d) {delete [] meta_taxa_1d, meta_taxa_1d = NULL;}
      if(arr_index_start) {delete [] arr_index_start, arr_index_start = NULL;}
      if(arr_rel_start) {delete [] arr_rel_start, arr_rel_start = NULL;}
    } else {
      for(uint sender_id = 0; sender_id < (uint)number_of_nodes; sender_id++) {
	if(sender_id != (uint)root_id) {
	  //! Receives the metadata, either to the root, of from the root to all the other processes:
	  struct uint_3x_obj obj_info = uint_3x_obj_t();
	  MPI_Recv(&obj_info, 1, mpi_type_uint_3x, sender_id, object_info, comm, MPI_STATUS_IGNORE);

	  //! Initiates local variables, making the code similar to the code sed sending the objects:
	  const uint cnt_taxa_pairs = obj_info.taxon_id; // 'taxon_id' represents there the total number of taxon-pairs.
	  const uint cnt_total_length_of_index_data = obj_info.total_length_of_index;
	  const uint cnt_total_length_of_rel_data   = obj_info.total_length_of_rel;

	  //! Allocates, and initiates the list for the lists (buffers), using the "lengths" received:
	  struct uint_3x_obj *meta_taxa_1d = new struct uint_3x_obj[cnt_taxa_pairs]; 
	  index_t *arr_index_start = new index_t[cnt_total_length_of_index_data];
	  T *arr_rel_start = new T[cnt_total_length_of_rel_data];
	  assert(meta_taxa_1d);
	  assert(arr_index_start);
	  assert(arr_rel_start);
	  //! Receives the length of the subsections the "main list" consists of, to the node(s) redeving data:      
	  MPI_Recv(meta_taxa_1d, cnt_taxa_pairs, mpi_type_uint_3x, sender_id, object_info_in_list, comm, MPI_STATUS_IGNORE);
	  //! Receives the index_list:
	  MPI_Recv(arr_index_start, cnt_total_length_of_index_data, mpi_type_uint_2x, sender_id, index_list, comm, MPI_STATUS_IGNORE);
	  //! Receives the data of type T
	  MPI_Recv(arr_rel_start, cnt_total_length_of_rel_data, mpi_type_T, sender_id, T_list, comm, MPI_STATUS_IGNORE);

	  //! Expands the list
	  merge_data(meta_taxa_1d, arr_index_start, arr_rel_start, obj_info);

	  //! Deallocates the objects at the end:
	  if(meta_taxa_1d) {delete [] meta_taxa_1d, meta_taxa_1d = NULL;}
	  if(arr_index_start) {delete [] arr_index_start, arr_index_start = NULL;}
	  if(arr_rel_start) {delete [] arr_rel_start, arr_rel_start = NULL;}
	}
      }
    }
  }
 public:
  /**
     @brief Sends- and retrieves the data, using functions defined in the MPI library.
  **/
  void mpi_make_data_consistent_accross_nodes(MPI_Comm comm, int myrank, int number_of_nodes, uint index_start, uint index_end_pluss_one, list_file_parse<T> *this_caller) {
    assert(this_caller);
    assert(number_of_nodes>0);
    assert(index_end_pluss_one);
    //! Collects all the data into myrank=0:
    const int root_id = 0;
    const uint taxon_length = this_caller->get_taxon_length();

    bool *list_of_myrank_taxa_responsilibties = NULL;
      //! Identifies taxa this object is responsible for, ie, the node this object belongs to; set from "mpi_transfer_list_file_parse" object:
    if(!initiate_taxa_responsibilities(list_of_myrank_taxa_responsilibties)) return;
    
    //! Each node generates a list-overview of the number of taxon-pairs for each taxa,
    //! followed by an an all-reduce operation, summing the total number of data for each taxon.        
    //! Note: 'list_of_data_to_send' holds the complete list of the length for each taxon.
    uint *list_of_data_to_send = NULL, *sum_for_all_nodes_length_of_taxa = NULL;
    tx_rx_complete_lengths_of_taxa(list_of_data_to_send, sum_for_all_nodes_length_of_taxa);
    
    //! Build the lists: The node with the highest frequency for a taxon, gets the responsibility; inititates.
    uint *responsibilities = NULL, *nodes_having_taxon = NULL, *maximum_lengths = NULL;
    set_myranks_taxa_responsibility(responsibilities,nodes_having_taxon, maximum_lengths, list_of_data_to_send, list_of_myrank_taxa_responsilibties);
    
    if(myrank != 0) {
      //!  Does preprocessing, to know the length of the lists, avoiding to resize:
      uint cnt_taxa_pairs = 0, cnt_total_length_of_rel_data = 0, cnt_total_length_of_index_data = 0;
      get_metadata_information(index_start, index_end_pluss_one, cnt_taxa_pairs, cnt_total_length_of_rel_data, cnt_total_length_of_index_data, this_caller);
      struct uint_3x_obj obj_info = uint_3x_obj_t(cnt_taxa_pairs, cnt_total_length_of_rel_data, cnt_total_length_of_index_data);
      //! Calls the transfer routine:
      mpi_transfer_data(obj_info, comm, myrank, root_id, number_of_nodes, index_start, index_end_pluss_one, this_caller, mpi_send_blocking);
    } else mpi_receive_data_and_merge_it_into_this(comm, myrank, number_of_nodes, index_start, index_end_pluss_one, this_caller);
    
    //! Transfers the data to all the nodes, doing the sending in "opposite" manner:
    if(myrank == 0) {
      //! Does preprocessing, to know the length of the lists, avoiding to resize:
      uint cnt_taxa_pairs = 0, cnt_total_length_of_rel_data = 0, cnt_total_length_of_index_data = 0;
      get_metadata_information(index_start, index_end_pluss_one, cnt_taxa_pairs, cnt_total_length_of_rel_data, cnt_total_length_of_index_data, this_caller);
      struct uint_3x_obj obj_info = uint_3x_obj_t(cnt_taxa_pairs, cnt_total_length_of_rel_data, cnt_total_length_of_index_data);
      //! Calls the transfer routine:
      mpi_transfer_data(obj_info, comm, myrank, root_id, number_of_nodes, index_start, index_end_pluss_one, this_caller, mpi_bcast);
    } else mpi_receive_data_and_merge_it_into_this(comm, myrank, number_of_nodes, index_start, index_end_pluss_one, this_caller);

    //!    Update the scheduling algorithm with the list (using macro variables making this a specific option if MPI is compiled into it).
    this_caller->set_list_of_nodes_taxa_responsibilities(list_of_myrank_taxa_responsilibties, taxon_length);
    assert(list_of_myrank_taxa_responsilibties == NULL);

      if(list_of_data_to_send) {delete [] list_of_data_to_send; list_of_data_to_send = NULL;}
      if(responsibilities) {delete [] responsibilities; responsibilities = NULL; nodes_having_taxon = NULL; maximum_lengths = NULL;}
      //    if(arr_lastPosInOwnData) {delete [] arr_lastPosInOwnData; arr_lastPosInOwnData = NULL;}
      //    if(nodelist_meta_taxa_1d) {delete [] nodelist_meta_taxa_1d, nodelist_meta_taxa_1d = NULL;}
      //    if(nodelist_arr_index_start) {delete [] nodelist_arr_index_start, nodelist_arr_index_start = NULL;}
      //    if(nodelist_arr_rel_start) {delete [] nodelist_arr_rel_start, nodelist_arr_rel_start = NULL;}	  
      //#ifndef NDEBUG
	if(sum_for_all_nodes_length_of_taxa) {delete [] sum_for_all_nodes_length_of_taxa; sum_for_all_nodes_length_of_taxa = NULL;}
  }

 private:
  //! @return false if furhter processing is not to be taken place
  bool initiate_taxa_responsibilities(bool *&list_of_myrank_taxa_responsilibties) {
    assert(taxon_length);
    assert(!list_of_myrank_taxa_responsilibties);
    list_of_myrank_taxa_responsilibties = new bool[taxon_length];    
    log_builder::test_memory_condition_and_if_not_abort(list_of_myrank_taxa_responsilibties!=NULL, __LINE__, __FILE__, __FUNCTION__);    
    memset(list_of_myrank_taxa_responsilibties, 0, sizeof(bool)*taxon_length);

    //! Only sends data to- and from if there are more than one node involved:
    if(number_of_nodes == 1) {
      for(uint i = 0; i < taxon_length; i++) list_of_myrank_taxa_responsilibties[i] = true;
      //!    Update the scheduling algorithm with the list (using macro variables making this a specific option if MPI is compiled into it).
      this_caller->set_list_of_nodes_taxa_responsibilities(list_of_myrank_taxa_responsilibties, taxon_length);
      assert(list_of_myrank_taxa_responsilibties == NULL);
      return false;
    } else return true;
  }

  /**
     Each node generates a list-overview of the number of taxon-pairs for each taxa,
     followed by an an all-reduce operation, summing the total number of data for each taxon.
  **/
  void tx_rx_complete_lengths_of_taxa(uint *&list_of_data_to_send, uint *&sum_for_all_nodes_length_of_taxa) {
    assert(number_of_nodes);
    assert(taxon_length);
    assert(!list_of_data_to_send);

    const uint total_size = number_of_nodes * taxon_length;
    list_of_data_to_send = new uint[total_size];
    assert(list_of_data_to_send);
    log_builder::test_memory_condition_and_if_not_abort(list_of_data_to_send!=NULL, __LINE__, __FILE__, __FUNCTION__);
    memset(list_of_data_to_send, 0, sizeof(uint)*total_size);
    //    const uint myrank_start_pos = myrank * taxon_length;
    uint sum_inserted = 0;
#ifndef NDEBUG
    //! To verify that the transfrring is correct:
    sum_for_all_nodes_length_of_taxa = new uint[taxon_length];
    assert(sum_for_all_nodes_length_of_taxa);
    memset(sum_for_all_nodes_length_of_taxa, 0, sizeof(uint)*taxon_length);
#endif

    for(uint i = 0; i < taxon_length; i++) {
      const uint myrank_start_pos = myrank * taxon_length;
      const uint index_this = i+myrank_start_pos;
      assert(index_this == get_1d_node_taxon_coordinate(myrank, i, taxon_length));
      list_of_data_to_send[index_this] = this_caller->get_total_number_of_pairs_for_taxon(i);
      sum_inserted += this_caller->get_total_number_of_pairs_for_taxon(i);
#ifndef NDEBUG
      sum_for_all_nodes_length_of_taxa[i] = this_caller->get_total_number_of_pairs_for_taxon(i);
#endif
    }
    assert(sum_inserted == this_caller->getTotalLengthOfData());
    
    //! The mpi sending procedure, summing the values using the same pointer for the task:
    // Note: The MPI_Op type to use does not matter, as the lists not set should be disjoint.
    MPI_Allreduce(MPI_IN_PLACE, list_of_data_to_send, total_size, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

#ifndef NDEBUG
    //! Update the container, to hold the information about the total lengths for the taxa:
    MPI_Allreduce(MPI_IN_PLACE, sum_for_all_nodes_length_of_taxa, taxon_length, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#endif
  }


  //! Build the lists: The node with the highest frequency for a taxon, gets the responsibility; inititates.
  void set_myranks_taxa_responsibility(uint *&responsibilities, uint *&nodes_having_taxon, uint *&maximum_lengths, uint *list_of_data_to_send, bool *list_of_myrank_taxa_responsilibties) {
    //! The first taxon_length elements holds the id, the second block of taxon_length holds
    //! the number of nodes having data for this taxon, and the last holds the maximum lengths.
    const uint responsibilities_length = taxon_length * 3; 
    responsibilities = new uint[responsibilities_length+1];// '+1' to use the last bit for testing.
    assert(responsibilities);
    log_builder::test_memory_condition_and_if_not_abort(responsibilities!=NULL, __LINE__, __FILE__, __FUNCTION__);
    memset(responsibilities, 0, (sizeof(uint)*responsibilities_length));
    responsibilities[responsibilities_length] = UINT_MAX; // In order to verify that our pointer-handling is correct.
    nodes_having_taxon = responsibilities + taxon_length;
    maximum_lengths = nodes_having_taxon + taxon_length;
    assert(maximum_lengths[taxon_length] == UINT_MAX);
   
    //! Uses the data container to merge:
    //! NOTE: Performed seperately (the same operation) for all the nodes:
    for(uint id = 0; id < (uint)number_of_nodes; id++) {
      const uint id_start_pos = id * taxon_length;	
      for(uint i = 0; i < taxon_length; i++) {
	if(list_of_data_to_send[id_start_pos+i]) nodes_having_taxon[i]++;
	if(list_of_data_to_send[id_start_pos+i] > maximum_lengths[i]) {
	  // Takes the responsilibty, and update the length:
	  responsibilities[i]   = id;
	  maximum_lengths[i] = list_of_data_to_send[id_start_pos+i];
	  if(id == (uint)myrank) {list_of_myrank_taxa_responsilibties[i] = true;
	  } else            list_of_myrank_taxa_responsilibties[i] = false;
	}
      }
    }
  }

  //! Prints the responsibliites for 'this node':
  void print_myrank_responsiblities(uint taxon_length, uint *responsibilities, uint *maximum_lengths) {
    assert(taxon_length);
    assert(responsibilities);
    assert(maximum_lengths);
    for(uint i = 0; i < taxon_length; i++) {
      printf("[myrank=%d]\t-\tnode_id[%d] is responsible for taxon[%u] with length(%u), at line %d in file %s\n", myrank, responsibilities[i], i, maximum_lengths[i], __LINE__, __FILE__);
    }
  }
  /**
     @brief Sends the metadata to the correct node, who receives it.
     @remarks:
     # Sends information about the memory needed to be allocated before receiving the data,.
     # Regarded as the preprocessing, ie, to know the length of the lists, avoiding to resize.
  **/
  void tx_rx_taxon_specific_lengths(struct uint_3x_obj *sender_obj_info, struct uint_3x_obj *receiver_obj_info, 
				    //struct uint_3x_obj **debug_sender_obj_info, struct uint_3x_obj **debug_receiver_obj_info, 
				    uint *sender_nodes_having_data_for, uint **&arr_lastPosInOwnData, uint *list_of_data_to_send, uint *sum_for_all_nodes_length_of_taxa, uint *cnt_received_rel_not_reciprocal, MPI_Request *receive_status, MPI_Request *receive_status_debug) {
    assert(taxon_length);
    assert(sender_obj_info);
    assert(receiver_obj_info);
    assert(sender_nodes_having_data_for);
    assert(receive_status);
    //! Initializes
    bool arr_lastPos_was_set_before_function_call = false;
    
    if(!arr_lastPosInOwnData) {
      arr_lastPosInOwnData = new uint*[number_of_nodes];
      log_builder::test_memory_condition_and_if_not_abort(arr_lastPosInOwnData!=NULL, __LINE__, __FILE__, __FUNCTION__);
    } else arr_lastPos_was_set_before_function_call = true;
    //! Initialises the containers:
    for(uint id_receiver = 0; id_receiver < (uint)number_of_nodes; id_receiver++) {
      sender_nodes_having_data_for[id_receiver] = 0;
      receiver_obj_info[id_receiver] = uint_3x_obj_t();
      if(!arr_lastPos_was_set_before_function_call) {
	//! Intialises it:
	arr_lastPosInOwnData[id_receiver] = NULL;
	cnt_received_rel_not_reciprocal[id_receiver] = 0;
	sender_obj_info[id_receiver] = uint_3x_obj_t();
      } else {
	assert(arr_lastPosInOwnData[id_receiver]); // Should be set
      }
    }
    
    //! Iterate thorugh the nodes, coubnting the number of nodes the receiver will wait for, before the operation is completed:
    for(uint id_receiver = 0; id_receiver < (uint)number_of_nodes; id_receiver++) {
      if(id_receiver != (uint)myrank) { // Sends the metadata:
	uint cnt_p_not_reciprocal = 0; // Value used in debug mode for the receiver, ensuring that the new vount is correct.
	assert(!arr_lastPosInOwnData[id_receiver]);
	assert(!cnt_received_rel_not_reciprocal[myrank]);	  
	arr_lastPosInOwnData[id_receiver] = get_metadata_information(id_receiver, list_of_data_to_send, taxon_length, number_of_nodes, sender_obj_info[id_receiver], this_caller, cnt_p_not_reciprocal);

	//! Sends the length of the lists to the node receiving data:	
	MPI_Request send_req; MPI_Isend(&sender_obj_info[id_receiver], 1, mpi_type_uint_3x, id_receiver, object_info, MPI_COMM_WORLD, &send_req);
#ifndef NDEBUG
	MPI_Request send_req_debug; MPI_Isend(&cnt_p_not_reciprocal, 1, MPI_UNSIGNED, id_receiver, 1, MPI_COMM_WORLD, &send_req_debug);
#endif
      }
    }

    //! Receives the metadata from the other processes:
    for(int sender_id = 0; sender_id < number_of_nodes; sender_id++) {
      if(sender_id != myrank) { // Does not communicate with "my self"
	MPI_Irecv(&receiver_obj_info[sender_id], 1, mpi_type_uint_3x, sender_id, object_info, MPI_COMM_WORLD, &receive_status[sender_id]);
#ifndef NDEBUG
	MPI_Irecv(&cnt_received_rel_not_reciprocal[sender_id],1,MPI_UNSIGNED,sender_id,1, MPI_COMM_WORLD, &receive_status_debug[sender_id]);
#endif
      } else cnt_received_rel_not_reciprocal[myrank] = 0;
    }
  }
   
  //! Validate that the correct number of data is received
  void validate_correct_metadata_received(int node_id, uint taxon_length, uint cnt_received_rel_not_reciprocal, uint *list_of_data_to_send, uint *sum_for_all_nodes_length_of_taxa) {
#ifndef NDEBUG
    uint cnt_this_rel = 0, expected_this_size = 0;
    for(uint i = 0; i < taxon_length; i++) {
      if(list_of_data_to_send[i+(node_id *taxon_length)]) {
	cnt_this_rel += list_of_data_to_send[i+(node_id *taxon_length)];
	expected_this_size += sum_for_all_nodes_length_of_taxa[i];
      }
    }
    assert((cnt_this_rel + cnt_received_rel_not_reciprocal) == expected_this_size);    
#endif
  }
  /**
     @brief Transmits the data using non-syncronous communication:
     @remarks Tests for this function are included in the caller of it.
  **/
  void rx_data(int myrank, int number_of_nodes, struct uint_3x_obj *receiver_obj_info,
	       T **&lst_arr_rel_start, index_t **&lst_arr_index_start, struct uint_3x_obj **&lst_meta_taxa_1d,
	       //	       MPI_Datatype &mpi_type_uint_3x, MPI_Datatype &mpi_type_uint_2x, MPI_Datatype &mpi_type_T, 
	       MPI_Request *receive_status, MPI_Request *meta_recv_req, MPI_Request *index_recv_req, MPI_Request *T_recv_req	       
	       ) {
    lst_meta_taxa_1d = new struct uint_3x_obj*[number_of_nodes];
    lst_arr_index_start= new index_t*[number_of_nodes];
    //    uint pos_recv = 0;	  
    lst_arr_rel_start = new T*[number_of_nodes];
    for(int sender_id = 0; sender_id < number_of_nodes; sender_id++) {
      if(sender_id != myrank) { // Does not communicate with "my self"
	//! Does not process until the data is received:
	MPI_Wait(&receive_status[sender_id],MPI_STATUS_IGNORE);

	//! Initiates local variables, making the code similar to the code sending the objects:
	const uint cnt_taxa_pairs = receiver_obj_info[sender_id].taxon_id; // 'taxon_id' represents there the total number of taxon-pairs.
	const uint cnt_total_length_of_index_data = receiver_obj_info[sender_id].total_length_of_index;
	const uint cnt_total_length_of_rel_data   = receiver_obj_info[sender_id].total_length_of_rel;
	if(cnt_total_length_of_rel_data > 0) {
	  //! Allocates, and initiates the list for the lists (buffers), using the "list_of_data_to_send" received:
	  lst_meta_taxa_1d[sender_id] = new struct uint_3x_obj[cnt_taxa_pairs]; 
	  lst_arr_index_start[sender_id] = new index_t[cnt_total_length_of_index_data];
	  lst_arr_rel_start[sender_id] = new T[cnt_total_length_of_rel_data];
	  assert(lst_meta_taxa_1d[sender_id]);
	  assert(lst_arr_index_start[sender_id]);
	  assert(lst_arr_rel_start[sender_id]);
	  //! Receives the length of the subsections the "main list" consists of, to the node(s) redeving data:      
	  MPI_Irecv(lst_meta_taxa_1d[sender_id], cnt_taxa_pairs, mpi_type_uint_3x, sender_id, object_info_in_list, MPI_COMM_WORLD, &meta_recv_req[sender_id]);
	  //! Receives the index_list:
	  MPI_Irecv(lst_arr_index_start[sender_id], cnt_total_length_of_index_data, mpi_type_uint_2x, sender_id, index_list, MPI_COMM_WORLD, &index_recv_req[sender_id]);
	  //! Receives the data of type T
	  MPI_Irecv(lst_arr_rel_start[sender_id], cnt_total_length_of_rel_data, mpi_type_T, sender_id, T_list, MPI_COMM_WORLD, &T_recv_req[sender_id]);
#ifndef NDEBUG
	  //! Receives metadata to be verified by the receiver:
	  assert(debug_sum_distances_total_received);
	  assert(debug_sum_distances_total_received[sender_id] == 0);
	  // As this is the smallest element to send, we assume it has arrived when the longer data has arrived:
	  MPI_Request temp_req;
	  MPI_Irecv(&debug_sum_distances_total_received[sender_id], 1, MPI_FLOAT, sender_id, debug_sum_distance, MPI_COMM_WORLD, &temp_req);	  
#endif

	}
      }
    }
  // TODO: Consider doing something here: While the messages are delivered, we could do computations before the "waitAll" is completed
    // TODO: Consider using a foor-loop with tests, iteratively doing the "merging", until all data are received.

  }



  //! Writes stdout the data received for 'this' node.
  void print_data_received(int myrank, int number_of_nodes, struct uint_3x_obj *receiver_obj_info) {
    for(uint id_receiver = 0; id_receiver < (uint)number_of_nodes; id_receiver++) {
      printf("myrank(%d) received the following from node(%u): pairs(%u), cnt_info(%u) and cnt_rel(%u) (at file %s)\n", myrank, id_receiver, receiver_obj_info[id_receiver].taxon_id, receiver_obj_info[id_receiver].total_length_of_index, receiver_obj_info[id_receiver].total_length_of_rel, __FILE__);
    }	
  }

  //! Receives- and merges the data received into the caller (ie, the list_file_parse object given).
  void rx_list_file_parse_data(struct uint_3x_obj *receiver_obj_info, MPI_Request *receive_status) {
    //! Starts the process of receiving:
    T **lst_arr_rel_start = NULL; 
    struct uint_3x_obj **lst_meta_taxa_1d = NULL;
    index_t **lst_arr_index_start = NULL;
    //! Receives the data using non-syncronous communication:
    MPI_Request meta_recv_req[number_of_nodes], index_recv_req[number_of_nodes],T_recv_req[number_of_nodes];
    rx_data(myrank, number_of_nodes, receiver_obj_info, lst_arr_rel_start, lst_arr_index_start, lst_meta_taxa_1d, 
	    //	    mpi_type_uint_3x, mpi_type_uint_2x, mpi_type_T,
	    receive_status, meta_recv_req, index_recv_req,T_recv_req
	    );
      
    //! Merges the data received into the caller (ie, the list_file_parse object given).
    for(int sender_id = 0; sender_id < number_of_nodes; sender_id++) {
      if(sender_id != myrank) { // Does not communicate with "my self"  
	const uint cnt_total_length_of_rel_data   = receiver_obj_info[sender_id].total_length_of_rel;
	if(cnt_total_length_of_rel_data > 0) {
	  //! Waits before merging is started:
	  MPI_Wait(&meta_recv_req[sender_id], MPI_STATUS_IGNORE);
	  MPI_Wait(&index_recv_req[sender_id], MPI_STATUS_IGNORE);
	  MPI_Wait(&T_recv_req[sender_id], MPI_STATUS_IGNORE);
#ifndef NDEBUG
	  const float sum_before_merging = this_caller->get_sum_of_pair_distances_in_memory();
#endif
	  //! Expands the list
	  merge_data(lst_meta_taxa_1d[sender_id], lst_arr_index_start[sender_id], lst_arr_rel_start[sender_id], receiver_obj_info[sender_id]);

#ifndef NDEBUG
	  const float sum_after_merging = this_caller->get_sum_of_pair_distances_in_memory();
	  const float distances_added = sum_after_merging - sum_before_merging;
	  log_builder::compare_floats(distances_added, debug_sum_distances_total_received[sender_id], cnt_total_length_of_rel_data, __LINE__, __FILE__, __FUNCTION__);
#endif
	  //! Deallocates the objects at the end:
	  if(lst_meta_taxa_1d[sender_id]) {delete [] lst_meta_taxa_1d[sender_id], lst_meta_taxa_1d[sender_id] = NULL;}
	  if(lst_arr_index_start[sender_id]) {delete [] lst_arr_index_start[sender_id], lst_arr_index_start[sender_id] = NULL;}
	  if(lst_arr_rel_start[sender_id]) {delete [] lst_arr_rel_start[sender_id], lst_arr_rel_start[sender_id] = NULL;}
	}
      }
    }
    //! Deallocates the data:
    if(lst_meta_taxa_1d) {delete [] lst_meta_taxa_1d, lst_meta_taxa_1d = NULL;}
    if(lst_arr_index_start) {delete [] lst_arr_index_start, lst_arr_index_start = NULL;}
    if(lst_arr_rel_start) {delete [] lst_arr_rel_start, lst_arr_rel_start = NULL;}
    
  }

 public:
  /**
     @brief Sends only selective data accross nodes.
     @todo Requires that at least the same number of taxa must be found in the file as the number of nodes used for this processing. 
     @remarks Following the procedure:
     -# Each node generates a list-overview of the number of taxon-pairs for each taxa.
     -# The data are exchanged
     -# The node with the highest frequency for a taxon, gets the responsibility; inititates.
     -# (a) The first taxon_length eleemnts holds the id, the second block of taxon_length holds
     -# Performs intra-node-communication, using non-blocking operations (MPI_Isend and MPI_Irecv, with a blocking until all data are received), performing an on-the-flow packing/unpacking of the data (hopefully the time spent doing the packing/unpacking would correspond (??) to the time sending the data).
     -# Transmit- and receive all the data appropriately: Let all the node (xmt 'id') send, and 'id' receive.
  **/
  void mpi_send_only_selective_data() {
#ifndef NDEBUG
    //! Below variable to validate that receiving of data is correct, and thereby conclude that any difference with our expectations is due to mering of data).
    assert(debug_sum_distances_total_received);
    debug_sum_distances_total_received[myrank] = this_caller->get_sum_of_pair_distances_in_memory();

    //! Sums the total value of pairs, in order to verify at the end of 'this operation', that this value is
    //! consistent with the new- and updated value:
    float sum_myrank = this_caller->get_sum_of_pair_distances_in_memory();
    float total_sum_before_merging = 0;
    MPI_Allreduce(&sum_myrank, &total_sum_before_merging, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
    assert(sum_myrank <= total_sum_before_merging);

    //! Sums the total number of pairs:
    uint cnt_myrank_before_merging = this_caller->getTotalLengthOfData();
    uint total_cnt_before_merging  = 0;
    MPI_Allreduce(&cnt_myrank_before_merging, &total_cnt_before_merging, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    assert(cnt_myrank_before_merging <= total_cnt_before_merging);
#endif
    //! Identifies taxa this object is responsible for, ie, the node this object belongs to; set from "mpi_transfer_list_file_parse" object:
    bool *list_of_myrank_taxa_responsilibties = NULL;
    if(!initiate_taxa_responsibilities(list_of_myrank_taxa_responsilibties)) return;

    //! Each node generates a list-overview of the number of taxon-pairs for each taxa,
    //! followed by an an all-reduce operation, summing the total number of data for each taxon.        
    //! Note: 'list_of_data_to_send' holds the complete list of the length for each taxon.
    uint *list_of_data_to_send = NULL, *sum_for_all_nodes_length_of_taxa = NULL;
    tx_rx_complete_lengths_of_taxa(list_of_data_to_send, sum_for_all_nodes_length_of_taxa);
    
    //! Build the lists: The node with the highest frequency for a taxon, gets the responsibility; inititates.
    uint *responsibilities = NULL, *nodes_having_taxon = NULL, *maximum_lengths = NULL;
    set_myranks_taxa_responsibility(responsibilities,nodes_having_taxon, maximum_lengths, list_of_data_to_send, list_of_myrank_taxa_responsilibties);

    // A list of receving- and sending data:
    struct uint_3x_obj sender_obj_info[number_of_nodes];
    struct uint_3x_obj receiver_obj_info[number_of_nodes];
    uint sender_nodes_having_data_for[number_of_nodes]; // Holds the number of nodes that have data to send for a particular node.
    uint **arr_lastPosInOwnData = NULL;
    //! Sends the metadata to the correct node, who receives it:
    uint cnt_received_rel_not_reciprocal[number_of_nodes]; // Value used in debug mode for the receiver, ensuring that the new vount is correct.
    MPI_Request receive_status[number_of_nodes]; // To avoid not processing data before it is received.
    MPI_Request receive_status_debug[number_of_nodes]; // As a safeguard, eveno though it in most cases will not be mandatory.
    tx_rx_taxon_specific_lengths(sender_obj_info, receiver_obj_info, sender_nodes_having_data_for, arr_lastPosInOwnData, list_of_data_to_send, sum_for_all_nodes_length_of_taxa, cnt_received_rel_not_reciprocal, receive_status, receive_status_debug);

    // TODO: The problem with this approach is the assumption that all the data sent will fit in memory, something that might not be the case if its' a huge chunk of "possible inparalogs" sent:
 
    // print_data_received(myrank, number_of_nodes, receiver_obj_info);

    //! Initiates variables for the sending:
    struct uint_3x_obj **nodelist_meta_taxa_1d = new struct uint_3x_obj*[number_of_nodes];
    index_t **nodelist_arr_index_start = new index_t*[number_of_nodes];
    T **nodelist_arr_rel_start = new T*[number_of_nodes];
    for(int i = 0; i < number_of_nodes; i++) {
      nodelist_meta_taxa_1d[i] = NULL, nodelist_arr_index_start[i] = NULL, nodelist_arr_rel_start[i] = NULL;
    }

    //! Variables ensuring safe memory deletion (when sending operation has completed).
    MPI_Request nodelist_meta_send_req[number_of_nodes], nodelist_index_send_req[number_of_nodes], nodelist_T_send_req[number_of_nodes]; 
    //! Every node sends, even if the data is empty, to avoid deadlock:
    for(uint id_receiver = 0; id_receiver < (uint)number_of_nodes; id_receiver++) {
      if(id_receiver != (uint)myrank) { 
	//! De-allocates all the 'memory containers' safely received:
	non_final_deallocate_data_containers_used_for_sending(myrank, number_of_nodes,
						    nodelist_meta_taxa_1d, nodelist_arr_index_start, nodelist_arr_rel_start,
						    nodelist_meta_send_req, nodelist_index_send_req, nodelist_T_send_req);
	//! Transfers the data:
	mpi_transfer_data(sender_obj_info[id_receiver], MPI_COMM_WORLD, myrank, id_receiver, number_of_nodes, 0, taxon_length,
			  this_caller, mpi_Isend_nonblocking, list_of_data_to_send, 		  
			  //! Note 'id_receiver' used below as index instead of 'myrank', as we have stored all the data sent to 'id_receiver' (from 'myrank') at this index.
			  arr_lastPosInOwnData[id_receiver],
			  nodelist_meta_taxa_1d[id_receiver], //mpi_type_uint_3x,
			  nodelist_arr_index_start[id_receiver], //mpi_type_uint_2x,
			  nodelist_arr_rel_start[id_receiver], //mpi_type_T,
			  nodelist_meta_send_req[id_receiver], nodelist_index_send_req[id_receiver], nodelist_T_send_req[id_receiver]
			  );

	// TODO: Validate that it is safe the below 'deletion':
	if(arr_lastPosInOwnData[id_receiver]) {delete [] arr_lastPosInOwnData[id_receiver]; arr_lastPosInOwnData[id_receiver] = NULL;}
      }
    }

    //! Removes content not regarded as being of future interest for this node:
    //! Note: Uses possible lag (wating time), therby this operation might not comsume running time.
    assert(list_of_data_to_send);
    for(uint i = 0; i < taxon_length; i++) {    
      for(uint out = 0; out < taxon_length; out++) {
	if(!(receive_node_is_interested_in_these_data(myrank, i, out, taxon_length, list_of_data_to_send))) {
 	  this_caller->free_memory_for_taxon_pair(i, out);
	}
      }
    }

    //! Receives- and merges the data received into the caller (ie, the list_file_parse object given).
    non_final_deallocate_data_containers_used_for_sending(myrank, number_of_nodes,
							  nodelist_meta_taxa_1d, nodelist_arr_index_start, nodelist_arr_rel_start,
							  nodelist_meta_send_req, nodelist_index_send_req, nodelist_T_send_req);

    rx_list_file_parse_data(receiver_obj_info, receive_status);

    //! Deallocates the remaning data containers used for sending:
    deallocate_data_containers_used_for_sending(myrank, number_of_nodes,
						nodelist_meta_taxa_1d, nodelist_arr_index_start, nodelist_arr_rel_start,
						
						nodelist_meta_send_req, nodelist_index_send_req, nodelist_T_send_req);    

#ifndef NDEBUG    
    //! Waits: As a safeguard, eveno though it in most cases will not be mandatory.    
    for(int i = 0; i < number_of_nodes; i++) {
      if(i != myrank) {MPI_Wait(&receive_status_debug[i],MPI_STATUS_IGNORE);}
    }
    //! Asserts the the correct data is received:
    //! Note: First sums them at this point, due to the asynchrounous way they were received
    uint cnt_received_rel_not_reciprocal_total = 0;
    for(int i = 0; i < number_of_nodes; i++) cnt_received_rel_not_reciprocal_total += cnt_received_rel_not_reciprocal[i];
    validate_correct_metadata_received(myrank, taxon_length, cnt_received_rel_not_reciprocal_total, list_of_data_to_send, sum_for_all_nodes_length_of_taxa);

    //! Builds an overview of the protein pairs expected to be found in the object ('this_caller') the
    //! received data is inserted into.
    //! Note: Executed in two seperate steps (beofre- and after the update of the 'mother object'),
    //! In order to (a) enable the usage and (b) avoid a more complicated opeartion using other type of lists.
    uint expected_this_size = 0;
    for(uint i = 0; i < taxon_length; i++) {
      if(list_of_myrank_taxa_responsilibties[i]) {
	expected_this_size += sum_for_all_nodes_length_of_taxa[i];
      }
    }
#endif

    //!    Update the scheduling algorithm with the list (using macro variables making this a specific option if MPI is compiled into it).
    this_caller->set_list_of_nodes_taxa_responsibilities(list_of_myrank_taxa_responsilibties, taxon_length);
    assert(list_of_myrank_taxa_responsilibties == NULL);

#ifndef NDEBUG
    //! Validate that the 'mother object' consists of all the data in total orignialy found distributed among the nodes:
    const uint pairs_found = this_caller->getTotalLengthOfData_for_myrank();
    if(!(pairs_found == expected_this_size)) {
      fprintf(stderr, "myrank[%d]\tpairs_found(%u), but was expectint_cnt(%u) at line %d in file %s\n", myrank, pairs_found, expected_this_size, __LINE__, __FILE__);
      assert(pairs_found == expected_this_size);
    }

    //! Sums the total number of pairs:
    //! Note: In order to be able comparing number of pairs, we must use the only those for 'myrank'
    uint cnt_myrank_after_merging = this_caller->getTotalLengthOfData_for_myrank();
    uint total_cnt_after_merging  = 0;
    MPI_Allreduce(&cnt_myrank_after_merging, &total_cnt_after_merging, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    assert(cnt_myrank_after_merging <= total_cnt_after_merging);

    //! Sums the total number of distances for the pairs:
    //! Note: In order to be able comparing number of pairs, we must use the only those for 'myrank'
    float sum_myrank_after = this_caller->get_sum_of_pair_distances_for_myrank_in_memory();
    float total_sum_after_merging = 0;
    MPI_Allreduce(&sum_myrank_after, &total_sum_after_merging, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );

    if(total_cnt_before_merging != total_cnt_after_merging) {
      fprintf(stderr, "!!\t [myrank=%d]\ttotal_cnt_before_merging(%u) != total_cnt_after_merging(%u), at %s:%d\n", myrank, total_cnt_before_merging, total_cnt_after_merging, __FILE__, __LINE__);
    }
    assert(sum_myrank_after <= total_sum_after_merging);
    log_builder::compare_floats(total_sum_before_merging, total_sum_after_merging, total_cnt_after_merging, __LINE__, __FILE__, __FUNCTION__);
#endif
    
    //! Deallocates the pointer used:
    //! Note: The below containers are first de-allocated here due to (a) their relatively low memory consumption and (b) to enable easy validation (after their intented usage):
    if(list_of_data_to_send) {delete [] list_of_data_to_send; list_of_data_to_send = NULL;}
    if(responsibilities) {delete [] responsibilities; responsibilities = NULL; nodes_having_taxon = NULL; maximum_lengths = NULL;}
    if(arr_lastPosInOwnData) {delete [] arr_lastPosInOwnData; arr_lastPosInOwnData = NULL;}
    if(nodelist_meta_taxa_1d) {delete [] nodelist_meta_taxa_1d, nodelist_meta_taxa_1d = NULL;}
    if(nodelist_arr_index_start) {delete [] nodelist_arr_index_start, nodelist_arr_index_start = NULL;}
    if(nodelist_arr_rel_start) {delete [] nodelist_arr_rel_start, nodelist_arr_rel_start = NULL;}	  
#ifndef NDEBUG
    if(sum_for_all_nodes_length_of_taxa) {delete [] sum_for_all_nodes_length_of_taxa; sum_for_all_nodes_length_of_taxa = NULL;}
#endif
  }

 private:
  //! Sends- and receives meatadata in wa manner similar to bcast:
  void tx_rx_equal_metadata_to_all(struct uint_3x_obj *rx_and_tx_obj_info, uint *sender_nodes_having_data_for, uint *&arr_lastPosInOwnData, uint *cnt_received_rel_not_reciprocal, MPI_Request *receive_status, MPI_Request *receive_status_debug) {
    // Does a processing before, building set of metadata equal for all the nodes to send data to:
    arr_lastPosInOwnData = new uint[number_of_nodes];
    log_builder::test_memory_condition_and_if_not_abort(arr_lastPosInOwnData!=NULL, __LINE__, __FILE__, __FUNCTION__);
    memset(arr_lastPosInOwnData, 0, sizeof(uint)*number_of_nodes);

    //! Get metadata information about 'this' node:
    uint this_cnt_value_not_reciprocal = 0;
    arr_lastPosInOwnData = get_metadata_information(rx_and_tx_obj_info[myrank], this_cnt_value_not_reciprocal);
    for(uint node_id = 0; node_id < (uint)number_of_nodes; node_id++) {	
      sender_nodes_having_data_for[node_id] = 0;
      if(node_id != (uint)myrank) {
	rx_and_tx_obj_info[node_id] = uint_3x_obj_t();
	//! Sends the length of the lists to the node receiving data:	
	MPI_Request send_req; MPI_Isend(&rx_and_tx_obj_info[myrank], 1, mpi_type_uint_3x, node_id, object_info, MPI_COMM_WORLD, &send_req);
#ifndef NDEBUG
	MPI_Request send_req_debug; MPI_Isend(&this_cnt_value_not_reciprocal,1,MPI_UNSIGNED, node_id, 1, MPI_COMM_WORLD, &send_req_debug);
#endif

	//! Receives from all the nodes:
	MPI_Irecv(&rx_and_tx_obj_info[node_id], 1, mpi_type_uint_3x, node_id, object_info, MPI_COMM_WORLD, &receive_status[node_id]);
#ifndef NDEBUG
	MPI_Irecv(&cnt_received_rel_not_reciprocal[node_id], 1, MPI_UNSIGNED, node_id, 1, MPI_COMM_WORLD, &receive_status_debug[node_id]);
#endif
      }
    }
  }

 public:
  /**
     @brief Sends all of 'this' data to al the other nodes, and receives the data from all the other nodes.
  **/
  void tx_rx_all_data_to_all_nodes() {
#ifndef NDEBUG
    //! Sums the total value of pairs, in order to verify at the end of 'this operation', that this value is
    //! consistent with the new- and updated value:    
    float sum_myrank = this_caller->get_sum_of_pair_distances_in_memory();
    float total_sum_before_merging = 0;
    MPI_Allreduce(&sum_myrank, &total_sum_before_merging, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
    assert(sum_myrank <= total_sum_before_merging);
    //    tx_rx_complete_lengths_of_taxa(list_of_data_to_send, sum_for_all_nodes_length_of_taxa);
    
    uint cnt_T_before_myrank = this_caller->getTotalLengthOfData();
    uint total_sum_T_before_merging = cnt_T_before_myrank;
    MPI_Allreduce(&cnt_T_before_myrank, &total_sum_T_before_merging, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    assert(cnt_T_before_myrank <= total_sum_T_before_merging);
#endif

    // A list of receving- and sending data:
    //    struct uint_3x_obj sender_obj_info; // = uint_3x_obj(); //[number_of_nodes];
    //! Below container used for both sending- and receiving
    struct uint_3x_obj rx_and_tx_obj_info[number_of_nodes];
    uint sender_nodes_having_data_for[number_of_nodes]; // Holds the number of nodes that have data to send for a particular node.
    //! Note: Stores all the data sent in the below variable
    uint *arr_lastPosInOwnData = NULL;
    //! Sends the metadata to the correct node, who receives it:
    uint cnt_received_rel_not_reciprocal[number_of_nodes]; // Value used in debug mode for the receiver, ensuring that the new vount is correct.
    MPI_Request receive_status[number_of_nodes]; // To avoid not processing data before it is received.
    MPI_Request receive_status_debug[number_of_nodes]; // As a safeguard, eveno though it in most cases will not be mandatory.
    tx_rx_equal_metadata_to_all(rx_and_tx_obj_info, sender_nodes_having_data_for, arr_lastPosInOwnData, cnt_received_rel_not_reciprocal, receive_status, receive_status_debug);


    //
    //! Initiates variables for the sending:
    struct uint_3x_obj *meta_taxa_1d = NULL;
    index_t *arr_index_start = NULL; 
    T *arr_rel_start = NULL; //new T*[number_of_nodes];

    //! Transforms the 2d-objects into 1d:
    uint cnt_index_inserted = 0, cnt_rel_inserted = 0, elements_used = 0;
    //    uint cnt_taxa_pairs = 0;
    assert(!meta_taxa_1d);
    assert(!arr_index_start);
    assert(!arr_rel_start);
    uint cnt_taxa_pairs = rx_and_tx_obj_info[myrank].taxon_id, cnt_total_length_of_rel_data = rx_and_tx_obj_info[myrank].total_length_of_rel, cnt_total_length_of_index_data = rx_and_tx_obj_info[myrank].total_length_of_index;

    transform_2d_into_1d(this_caller,
			 listTaxa, this_caller->get_buffer(), 0, taxon_length, cnt_taxa_pairs, mpi_Isend_nonblocking,
			 taxon_length, arr_lastPosInOwnData, elements_used, meta_taxa_1d,
			 UINT_MAX, NULL,
			 // Data for the T type:
			 cnt_rel_inserted, arr_rel_start, cnt_total_length_of_rel_data,
			 // Data for the index_t type:
			 cnt_index_inserted, arr_index_start, cnt_total_length_of_index_data);
    assert(elements_used == cnt_taxa_pairs); // The number of inserted elements should correspond to the reserved length.

    //! Variables ensuring safe memory deletion (when sending operation has completed).
    MPI_Request meta_send_req[number_of_nodes], index_send_req[number_of_nodes], T_send_req[number_of_nodes]; 
    //! Every node sends, even if the data is empty, to avoid deadlock:
    for(uint id_receive_node = 0; id_receive_node < (uint)number_of_nodes; id_receive_node++) {
      if(id_receive_node != (uint)myrank) { 
	//! Transfers the data:
	MPI_Isend(meta_taxa_1d, cnt_taxa_pairs, mpi_type_uint_3x, id_receive_node, object_info_in_list,
		  MPI_COMM_WORLD, &meta_send_req[id_receive_node]);
	MPI_Isend(arr_index_start, cnt_total_length_of_index_data, mpi_type_uint_2x, id_receive_node, index_list,
		  MPI_COMM_WORLD, &index_send_req[id_receive_node]);
	MPI_Isend(arr_rel_start, cnt_total_length_of_rel_data, mpi_type_T, id_receive_node, T_list,
		  MPI_COMM_WORLD, &T_send_req[id_receive_node]);
      }
    }

    //! Receives- and merges the data received into the caller (ie, the list_file_parse object given).
    rx_list_file_parse_data(rx_and_tx_obj_info, receive_status);

    //! Waits: As a safeguard, eveno though it in most cases will not be mandatory.    
    for(int i = 0; i < number_of_nodes; i++) {
      if(i != myrank) {
	MPI_Wait(&meta_send_req[i],  MPI_STATUS_IGNORE);
	MPI_Wait(&index_send_req[i], MPI_STATUS_IGNORE);
	MPI_Wait(&T_send_req[i],     MPI_STATUS_IGNORE);
      }
    }
    if(meta_taxa_1d) {delete [] meta_taxa_1d, meta_taxa_1d = NULL;}
    if(arr_index_start) {delete [] arr_index_start, arr_index_start = NULL;}
    if(arr_rel_start) {delete [] arr_rel_start, arr_rel_start = NULL;}	  

#ifndef NDEBUG    
    //! Waits: As a safeguard, eveno though it in most cases will not be mandatory.    
    for(int i = 0; i < number_of_nodes; i++) {
      if(i != myrank) {MPI_Wait(&receive_status[i],MPI_STATUS_IGNORE);}
    }
    
    //! Validates the contents of the metadata containers:
    uint total_sum_T_after_merging = 0;
    for(uint node_id = 0; node_id < (uint)number_of_nodes; node_id++) {
      total_sum_T_after_merging += rx_and_tx_obj_info[node_id].total_length_of_rel;
    }
    if(total_sum_T_before_merging != total_sum_T_after_merging) {
      fprintf(stderr, "!!\t[myrank=%d]\t total_sum_T_before_merging(%u) and total_sum_T_after_merging(%u) at line %d in file %s\n", myrank, total_sum_T_before_merging, total_sum_T_after_merging, __LINE__, __FILE__);
      assert(total_sum_T_before_merging == total_sum_T_after_merging);
    }

    uint cnt_T_after_merge_myrank = this_caller->getTotalLengthOfData();
    uint total_sum_T_after_merge_merging = cnt_T_after_merge_myrank;
    MPI_Allreduce(&cnt_T_after_merge_myrank, &total_sum_T_after_merge_merging, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    assert(cnt_T_after_merge_myrank <= total_sum_T_after_merge_merging);
    
    const uint val_mult = (number_of_nodes*total_sum_T_after_merging);
    //    printf( "--\t[myrank=%d]\t total_sum_T_after_merge_merging(%u) and %u=(nodes(%u)*total_sum_T_after_merging(%u)) at line %d in file %s\n", myrank, total_sum_T_after_merge_merging, val_mult, (uint)number_of_nodes, total_sum_T_after_merging, __LINE__, __FILE__);
    if(total_sum_T_after_merge_merging != val_mult) {
      fprintf(stderr, "!!\t[myrank=%d]\t total_sum_T_after_merge_merging(%u) and %u=(nodes(%u)*total_sum_T_after_merging(%u)) at line %d in file %s\n", myrank, total_sum_T_after_merge_merging, val_mult, (uint)number_of_nodes, total_sum_T_after_merging, __LINE__, __FILE__);
      assert(total_sum_T_after_merge_merging == val_mult);
    }
#endif

    //! De-allocates the memory reserved
    if(arr_lastPosInOwnData) {delete [] arr_lastPosInOwnData; arr_lastPosInOwnData = NULL;}
  }

  //! Deallocates memory
  void free_memory() {
    MPI_Type_free(&mpi_type_uint_3x);
    MPI_Type_free(&mpi_type_uint_2x);    
    MPI_Type_free(&mpi_type_T); 
#ifndef NDEBUG
    if(debug_sum_distances_total_received) {delete [] debug_sum_distances_total_received; debug_sum_distances_total_received = NULL;}
    if(debug_sum_distances_total_to_send) {delete [] debug_sum_distances_total_to_send; debug_sum_distances_total_to_send = NULL;}
#endif
  }

  //! The constructor:
  mpi_transfer_list_file_parse(list_file_parse<T> *_this_caller) :
    myrank(0), number_of_nodes(0), this_caller(_this_caller), taxon_length(0) 
#ifndef NDEBUG
    , debug_sum_distances_total_to_send(NULL), debug_sum_distances_total_received(NULL)
#endif
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
      MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes);
      assert(this_caller);
      //      taxon_length = (uint)this_caller->get_taxon_length();
      //! Builds mpi datatype to ensure consistencies accross different hardware platforms:
      mpi_build_datatype_3_uint(mpi_type_uint_3x);
      mpi_build_datatype_2_uint(mpi_type_uint_2x);    
      T *temp = new T; build_datatype_T(temp, mpi_type_T); if(temp) {delete temp; temp = NULL;}    
      int temp_taxon_length = 0; listTaxa = this_caller->get_listTaxa(temp_taxon_length);
      taxon_length = (uint)temp_taxon_length;
#ifndef NDEBUG
      assert(number_of_nodes);
      debug_sum_distances_total_to_send = new float[number_of_nodes];
      assert(debug_sum_distances_total_to_send);
      memset(debug_sum_distances_total_to_send, 0, sizeof(uint)*number_of_nodes);	
      debug_sum_distances_total_received = new float[number_of_nodes];
      assert(debug_sum_distances_total_received);
      memset(debug_sum_distances_total_received, 0, sizeof(uint)*number_of_nodes);	
#endif
    }; 
    
};

#endif
#endif
