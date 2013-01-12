#ifndef mpi_transfer_list_norm_h
#define mpi_transfer_list_norm_h
#ifdef USE_MPI
/**
   @file
   @brief Holds the enitre list_norm object
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
#include "mpi_transfer_norm.h"
/**
   @struct mpi_transfer_list_norm
   @brief Holds the enitre list_norm object
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
struct mpi_transfer_list_norm {    
  // -------------------------------------------------
  // The below values are to be sent in the first 'run': 
  //
  //! max_input_value to be found the max of.
  float max_input_value;
  //! taxon_length to be used converting the 1d-indexes into 2d-indexes.
  uint taxon_length;    
  //! The size of the list to transfer.
  uint size_list;
  // -------------------------------------------------

  //! The list of norm objects to transfer
  mpi_transfer_norm_t *list;
  
  private:
  enum send_types {object_info, object_lists};
  MPI_Datatype mpi_type_info;
  MPI_Datatype mpi_type_norm_array;
  public:
  // Prints info about the settings for this object:
  void print_info(int myrank) {
    printf("[myrank=%d]\t max_input_value(%f), taxon_length(%u), size_list(%u) at line %d in file %s\n", myrank, max_input_value, taxon_length, size_list, __LINE__, __FILE__);
  }
  
  //! Gets the 2d-coordinate given the id.
  void get_2d_coordinate(uint id, uint &in, uint &out, uint taxon_length) {
    if(id) {
      assert(taxon_length);
      out = id % taxon_length;
      //            in = 0; printf("--\t id(%u) --> in(%u) and out(%u) at line %d in file %s\n", id, in, out, __LINE__, __FILE__);
      assert(id >= out);
      in = (id - out) / taxon_length;
    } else {in = 0, out = 0;}
    //    printf("--\t id(%u) --> in(%u) and out(%u), for taxon_length(%u) at line %d in file %s\n", id, in, out, taxon_length, __LINE__, __FILE__);
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


  private:
  void initiate_list_to_mimimum_size(uint min_size) {
    assert(taxon_length);
    assert(!list);
    if(!list) {
      size_list = min_size;
      assert(size_list);
      list = mpi_transfer_norm_t::init(size_list);
    }// else if (size_list < min_size
    assert(list);
  }

  public:
  /**
     @brief Builds the list:
     @remarks
     # The list could be huge, therefore tries making it as compact as possible
  **/
  void build_list(float _max_input, uint _t_length, norm_t ***norm_list) {
    max_input_value = _max_input; taxon_length = _t_length;       
    assert(taxon_length);
    assert(norm_list);
    size_list = taxon_length*taxon_length;
    if(taxon_length > 1000) size_list = (uint)(size_list*0.1); // Estimates an filling degree of 10%
    assert(taxon_length);
    if(!list) {
      assert(size_list);
      list = mpi_transfer_norm_t::init(size_list);
    }
    uint elements_used = 0;
    for(uint i = 0; i < taxon_length; i++) {
      if(norm_list[i]) {
	for(uint out = 0; out < taxon_length; out++) {
	  if(elements_used == size_list) {
	    mpi_transfer_norm_t::enlarge(list, size_list);
	    assert(size_list > elements_used);
	  }
	  if(norm_list[i][out]) {
	    uint id = get_1d_coordinate(i, out, taxon_length);
	    list[elements_used] = mpi_transfer_norm(id, norm_list[i][out]->get_cnt(), norm_list[i][out]->get_sum(), norm_list[i][out]->get_cnt_zeros());
	    elements_used++;		
	  }
	}
      }
    }
    size_list = elements_used; // "Resizes"
#ifndef NDEBUG
    //! Validates that the new structure correspons to the old one:
    if(norm_list) {
      norm_t ***recip = get_norm_2d_list(list, taxon_length);
      assert(recip);
      for(uint i = 0; i < taxon_length; i++) {
	if(norm_list[i]) {
	  for(uint out = 0; out < taxon_length; out++) {
	    if(norm_list[i][out]) {
	      norm_t *element = recip[i][out]; ////this_obj->get_list_element(i, out);
	      assert(element);
	      assert(norm_list[i][out]->is_equal(element));
	    }
	  }
	}
      }
      //! Deallocates the temporary norm object:
      uint temp_length = taxon_length;
      norm::close(recip, temp_length);
    }  	   
#endif
  }
  private:
  //! @return a transformation of the internal 1d-list into an matrix of norm objects.  
  norm_t ***get_norm_2d_list(mpi_transfer_norm *list, uint taxon_length) {
    if(taxon_length && list) {
      norm ***recip_list = norm::init(taxon_length);
      for(uint elements_used = 0; elements_used < size_list; elements_used++) {
	const uint id = list[elements_used].id;
	uint in, out; get_2d_coordinate(id, in, out, taxon_length);
	recip_list[in][out] = new norm(list[elements_used].cnt, list[elements_used].sum, list[elements_used].cnt_zero);
      }
      return recip_list;
    } else return NULL;
  }
  //! Does a transformation of the internal 1d-list into an matrix of norm objects.  
  void update_norm_2d_list(norm_t ***recip_list, mpi_transfer_norm *list, uint taxon_length) {
    assert(recip_list);
    if(taxon_length && list) {
      for(uint elements_used = 0; elements_used < size_list; elements_used++) {
	const uint id = list[elements_used].id;
	uint in, out; get_2d_coordinate(id, in, out, taxon_length);
	assert(recip_list[in]);
	if(recip_list[in][out]) recip_list[in][out]->merge(list[elements_used].cnt, list[elements_used].sum, list[elements_used].cnt_zero);
	else recip_list[in][out] = new norm(list[elements_used].cnt, list[elements_used].sum, list[elements_used].cnt_zero);
      }
    }
  }
  //! Does a transformation of the internal 1d-list into an matrix of norm objects.  
  void update_2d_list_ignore_old_values(norm_t ***recip_list, mpi_transfer_norm *list, uint taxon_length) {
    assert(recip_list);
    if(taxon_length && list) {
      for(uint elements_used = 0; elements_used < size_list; elements_used++) {
	const uint id = list[elements_used].id;
	uint in, out; get_2d_coordinate(id, in, out, taxon_length);
	assert(recip_list[in]);
	if(recip_list[in][out]) *recip_list[in][out] = norm(list[elements_used].cnt, list[elements_used].sum, list[elements_used].cnt_zero);
	else recip_list[in][out] = new norm(list[elements_used].cnt, list[elements_used].sum, list[elements_used].cnt_zero);
      }
    }
  }
  public:

  //! Sends the content of this object to the root node.
  void send_data_to_root(MPI_Comm comm, int myrank, int number_of_nodes) {
    assert(myrank != 0);
    const int receive_id = 0;
    assert(myrank != 0);
    enum send_types tag = object_info;
    int ret_val = MPI_Send(this, 1, mpi_type_info, receive_id, tag, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__);
    assert(list);
    tag = object_lists;
    ret_val = MPI_Send(list, size_list, mpi_type_norm_array, receive_id, tag, comm);
    if(ret_val) printf("[myrank=%d]\tgot ret_value(%d) after doing a sending at line %d in file %s\n", myrank, ret_val, __LINE__, __FILE__);
  }

  //! Gets the content of this object, and merges it into the exsisting list_norm object
  void get_and_merge_data(list_norm *this_obj, MPI_Comm comm, int myrank, int number_of_nodes) {
    assert(myrank == 0);
    assert(taxon_length);
    norm_t ***norm_2d = norm::init(taxon_length);
    //! Receives the length of the lists:
    for(int sender_id = 1; sender_id < (int)number_of_nodes; sender_id++) {
      //! Gets the length of the lists:
      mpi_transfer_list_norm t_settings = mpi_transfer_list_norm(taxon_length);
      MPI_Status status_object_info;
      enum send_types tag = object_info;
      MPI_Recv(&t_settings, 1, mpi_type_info, sender_id, tag, comm, &status_object_info) ;
      //! Update the max-input-value:
      //      if(t_settings.max_input_value > max_input_value) max_input_value = t_settings.max_input_value;

      //! Prepares for receiving the list:
      if(t_settings.size_list > 0) {
	initiate_list_to_mimimum_size(t_settings.size_list); size_list = t_settings.size_list;
	assert(list);
      } else size_list = 0;
      //! Receives the lists, and merges them:
      tag = object_lists;

      MPI_Recv(list, size_list, mpi_type_norm_array, sender_id, tag, comm, &status_object_info) ;
      //! Includes the data got into the list:
      update_norm_2d_list(norm_2d, list, t_settings.taxon_length);
      // Prepares for the next run by deallocating temporary containers:
      if(list) {delete [] list; list = NULL;}
      t_settings.free_mem();      
    }
    //! Merges the list with the exsisting:
    if(norm_2d) {
      list_norm temp = list_norm();
      temp.set_data(norm_2d, taxon_length, max_input_value);
      this_obj->merge_basis(temp);
      norm::close(norm_2d, taxon_length);
    }
  }

  //! Gets the content of this object, and merges it into the exsisting list_norm object
  void broadcast_get_and_merge_data(list_norm *this_obj, MPI_Comm comm, int myrank, int number_of_nodes) {
    assert(myrank != 0);
    assert(myrank != 0);
    const int root_id = 0;
    mpi_transfer_list_norm t_settings = mpi_transfer_list_norm(taxon_length);
    MPI_Bcast(&t_settings, 1, mpi_type_info, root_id, comm);
    //! Prepares for receiving the list:
    if(list){    delete [] list; list = NULL; size_list = 0;}
    initiate_list_to_mimimum_size(t_settings.size_list); size_list = t_settings.size_list;
    MPI_Bcast(list, size_list, mpi_type_norm_array, root_id, comm) ;
    assert(list);

    //! Gets the list to work on:
    norm_t ***norm_2d; this_obj->get_norm_2d(norm_2d);
    assert(norm_2d);
    //! Merges the list with the exsisting:
    update_2d_list_ignore_old_values(norm_2d, list, t_settings.taxon_length);
      // Prepares for the next run by deallocating temporary containers:
    if(list) {delete [] list; list = NULL;}
    t_settings.free_mem();      
  }
  void assert_result(list_norm *this_obj, MPI_Comm comm, int myrank, int number_of_nodes) {
#ifndef NDEBUG
    loint my_cnt_ = 0; uint my_cnt=0; float my_dist_sum = 0; uint my_zeros = 0;
    uint global_cnt=0; float global_dist_sum = 0; uint global_zeros = 0;
    this_obj->get_object_sum(my_cnt_, my_dist_sum, my_zeros); my_cnt = (uint)my_cnt;
    const uint my_cnt_mult       = my_cnt * number_of_nodes;
    const float my_dist_sum_mult = my_dist_sum * number_of_nodes;
    const uint my_zeros_mult     = my_zeros * number_of_nodes;
    MPI_Allreduce(&my_cnt, &global_cnt, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&my_dist_sum, &global_dist_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&my_zeros, &global_zeros, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    if(my_cnt_mult != global_cnt) {
      fprintf(stderr, "!!\t[myrank=%d]\t my_cnt_mult(%u) != global_cnt(%u) at line %d in file %s\n", myrank, my_cnt_mult, global_cnt, __LINE__, __FILE__);
      assert(my_cnt_mult == global_cnt);
    }
    log_builder::compare_floats(my_dist_sum_mult, global_dist_sum, global_cnt, __LINE__, __FILE__, __FUNCTION__);
/*     if((int)my_dist_sum_mult != (int)global_dist_sum) { */
/*       fprintf(stderr, "!!\t[myrank=%d]\t my_dist_sum_mult(%f) != global_dist_sum(%f) at line %d in file %s\n", myrank, my_dist_sum_mult, global_dist_sum, __LINE__, __FILE__); */
/*       assert(my_dist_sum_mult == global_dist_sum); */
/*     } */
    if(my_zeros_mult != global_zeros) {
      fprintf(stderr, "!!\t[myrank=%d]\t my_zeros_mult(%u) != global_zeros(%u) at line %d in file %s\n", myrank, my_zeros_mult, global_zeros, __LINE__, __FILE__);
      assert(my_zeros_mult == global_zeros);
    }
#endif
  }

  //! Broadcasts the merge data to all the nodes.
  void broadcast_send_merged_data(MPI_Comm comm, int myrank, int number_of_nodes) {
    assert(myrank == 0);
    const int root_id = 0;
    MPI_Bcast(this, 1, mpi_type_info, root_id, comm);
    assert(list);
    MPI_Bcast(list, size_list, mpi_type_norm_array, root_id, comm) ;
  }

  //! Deallocates the internal list.
  void deallocate_internal_list(){
    if(list) mpi_transfer_norm::close(list, size_list);
  }
  void free_mem() {
    if(list) mpi_transfer_norm::close(list, size_list);
    MPI_Type_free(&mpi_type_info);
    MPI_Type_free(&mpi_type_norm_array);    
  }
  private:
  void build_Datatype_info() {
    int block_lengths[3] = {1,1,1};
    MPI_Datatype old_types[3] = {MPI_FLOAT, MPI_UNSIGNED, MPI_UNSIGNED};
    //! The location of each element;
    MPI_Aint indices[3] = {0, sizeof(float), sizeof(uint)};
    indices[2] += indices[1]; 
    //! Builds the mpi representation of the structure:
    MPI_Type_struct(3, block_lengths, indices, old_types, &mpi_type_info);
    MPI_Type_commit(&mpi_type_info);
  }

  void build_Datatype_norm_array() {
    int block_lengths[4] = {1,1,1, 1};
    MPI_Datatype old_types[4] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_FLOAT, MPI_UNSIGNED};
    //! The location of each element;
    // TODO: Is the below correcT: Should it not be "   {0, sizeof(uint), sizeof(uint), sizeof(float)} " ??;
    MPI_Aint indices[4] = {0, sizeof(uint), sizeof(uint), sizeof(float)};
    //    MPI_Aint indices[4] = {0, sizeof(uint), sizeof(float), sizeof(uint)};
    for(uint i = 1; i < 4;i++) {indices[i] +=indices[i-1];}
    //! Builds the mpi representation of the structure:
    MPI_Type_struct(4, block_lengths, indices, old_types, &mpi_type_norm_array);
    MPI_Type_commit(&mpi_type_norm_array);
  }
  public:
  //! The constructor:
  mpi_transfer_list_norm(uint _taxon_length) : max_input_value(0.0), taxon_length(_taxon_length), size_list(0), list(NULL)
  {
    build_Datatype_info();
    build_Datatype_norm_array();
  }
    
  };

/**
   @brief Holds the enitre list_norm object
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
typedef struct mpi_transfer_list_norm mpi_transfer_list_norm_t;
#endif
#endif
