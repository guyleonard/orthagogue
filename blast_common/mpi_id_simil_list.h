#ifndef mpi_id_simil_list_h
#define mpi_id_simil_list_h
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
#include "id_simil_list.h"
#include "mpi_transfer_rel.h"
/**
   @class mpi_id_simil_list
   @brief Contains methods performing send-, receive- and merge operationd for the class id_simil_list.
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
class mpi_id_simil_list {    
 private:
  int myrank;
  int number_of_nodes;
  int *total_receive_size_for_the_nodes;
  uint internal_list_1d_size;
  rel_t *internal_list_1d;
  //! 'map_internal_list' is the list holding the specific sizes for each protein:
  uint *map_internal_list;
  rel_t *global_1d_list;
  uint global_1d_list_size;
  uint *global_map_list;
  taxa *listTaxa;
  int taxon_length;
  //! @return the number of elements for the given object, using class id_simil_list internal syntax
  uint get_size(rel_t *obj) {
    const uint sum = id_simil_list::get_number_of_pairs(obj);
    //! Note: '+1' to add place for the metadata object
    if(sum) {return sum+1;}
    else return 0;
  }
  //! @return the number of elements for the given list, using class id_simil_list internal syntax
  uint get_size(rel_t **arr, uint arr_size) {    
    uint total_size = 0;
    if(arr && arr_size) {
      for(uint i = 0; i < arr_size; i++) {
	total_size += get_size(arr[i]);
      }
    }
    return total_size;    
  }

  //! @return the number of total number of pairs for the given list, using class id_simil_list internal syntax
  uint get_total_number_of_pairs(rel_t **arr, uint arr_size) {    
    uint total_size = 0;
    if(arr && arr_size) {
      for(uint i = 0; i < arr_size; i++) {
	total_size += id_simil_list::get_number_of_pairs(arr[i]);
      }
    }
    return total_size;    
  }


  /**
     @return the total number of pairs, given the input
  **/
  float get_total_number_of_distances(rel *node_1d_list, uint node_1d_list_size, uint *node_1d_map, uint arr_size) {
    uint cnt_pairs=0; float cnt_sum = 0;
    find_total_number_of_values(node_1d_list, node_1d_list_size, node_1d_map, arr_size, cnt_pairs, cnt_sum);
    return cnt_sum;
  }

  /**
     @return the total number of pairs, given the input
  **/
  uint get_total_number_of_pairs(rel *node_1d_list, uint node_1d_list_size, uint *node_1d_map, uint arr_size) {
    uint cnt_pairs=0; float cnt_sum = 0;
    find_total_number_of_values(node_1d_list, node_1d_list_size, node_1d_map, arr_size, cnt_pairs, cnt_sum);
    return cnt_pairs;
  }
  /**
     @brief Sets the total number of pairs- and distances, given the input
  **/
  void find_total_number_of_values(rel *node_1d_list, uint node_1d_list_size, uint *node_1d_map, uint arr_size, uint &cnt_elements_added, float &sum_distances) {
    assert(!cnt_elements_added);
    assert(!sum_distances);
    if(node_1d_list && node_1d_list_size && node_1d_map && arr_size) {
      rel *node_1d_list_end = node_1d_list + node_1d_list_size;
      uint total_sum_of_lengths_added = 0;
      for(uint i = 0; i < arr_size; i++) {
	if(node_1d_map[i]) {
	  rel *current = node_1d_list + total_sum_of_lengths_added;
	  assert(current < node_1d_list_end);
	  if(current) {
	    const uint size = get_size(current); //id_simil_list::get_number_of_pairs(temp[i]);
	    if(size) {
	      cnt_elements_added += size-1;
	      for(uint out = 1; out < size; out++) {
		sum_distances += current[out].distance;
	      }
	    }
	  }
	  total_sum_of_lengths_added += node_1d_map[i];
	}
      }
    }
  }

  //! @return the number of total sum of pairs, given the global list:
  float get_total_number_of_global_distances(uint arr_size) {
    uint cnt_pairs=0; float cnt_sum = 0;
    find_total_number_of_global_values(arr_size, cnt_pairs, cnt_sum);
    return cnt_sum;
  }

  //! @return the number of total number of pairs, given the global list:
  uint get_total_number_of_global_pairs(uint arr_size) {
    uint cnt_pairs=0; float cnt_sum = 0;
    find_total_number_of_global_values(arr_size, cnt_pairs, cnt_sum);
    return cnt_pairs;
  }
  //! Sets the number of total number of pairs, and their sum, given the global list:
  void find_total_number_of_global_values(uint arr_size, uint &cnt_relations_inserted, float &cnt_distances) {
    assert(!cnt_relations_inserted);
    assert(!cnt_distances);
    for(int node_id = 0; node_id < number_of_nodes; node_id++) {
     rel *node_1d_list = NULL; uint node_1d_list_size = 0; uint *node_1d_map = NULL;
     uint elements_inserted_this_node = 0; float cnt_distances_this = 0;
      if(get_node_1d_list_settings(node_id, node_1d_list, node_1d_list_size, node_1d_map, arr_size)) {
	find_total_number_of_values(node_1d_list, node_1d_list_size, node_1d_map, arr_size, elements_inserted_this_node, cnt_distances_this);
	cnt_relations_inserted += elements_inserted_this_node;
	cnt_distances += cnt_distances_this;
	
      }
#ifndef NDEBUG
      rel **temp = convert_global_1d_into_2d_rel(node_id, arr_size);
      uint temp_cnt_relations_inserted = 0;
      //! Note: Not all nodes have data; arises when #(nodes) > #(taxa).
      if(temp) { 
	for(uint i = 0; i < arr_size; i++) {
	  if(temp[i]) {
	    uint size = get_size(temp[i]); //id_simil_list::get_number_of_pairs(temp[i]);
	    if(size)    temp_cnt_relations_inserted += size-1;
	  }
	}
	de_allocate(temp, arr_size);
      }
      if(!(temp_cnt_relations_inserted == elements_inserted_this_node)) {
	fprintf(stderr, "!!\t[myrank=%d]\t temp_cnt_relations_inserted(%u) != elements_inserted_this_node(%u), at line %d in file %s\n", myrank, temp_cnt_relations_inserted, elements_inserted_this_node, __LINE__, __FILE__);
	assert(temp_cnt_relations_inserted == elements_inserted_this_node);
      }
#endif
    }
  }
  //! Sets the given size.
  void set_size(rel_t *&arr, const uint new_size) {
    assert(arr);
    assert(arr[0].ind_out < new_size);
    arr[0].distance = new_size;
  }
  //! Makes the list tight by remvoing information about empty space:
  void compress_list(rel_t **&arr, uint &arr_size) {
    if(arr && arr_size) {
      for(uint i = 0; i < arr_size; i++) {
	const uint new_size = get_size(arr[i]);
	if(new_size) {
	  set_size(arr[i], new_size);
	}
      }
    }
  }

  /**
     Converts the global list, but does not deallocate data:
     @remarks Include tests, comparing the new- with the old.
  **/
  rel **convert_global_1d_into_2d_rel(int node_id, int *displs, uint arr_size) {
    return convert_1d_into_2d_rel(global_1d_list+displs[node_id], total_receive_size_for_the_nodes[node_id], global_map_list+(node_id*arr_size), arr_size);
  }

  /**
     @return True if mapping was possible,
     @remarks Performs the mapping (ie, if the data is not set, false is returned).
  **/
  bool get_node_1d_list_settings(int node_id, rel *&node_1d_list, uint &node_1d_list_size, uint *&node_1d_map, uint arr_size) {    
    if(arr_size && global_1d_list_size) {
      assert(!node_1d_list);
      assert(!node_1d_list_size);
      assert(!node_1d_map);
      assert(arr_size);
      assert(global_1d_list);    
      assert(global_1d_list_size);
      assert(global_map_list);
      loint offset = 0; for(int i=0; i<node_id; i++) {offset+=total_receive_size_for_the_nodes[i];}
      node_1d_list = global_1d_list + offset;
      node_1d_list_size = (uint)total_receive_size_for_the_nodes[node_id];
      node_1d_map = global_map_list + (node_id*arr_size);
      return true;
    } else return false;
  }
  /**
     Converts the global list, but does not deallocate data:
     @remarks Include tests, comparing the new- with the old.
  **/
  rel **convert_global_1d_into_2d_rel(int node_id, uint arr_size) {
    rel *node_1d_list = NULL; uint node_1d_list_size = 0; uint *node_1d_map = NULL;
    if(get_node_1d_list_settings(node_id, node_1d_list, node_1d_list_size, node_1d_map, arr_size)) {
#ifndef NDEBUG
      loint offset = 0; for(int i=0; i<node_id; i++) {offset+=total_receive_size_for_the_nodes[i];}
      assert(node_1d_list == (global_1d_list + offset));
      assert(node_1d_list_size == (uint)total_receive_size_for_the_nodes[node_id]);
      assert(node_1d_map == global_map_list+(node_id*arr_size));
#endif
      return convert_1d_into_2d_rel(node_1d_list, node_1d_list_size, node_1d_map, arr_size);
    } else return NULL;
  }
  /**
     Converts the list, but does not deallocate data:
     @remarks Include tests, comparing the new- with the old.
  **/
  rel **convert_1d_into_2d_rel(rel *internal_list_1d, uint internal_list_1d_size, uint *map_internal_list, uint arr_size) {
    if(internal_list_1d && internal_list_1d_size && map_internal_list) {
      assert(map_internal_list);
      assert(internal_list_1d);
      assert(internal_list_1d_size);      
      rel *internal_list_1d_end = internal_list_1d + internal_list_1d_size;
      rel **temp = new rel*[arr_size];
      assert(temp);
      uint total_sum_of_lengths_added = 0;
      uint cnt_elements_added = 0;
      for(uint i = 0; i < arr_size; i++) {
	temp[i] = NULL;
	if(map_internal_list[i]) {
	  rel *current = internal_list_1d + total_sum_of_lengths_added;
	  assert(current < internal_list_1d_end);
	  temp[i] = new rel[map_internal_list[i]];
	  assert(temp[i]);
	  memcpy(temp[i], current, sizeof(rel)*map_internal_list[i]);
#ifndef NDEBUG
	  cnt_elements_added += get_size(temp[i]);
	  //! Compares with the 1d-list
	  for(uint out = 0; out < map_internal_list[i]; out++) {
	    assert(temp[i][out].ind_out == current[out].ind_out);
	    assert(temp[i][out].distance == current[out].distance);
	  }
#endif
	  total_sum_of_lengths_added += map_internal_list[i];
	}
      }
#ifndef NDEBUG
      assert(get_size(temp, arr_size) == cnt_elements_added);
      assert(id_simil_list::get_total_number_of_pairs(temp, arr_size) == get_total_number_of_pairs(temp, arr_size));
#endif
      return temp;
    } else return NULL; 
  }
  //! Translate the entire list into a 1d-representation.
  void translate_into_1d_list(id_simil_list *this_obj,rel_t **arr, uint arr_size) {
    if(!arr || !arr_size) {
#ifndef NDEBUG
      fprintf(stderr, "[myrank=%d]!!\tInput-putative-ortholog-list not set at line %d in file %s\n", myrank, __LINE__, __FILE__);
#endif
      return;
    }
    //! Fills the list:
    assert(total_receive_size_for_the_nodes);
    total_receive_size_for_the_nodes[myrank] = (int)get_size(arr, arr_size);
    
    //! Allocates- and initiates the 1d-list:
    assert(!internal_list_1d);
    assert(!internal_list_1d_size);
    //! Note: The size estimate also include the first object holding the metadata information.
    internal_list_1d_size = (uint)total_receive_size_for_the_nodes[myrank];
    loint temp_pos = 0;
    internal_list_1d = rel::init_list(internal_list_1d_size, temp_pos);
    rel *internal_list_1d_curr = internal_list_1d;
    rel *internal_list_1d_end = internal_list_1d + internal_list_1d_size;
    //! Allocates for the list holding the specific sizes for each protein (and initialises on-the-fly during the for-loop):
    assert(!map_internal_list);
    map_internal_list = new uint[arr_size];
    log_builder::test_memory_condition_and_if_not_abort(map_internal_list !=NULL, __LINE__, __FILE__, __FUNCTION__);
    uint total_sum_of_lengths_added = 0;
    for(uint i = 0; i < arr_size; i++) {
      map_internal_list[i] = get_size(arr[i]);
      uint temp_size = get_size(arr[i]);
      if(map_internal_list[i]) {
	assert(temp_size == map_internal_list[i]);
	assert((internal_list_1d_curr+map_internal_list[i]) <= internal_list_1d_end);
	memcpy(internal_list_1d_curr, arr[i], sizeof(rel)*map_internal_list[i]);
	internal_list_1d_curr += map_internal_list[i];
	assert(internal_list_1d_curr <= internal_list_1d_end);
#ifndef NDEBUG
	rel *current = internal_list_1d + total_sum_of_lengths_added;
	assert(current < internal_list_1d_end);
	for(uint out = 0; out < map_internal_list[i]; out++) {
	  assert(arr[i][out].ind_out == current[out].ind_out);
	  assert(arr[i][out].distance == current[out].distance);
	}
	total_sum_of_lengths_added += map_internal_list[i];
#endif
      }
    }

#ifndef NDEBUG
    
    if(total_receive_size_for_the_nodes[myrank]) { // Only validates if data is set:

      //! Validate that the number of sum_of_distances are the same in both containers (the temp and the real)
      float sum_of_distances_before = get_total_number_of_distances(internal_list_1d, internal_list_1d_size, map_internal_list, arr_size);
      float sum_of_distances_after  = this_obj->get_total_sum_of_distances_for_pairs();
      if((int)sum_of_distances_before != (int)sum_of_distances_after) {
	fprintf(stderr, "!!\t[myrank=%d]\t sum_of_distances_before(%f) and sum_of_distances_after(%f) at line %d in file %s\n", myrank, sum_of_distances_before, sum_of_distances_after, __LINE__, __FILE__);
	assert(sum_of_distances_before == sum_of_distances_after);
      }
      
      //! Compare the 1d-list with the given 
      total_sum_of_lengths_added = 0;
      for(uint i = 0; i < arr_size; i++) {
	if(map_internal_list[i]) {
	  rel *current = internal_list_1d + total_sum_of_lengths_added;
	  assert(current < internal_list_1d_end);
	  for(uint out = 0; out < map_internal_list[i]; out++) {
	    assert(arr[i][out].ind_out == current[out].ind_out);
	    assert(arr[i][out].distance == current[out].distance);
	  }
	  total_sum_of_lengths_added += map_internal_list[i];
	}
      }

      //! Do a back-translation, and validate that tehy are equal.
      rel **temp = convert_1d_into_2d_rel(internal_list_1d, internal_list_1d_size, map_internal_list, arr_size); 
      //! Compares with the original 2d-list (given as input)
      for(uint i = 0; i < arr_size; i++) {
	assert(get_size(temp[i]) == get_size(arr[i]));
	for(uint out = 0; out < get_size(temp[i]); out++) {
	  assert(temp[i][out].ind_out  == temp[i][out].ind_out);
	  assert(temp[i][out].distance == temp[i][out].distance);
	}
      }
  
      //! Validate that the number of relations are the same in both containers (the temp and the real)
      const uint relations_before = id_simil_list::get_total_number_of_pairs(arr, arr_size);
      const uint relations_after =  id_simil_list::get_total_number_of_pairs(temp, arr_size);
      if(relations_before != relations_after) {
	fprintf(stderr, "!!\t[myrank=%d]\t relations_before(%u) and relations_after(%u) at line %d in file %s\n", myrank, relations_before, relations_after, __LINE__, __FILE__);
	assert(relations_before == relations_after);
      }

      //! Validate that the number of sum_of_distances are the same in both containers (the temp and the real)
      sum_of_distances_before = id_simil_list::get_total_sum_of_distances_for_pairs(arr, arr_size);
      sum_of_distances_after =  id_simil_list::get_total_sum_of_distances_for_pairs(temp, arr_size);
      if(sum_of_distances_before != sum_of_distances_after) {
	fprintf(stderr, "!!\t[myrank=%d]\t sum_of_distances_before(%f) and sum_of_distances_after(%f) at line %d in file %s\n", myrank, sum_of_distances_before, sum_of_distances_after, __LINE__, __FILE__);
	assert(sum_of_distances_before == sum_of_distances_after);
      }

    
      //! De-allocates the temporary list:
      de_allocate(temp, arr_size);
    }
#endif
  }

  /**
     @brief Send metadata using “blocking communication”
     @remarks: Include tests validating the operation.
  **/
  void send_metadata_with_blocking(rel_t **arr, uint arr_size) {    
#ifndef NDEBUG
    //! The mpi sending procedure, summing the values using the same pointer for the task:
    assert(total_receive_size_for_the_nodes);
    uint my_size = (uint)total_receive_size_for_the_nodes[myrank];
    uint total_size_of_data = my_size;
    MPI_Allreduce(MPI_IN_PLACE, &total_size_of_data, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
    assert(total_size_of_data >= my_size);
      
#endif
    //! Sends- and receives the sizes of the lists:
    uint temp = (uint)total_receive_size_for_the_nodes[myrank];
    MPI_Allgather(&temp, 1, MPI_UNSIGNED,
		  total_receive_size_for_the_nodes, /*elements sent by each node is */1, MPI_UNSIGNED, MPI_COMM_WORLD);
#ifndef NDEBUG    
    uint sum = 0; for(int i = 0; i < number_of_nodes; i++) {sum+=(uint)total_receive_size_for_the_nodes[i];}
    if(sum != total_size_of_data) {
      printf("[myrank=%d]\tThe sending did not work: Total_sum_before_transfer=%u, and after_trasfer=%u, at line %d in file %s\n", myrank, total_size_of_data, sum, __LINE__, __FILE__); 
      assert(sum == total_size_of_data);
    }
#endif    
  }

  //! Deallocates the object given:
  void de_allocate(rel_t **&temp, uint arr_size) {
    if(temp) {
      for(uint i = 0; i < arr_size; i++) {
	if(temp[i]) {delete [] temp[i]; temp[i] = NULL;}
      }
      delete [] temp; temp = NULL;
    }
  }

  //! Send the internal list to all nodes, using "blocking communication"
  void tx_rx_pairs_lists_with_blocking(rel_t **arr, uint arr_size) {
    //! What we expect: If we have data, we expect both the container- and the size to be set:
    if(internal_list_1d || internal_list_1d_size) {
      assert(internal_list_1d_size);
      assert(internal_list_1d);
    }


    //! Sends- and receives the specifications for each pair:
    assert(!global_map_list);
    global_map_list = new uint[arr_size * number_of_nodes];
    log_builder::test_memory_condition_and_if_not_abort(global_map_list !=NULL, __LINE__, __FILE__, __FUNCTION__);
    memset(global_map_list, 0, sizeof(uint)*number_of_nodes*arr_size);
    MPI_Allgather(map_internal_list, arr_size, MPI_UNSIGNED,
		  global_map_list, /*elements sent by each node is */arr_size, MPI_UNSIGNED, MPI_COMM_WORLD);


    //! The lists should not have been initiated; then initiates them:
    assert(!global_1d_list);
    assert(!global_1d_list_size);
    global_1d_list_size = 0; for(int i = 0; i < number_of_nodes; i++) {global_1d_list_size+=(uint)total_receive_size_for_the_nodes[i];}
    assert(global_1d_list_size);
#ifndef NDEBUG
    //! Validates 'global_map_list': Get the total sum of lengths before- and after the 'merge':
    uint my_lengths_before = 0; for(uint i = 0; i < arr_size; i++) my_lengths_before += map_internal_list[i];
    uint lengths_before = my_lengths_before;
    MPI_Allreduce(MPI_IN_PLACE, &lengths_before, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
    assert(lengths_before >= my_lengths_before);
    uint lengths_after = 0; for(uint i = 0; i < (arr_size*number_of_nodes); i++) lengths_after += global_map_list[i];

    //! Validate that the number of lengths are the same in both containers (the temp and the real)
    if(lengths_before != lengths_after) {
      fprintf(stderr, "!!\t[myrank=%d]\t lengths_before(%u) and lengths_after(%u) at line %d in file %s\n", myrank, lengths_before, lengths_after, __LINE__, __FILE__);
      assert(lengths_before == lengths_after);
    }      
#endif

    //! Sends- and receives the sizes of the lists:
    loint temp_pos; global_1d_list = rel::init_list(global_1d_list_size, temp_pos);    
    //! Builds the specific datatype:
    MPI_Datatype mpi_type_rel; mpi_transfer_rel::mpi_build_datatype(mpi_type_rel);
    //    int global_1d_list_size_expected = global_1d_list_size;
    int recv_count[number_of_nodes]; memset(recv_count, 0, sizeof(int)*number_of_nodes);
    int displs[number_of_nodes]; memset(recv_count, 0, sizeof(int)*number_of_nodes);

    //! Sets the  recv_counts and displacements:
    for (int i=0; i<number_of_nodes; i++) {
      recv_count[i] = total_receive_size_for_the_nodes[i];
      if(i>0) {
	displs[i] = displs[i-1] + recv_count[i-1];
      } else displs[i] = 0;
    }

    MPI_Allgatherv(internal_list_1d, internal_list_1d_size,
		   mpi_type_rel, global_1d_list, recv_count,
		   displs,
		   //&global_1d_list_size_expected, 
		   mpi_type_rel, MPI_COMM_WORLD);
    MPI_Type_free(&mpi_type_rel);
   
    /*    
    if(true) { // Remove the contents of this block!
      static bool data_is_written = false; // Implying that only the first time data is written to the file:

      if(!data_is_written) {
	if(true) {
	  char cmd[100]; memset(cmd, '\0', 100); sprintf(cmd, "temp/before_%d_n%d.txt", myrank, number_of_nodes);
	  FILE *file = fopen(cmd, "w");
	  assert(file);
	  const uint length = internal_list_1d_size;
	  for(uint i = 0; i < length; i++) {
	    fprintf(file, "%d\t%d\t%f\n", i, internal_list_1d[i].ind_out, internal_list_1d[i].distance);
	  }
	  fclose(file);
	}
	char cmd[100]; memset(cmd, '\0', 100); sprintf(cmd, "temp/after_%d_n%d.txt", myrank, number_of_nodes);
	FILE *file = fopen(cmd, "w");
	assert(file);
	for(uint node_id = 0; node_id < (uint)number_of_nodes; node_id++) {
	  const uint length = displs[node_id] + recv_count[node_id];
	  for(uint i = displs[node_id]; i < length; i++) {
	    fprintf(file, "%d\t%d\t%f\n", i, global_1d_list[i].ind_out, global_1d_list[i].distance);
	  }
	}
	fclose(file);
      }

      data_is_written = true; 
    }
*/
#ifndef NDEBUG
    //! Get the total number of relations:
    const uint my_relations_before = id_simil_list::get_total_number_of_pairs(arr, arr_size);
    uint relations_before = my_relations_before;
    MPI_Allreduce(MPI_IN_PLACE, &relations_before, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
    assert(relations_before >= my_relations_before);
    
    //! Converts 'my object' into a 2d-list, and verifies it has the same number of relations as those myrank expect.:
    rel **temp = convert_1d_into_2d_rel(internal_list_1d, internal_list_1d_size, map_internal_list, arr_size);
    const uint my_relations_before_sent = id_simil_list::get_total_number_of_pairs(temp, arr_size);
    assert(my_relations_before_sent == my_relations_before);
    de_allocate(temp, arr_size);

    //! Verifies that the "merged object" gives the same result:
    temp = convert_global_1d_into_2d_rel(myrank, displs, arr_size);
    const uint my_relations_after_sent =  id_simil_list::get_total_number_of_pairs(temp, arr_size);
    assert(my_relations_after_sent == my_relations_before_sent);
    de_allocate(temp, arr_size);

    //! Get the total sum-of-distances:
    const float my_sum= id_simil_list::get_total_sum_of_distances_for_pairs(arr, arr_size);
    float sum_of_distances_before = my_sum;
    MPI_Allreduce(MPI_IN_PLACE, &sum_of_distances_before, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
    assert(sum_of_distances_before >= my_sum);


    //! Convert the received lists to 2d, and compare their
    uint relations_after =  0; 
    float sum_of_distances_after = 0; 
    rel_t *nodes_given_data = global_1d_list;
    for(int node_id = 0; node_id < number_of_nodes; node_id++) {
      const uint nodes_given_data_size = (uint)total_receive_size_for_the_nodes[node_id];
      //! Verify that the relative pouinnter corresponds to teh deplacement vector:
      assert(nodes_given_data == (global_1d_list + displs[node_id]));

      //! Update the container for current- and later validation:
      rel **temp = convert_global_1d_into_2d_rel(node_id, arr_size);
      //      rel **temp = convert_global_1d_into_2d_rel(node_id, displs, arr_size);
      relations_after        += id_simil_list::get_total_number_of_pairs(temp, arr_size);
      //! Validate that the number of relations are the same in both containers (the temp and the real)
      if(!(relations_before >= relations_after)) {
	fprintf(stderr, "!!\t[myrank=%d]\t relations_before(%u) and relations_after(%u) at line %d in file %s\n", myrank, relations_before, relations_after, __LINE__, __FILE__);
	assert(relations_before >= relations_after);
      }

      sum_of_distances_after += id_simil_list::get_total_sum_of_distances_for_pairs(temp, arr_size);
      assert((int)sum_of_distances_after <= (int)sum_of_distances_before);

      //! Updates the pointer for the next node (to be processed):
      nodes_given_data += nodes_given_data_size;

      //! Deallocates the temporary object:
      de_allocate(temp, arr_size);
    }

    //! Validate that the number of relations are the same in both containers (the temp and the real)
    if(relations_before != relations_after) {
      fprintf(stderr, "!!\t[myrank=%d]\t relations_before(%u) and relations_after(%u) at line %d in file %s\n", myrank, relations_before, relations_after, __LINE__, __FILE__);
      assert(relations_before == relations_after);
    }

    //! Validate that the number of sum_of_distances are the same in both containers (the temp and the real)
    if(sum_of_distances_before != sum_of_distances_after) {
      fprintf(stderr, "!!\t[myrank=%d]\t sum_of_distances_before(%f) and sum_of_distances_after(%f) at line %d in file %s\n", myrank, sum_of_distances_before, sum_of_distances_after, __LINE__, __FILE__);
      assert(sum_of_distances_before == sum_of_distances_after);
    }

    //! Merge the lists into a single object, validating that the new object consists of the sum of the earlier objects!
    class id_simil_list id_simil = id_simil_list(arr_size, NULL);
    assert(id_simil.has_data());
    nodes_given_data = global_1d_list;
    uint cnt_relations_inserted = 0;
    assert(0==id_simil.get_total_number_of_pairs()); // Should be empty.
    uint *global_map_list_curr = global_map_list;
    uint *global_map_list_end  = global_map_list + (number_of_nodes*arr_size);
    for(int node_id = 0; node_id < number_of_nodes; node_id++) {
      const uint nodes_given_data_size = (uint)total_receive_size_for_the_nodes[node_id];
      rel **temp = convert_global_1d_into_2d_rel(node_id, arr_size);
      //! Note: Not all nodes have data; arises when #(nodes) > #(taxa).
      if(temp) {
	for(uint i = 0; i < arr_size; i++) {
	  if(temp[i]) {
	    uint size = get_size(temp[i]); //id_simil_list::get_number_of_pairs(temp[i]);
	    assert(size == global_map_list_curr[i]);
	    rel_t *row =  temp[i]; 
	    //! Note: Adds '+1' as the length corresponds to the number of object, but we
	    //!       start on object '1' due to the first is a metadata-object.
	    for(uint out = 1; out < size; out++) {
	      id_simil.insertGlobalElement(i, row[out].ind_out, row[out].distance);
	      cnt_relations_inserted++;
	    }
	  }
	}
	de_allocate(temp, arr_size);
      }
      //! Mapping: Updates the pointerfor the next node (to be processed):
      global_map_list_curr += arr_size;
      assert(global_map_list_curr <= global_map_list_end);
      //! Content: Updates the pointer for the next node (to be processed):
      nodes_given_data += nodes_given_data_size;
    }

    //! Validate that the number of relations inserted corresponds the what the function call gives us:
    relations_after = id_simil.get_total_number_of_pairs();
    if(cnt_relations_inserted != relations_after) {
      fprintf(stderr, "!!\t[myrank=%d]\t cnt_relations_inserted(%u) and relations_after(%u) at line %d in file %s\n", myrank, cnt_relations_inserted, relations_after, __LINE__, __FILE__);
      assert(cnt_relations_inserted == relations_after);
    }

    assert(get_total_number_of_global_pairs(arr_size)==cnt_relations_inserted);
    //! Validate that the number of relations are the same in both containers (the temp and the real)
    if(relations_before != relations_after) {
      fprintf(stderr, "!!\t[myrank=%d]\t relations_before(%u) and relations_after(%u) at line %d in file %s\n", myrank, relations_before, relations_after, __LINE__, __FILE__);
      assert(relations_before == relations_after);
    }
    sum_of_distances_after = id_simil.get_total_sum_of_distances_for_pairs();
    //! Validate that the number of sum_of_distances are the same in both containers (the temp and the real)
    log_builder::compare_floats(sum_of_distances_before, sum_of_distances_after, relations_before, __LINE__, __FILE__, __FUNCTION__);
    id_simil.free_memory();   
#endif
  }

  //! Merge the received content into this array:
  void merge_received_content(id_simil_list *this_obj, uint total_number_of_proteins, uint &arr_size) {
    assert(this_obj);
    assert(global_map_list);
    if(!this_obj->has_data()) {
      arr_size = total_number_of_proteins;
      loint total = 0; for(int i=0; i<number_of_nodes; i++) {total+=total_receive_size_for_the_nodes[i];}
      const uint index_size = (total/arr_size)+1;
      this_obj->initGlobalArr(arr_size, index_size);
      //      id_simil_list::initLocalArr(arr, arr_size, index_size);
    }
    assert(arr_size);

#ifndef NDEBUG
    //! 
    //! Validate that the number of relations inserted corresponds the what the function call gives us:
    uint relations_before = get_total_number_of_pairs(internal_list_1d, internal_list_1d_size, map_internal_list, arr_size);
    uint relations_after = this_obj->get_total_number_of_pairs();
    if(relations_before != relations_after) {
      fprintf(stderr, "!!\t[myrank=%d]\t relations_before(%u) and relations_after(%u) at line %d in file %s\n", myrank, relations_before, relations_after, __LINE__, __FILE__);
      assert(relations_before == relations_after);
    }

    //! Validate that the number of sum_of_distances are the same in both containers (the temp and the real)
    float sum_of_distances_before = get_total_number_of_distances(internal_list_1d, internal_list_1d_size, map_internal_list, arr_size);
    float sum_of_distances_after  = this_obj->get_total_sum_of_distances_for_pairs();
    log_builder::compare_floats(sum_of_distances_before, sum_of_distances_after, relations_before, __LINE__, __FILE__, __FUNCTION__);

    //! Prepares the comparsion
    relations_before = get_total_number_of_global_pairs(arr_size);
    sum_of_distances_before = get_total_number_of_global_distances(arr_size);

#endif

    //! The setting of values into 'this':
    uint cnt_relations_inserted = 0;
    for(int node_id = 0; node_id < number_of_nodes; node_id++) {
      if(node_id != myrank) // If debug 
	{ // Does not add the relations of myrank, as they already resides in both the list and the debug counters:
	rel *node_1d_list = NULL; uint node_1d_list_size = 0; uint *node_1d_map = NULL;
	uint elements_inserted_this_node = 0;
	if(get_node_1d_list_settings(node_id, node_1d_list, node_1d_list_size, node_1d_map, arr_size)) {
	  rel *node_1d_list_end = node_1d_list + node_1d_list_size;
	  uint total_sum_of_lengths_added = 0;
	  uint found_sum_pairs = 0;
	  //	  uint cnt_elements_added = 0;
	  for(uint i = 0; i < arr_size; i++) {
	    if(node_1d_map[i]) {
	      rel *current = node_1d_list + total_sum_of_lengths_added;
	      assert(current < node_1d_list_end);
	      if(current) {
		const uint size = get_size(current); //id_simil_list::get_number_of_pairs(temp[i]);
		if(size) {
		  for(uint out = 1; out < size; out++) {
		    this_obj->insertGlobalElement(i, current[out].ind_out, current[out].distance);
#ifndef NDEBUG
		    sum_of_distances_after += current[out].distance;
		    relations_after++;
		    found_sum_pairs++;
#endif
		  }
		}
	      }
	      //! In order to place the pointer, we updates:
	      total_sum_of_lengths_added += node_1d_map[i];
	    }
	  }
#ifndef NDEBUG
	  const uint this_sum_pairs_actual = get_total_number_of_pairs(node_1d_list, node_1d_list_size, node_1d_map, arr_size);
	  assert(this_sum_pairs_actual == found_sum_pairs);
	  
#endif
	  //elements_inserted_this_node = get_total_number_of_pairs(node_1d_list, node_1d_list_size, node_1d_map, arr_size);;
	  cnt_relations_inserted += elements_inserted_this_node;

	}
      }
    }
#ifndef NDEBUG


    //! Validate that the number of relations are the same in both containers (the temp and the real)
    if(relations_before != relations_after) {
      fprintf(stderr, "!!\t[myrank=%d]\t relations_before(%u) and relations_after(%u) at line %d in file %s\n", myrank, relations_before, relations_after, __LINE__, __FILE__);
      assert(relations_before == relations_after);
    }

    //! Validate that the number of sum_of_distances are the same in both containers (the temp and the real)
    log_builder::compare_floats(sum_of_distances_before, sum_of_distances_after, relations_before, __LINE__, __FILE__, __FUNCTION__);

    //! The container we've set data into: Validating the number of relations are the same in both containers (the temp and the real)
    relations_after = this_obj->get_total_number_of_pairs();
    if(relations_before != relations_after) {
      fprintf(stderr, "!!\t[myrank=%d]\t relations_before(%u) and relations_after(%u) at line %d in file %s\n", myrank, relations_before, relations_after, __LINE__, __FILE__);
      assert(relations_before == relations_after);
    }
    //! The container we've set data into: Validating the number of sum_of_distances are the same in both containers (the temp and the real)
    sum_of_distances_after  = this_obj->get_total_sum_of_distances_for_pairs();
    log_builder::compare_floats(sum_of_distances_before, sum_of_distances_after, relations_before, __LINE__, __FILE__, __FUNCTION__);
#endif
  }
  //! Prints the list of responsibilities:
  void print_responsibility_list(bool *list_of_this_node_taxa_responsilibties, taxa* listTaxa, int taxon_length) {
    assert(listTaxa);
    assert(taxon_length);
    assert(list_of_this_node_taxa_responsilibties);
    for(int taxon_id = 0; taxon_id < taxon_length; taxon_id++) {
      if(list_of_this_node_taxa_responsilibties[taxon_id]) {
	printf("[taxon=%u] has as responsible myrank[%d]\n", taxon_id, myrank);
      }
    }
  }

  //! @return a global list of the elements given, except 'this':
  id_simil_list *build_id_simil_list(stack_rel *stackRel,bool *list_of_this_node_taxa_responsilibties, uint *&lst_taxa_lengths) {
    assert(stackRel);
    assert(listTaxa);
    assert(taxon_length);
    assert(list_of_this_node_taxa_responsilibties);
#ifndef NDEBUG
    //! Intialising, later to be used verifying that if a node has data about a taxon,
    //! it corresponds to the sum of of it found at all involved nodes:
    lst_taxa_lengths = new uint[taxon_length];
    log_builder::test_memory_condition_and_if_not_abort(lst_taxa_lengths !=NULL, __LINE__, __FILE__, __FUNCTION__);
    memset(lst_taxa_lengths, 0, sizeof(uint)*taxon_length);    
#endif
    //! Some counters:
    loint total_cnt_relations = 0, myrank_cnt_relations = 0, elements_inserted_into_temp_object = 0;
    const uint total_number_of_proteins = listTaxa[taxon_length-1].rel_end;
    assert(total_number_of_proteins);
    class id_simil_list *temp_list = new id_simil_list(total_number_of_proteins, listTaxa);
    for(int taxon_id = 0; taxon_id < taxon_length; taxon_id++) {
      if(list_of_this_node_taxa_responsilibties[taxon_id]) {
	//! Ignore, but if in debug count the relations:
#ifndef NDEBUG
	for(int world_index_in =listTaxa[taxon_id].rel_start; world_index_in< listTaxa[taxon_id].rel_end;world_index_in++) {
	  myrank_cnt_relations +=stackRel[world_index_in].unsafe_size();	
	  total_cnt_relations +=stackRel[world_index_in].unsafe_size();
	  lst_taxa_lengths[taxon_id] +=stackRel[world_index_in].unsafe_size();
	}
#endif
      } else {
	//! They do not belong to this node ('myrank'): We therefore pop- and send them away:
	for(int world_index_in =listTaxa[taxon_id].rel_start; world_index_in< listTaxa[taxon_id].rel_end;world_index_in++) {  
#ifndef NDEBUG
	  total_cnt_relations +=stackRel[world_index_in].unsafe_size();
	  elements_inserted_into_temp_object += stackRel[world_index_in].unsafe_size();
	  lst_taxa_lengths[taxon_id] +=stackRel[world_index_in].unsafe_size();
#endif
	  //! Insert them into the bunch to be sent:
	  o_rel_t out_pair;
	  while(stackRel[world_index_in].try_pop(out_pair)) {
	    const int world_out = out_pair.ind_out;
	    temp_list->insertGlobalElement(world_index_in, world_out, out_pair.distance);
	  }	
	}
      }
    }
#ifndef NDEBUG
    //! Verifies correct insertion in the temporary id_simil_list object:
    const loint total_relations = myrank_cnt_relations + temp_list->get_total_number_of_pairs();
    if(total_cnt_relations != total_relations) {
      printf("!!\t[myrank=%d]\t with co-orthologs(%lld) != (myrank_cnt_relations(%lld) + temp-sending-pairs(%u)), at line %d in file %s\n", myrank, total_cnt_relations, myrank_cnt_relations, temp_list->get_total_number_of_pairs(), __LINE__, __FILE__);
      assert(total_cnt_relations == total_relations);
    }
    assert(temp_list->get_total_number_of_pairs() == elements_inserted_into_temp_object);
#endif
    return temp_list;
  }

  //! @return the total number of relations for the given stackRel
  uint get_total_number_of_relations(stack_rel *stackRel) {
    assert(listTaxa);
    assert(taxon_length);
    uint total_cnt_relations = 0;
    for(int world_index_in =0; world_index_in< listTaxa[taxon_length-1].rel_end;world_index_in++) {
      total_cnt_relations += (uint)stackRel[world_index_in].unsafe_size();
    }
    return total_cnt_relations;
  }
  /**
     @brief Include the elements received into the stack:
     @param <stackRel> The element to insert data into
     @param <temp_obj> The element holding the data to be inserted into the list (ie, has the data we received).
     @param <list_of_this_node_taxa_responsilibties> To know what elements belong to us
     @param <lst_taxa_lengths> If in debug deallocates the list (used for verification).
     @remarks Verifying that if a node has data about a taxon, it corresponds to the sum of of it found at all involved nodes, ie we have the two cases: 
     # if a node has data about a taxon, it corresponds to the sum of of it found at all involved nodes and 
     # else, the taxon-ids elements should have removed all of those relations it does not bear responsibility for.
  **/
  void include_elements_received_into_stack(stack_rel *stackRel, id_simil_list *temp_obj, bool *list_of_this_node_taxa_responsilibties, uint *&lst_taxa_lengths) {
    assert(listTaxa);
    assert(stackRel);
    assert(taxon_length);
    assert(list_of_this_node_taxa_responsilibties);
    assert(temp_obj);
#ifndef NDEBUG
    assert(lst_taxa_lengths);
    //! Intialising, later to be used verifying that if a node has data about a taxon,
    //! it corresponds to the sum of of it found at all involved nodes:
    uint *lst_after_taxa_lengths = new uint[taxon_length];
    log_builder::test_memory_condition_and_if_not_abort(lst_after_taxa_lengths !=NULL, __LINE__, __FILE__, __FUNCTION__);
    memset(lst_after_taxa_lengths, 0, sizeof(uint)*taxon_length);    
#else
    assert(!lst_taxa_lengths);
#endif

    const uint total_number_of_proteins = listTaxa[taxon_length-1].rel_end;
    assert(total_number_of_proteins);
    //    class id_simil_list *temp_list = new id_simil_list(total_number_of_proteins, listTaxa);
    int inserted_cnt_elements = 0;
    for(int taxon_id = 0; taxon_id < taxon_length; taxon_id++) {
      if(list_of_this_node_taxa_responsilibties[taxon_id]) {
	for(int world_index_in =listTaxa[taxon_id].rel_start; world_index_in< listTaxa[taxon_id].rel_end;world_index_in++) {  
	  //! Using the object is safe, as we did not insert any elements at our position:
	  uint size = 0; rel_t* row = temp_obj->getGlobalRow(world_index_in, size);
	  if(size && (row != NULL)) {
	    for(uint out = 0; out< size; out++) {
	      stackRel[world_index_in].push(o_rel(row[out].ind_out, row[out].distance));
	      inserted_cnt_elements++;
	    }
	  }
#ifndef NDEBUG
	  lst_after_taxa_lengths[taxon_id] +=stackRel[world_index_in].unsafe_size();
#endif
	}
      }
    }
#ifndef NDEBUG
    //! Gets an updated sum of the taxa-collection:
    //! Verifies that if a node has data about a taxon, it corresponds to the sum of of it found at all involved nodes:    
    
    //! Receives the global list:
    MPI_Allreduce(MPI_IN_PLACE, lst_taxa_lengths, taxon_length, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
    //! If a node has responisiblity for a taxon, the number of elements in "lst_after_taxa_lengths' should correspond to the sum received:
    for(int taxon_id = 0; taxon_id < taxon_length; taxon_id++) {
      if(list_of_this_node_taxa_responsilibties[taxon_id]) {
	//! The total number for this taxon should correspond to those residing in this node:
	assert(lst_after_taxa_lengths[taxon_id] == lst_taxa_lengths[taxon_id]);
      } else {
	//! This elements should have removed all of those relations it does not bear responsibility for:
	assert(lst_after_taxa_lengths[taxon_id] == 0);
      }
    }



    //! De-allocates the list used:
    if(lst_after_taxa_lengths) {delete [] lst_after_taxa_lengths; lst_after_taxa_lengths = NULL;}
    if(lst_taxa_lengths) {delete [] lst_taxa_lengths; lst_taxa_lengths = NULL;}
#endif
  }
 public:

  //! Perform the the exchange, letting all nodes get a complete list of the orthologs.
  void send_co_orthologs(taxa *_listTaxa, int _taxon_length, stack_rel *stackRel, list_file_parse<rel> *listStructData) {
    //! Set some global variables    
    listTaxa = _listTaxa;    taxon_length = _taxon_length;
    assert(listStructData);
    bool *list_of_this_node_taxa_responsilibties = listStructData->get_list_of_nodes_taxa_responsilibties();
    assert(list_of_this_node_taxa_responsilibties);
    //! Builds a global list of the elements given, except 'this':
    uint *lst_taxa_lengths = NULL;
    id_simil_list *temp_obj = build_id_simil_list(stackRel, list_of_this_node_taxa_responsilibties, lst_taxa_lengths);
    assert(temp_obj);
    uint arr_size=0; rel_t **arr=temp_obj->get_internal_list(arr_size);
    tx_rx_before_recip(temp_obj, arr, arr_size, arr_size);
    //! Include the elements received into the stack; if in debug deallocates the list (used for verification).
    include_elements_received_into_stack(stackRel, temp_obj, list_of_this_node_taxa_responsilibties, lst_taxa_lengths);
    if(temp_obj) {temp_obj->free_memory(); delete temp_obj; temp_obj = NULL;}
  }

  //! Perform the the exchange, letting all nodes get a complete list of the orthologs.
  void tx_rx_before_recip(id_simil_list *this_obj, rel_t **&arr, uint &arr_size, uint total_number_of_proteins) {
    //! Makes the list tight by remvoing information about empty space:
    compress_list(arr, arr_size);
    //! Translate the entire list into a 1d-representation.
    translate_into_1d_list(this_obj, arr, arr_size);

    //! Send metadata using “blocking communication”
    send_metadata_with_blocking(arr, arr_size);
    
    //! Send the internal list to all nodes, using "blocking communication"
    tx_rx_pairs_lists_with_blocking(arr, arr_size);

    //! Merge the received content into this array:
    merge_received_content(this_obj, total_number_of_proteins, arr_size);
    //    merge_received_content(this_obj, arr, arr_size, total_number_of_proteins);
  }

  //! De-allocates the memory reserved.
  void free_memory() {
    if(total_receive_size_for_the_nodes) {
      delete [] total_receive_size_for_the_nodes;
      total_receive_size_for_the_nodes = NULL;
    }
    if(internal_list_1d) {
      delete [] internal_list_1d;
      internal_list_1d = NULL;
      internal_list_1d_size = 0;
    }
    if(map_internal_list) {
      delete [] map_internal_list;
      map_internal_list = NULL;
    }
    if(global_1d_list) {
      delete [] global_1d_list;
      global_1d_list = NULL;
      global_1d_list_size = 0;
    }
    if(global_map_list) {
      delete [] global_map_list;
      global_map_list = NULL;
    }
  }
  //! The constructor:
  mpi_id_simil_list() :
    myrank(0), number_of_nodes(0), total_receive_size_for_the_nodes(NULL),
    internal_list_1d_size(0), internal_list_1d(NULL), map_internal_list(NULL), global_1d_list(NULL),
    global_1d_list_size(0), global_map_list(NULL), listTaxa(NULL), taxon_length(0)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
      MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes);
      assert(number_of_nodes);
      total_receive_size_for_the_nodes = new int[number_of_nodes];
      log_builder::test_memory_condition_and_if_not_abort(total_receive_size_for_the_nodes!=NULL, __LINE__, __FILE__, __FUNCTION__);
      memset(total_receive_size_for_the_nodes, 0, sizeof(int)*number_of_nodes);
    }; 
    
};

#endif
#endif
