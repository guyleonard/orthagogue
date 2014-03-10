#ifndef relations_list_h
#define relations_list_h
/**
   @file
   @brief Includes a set of rel objects.
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
v * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with orthAgogue
 . If not, see <http://www.gnu.org/licenses/>.
*/

#include "types.h"
#include "p_rel.h"
#include "rel.h"
#include "index_list.h"
#include "taxa.h"
#include "basket_parse.h"
/**
   @class relations_list
   @ingroup blastfile_container 
   @brief Includes a set of rel objects.
   @tparam <T> A rel class.
   @remarks  Include methods handling insertioon, copying, merging. etc.
   - Holds the relations both those in memory (the list) and in a binary file.
   @author Ole Kristian Ekseth (oekseth)
**/
template<class T>
class relations_list {
 private:
  //! The temprary buffer (storing the file-data) holding in-place data.
  T *list;  //  p_rel_t *buffer; 
  //! The number of objects written to the file.
  loint file_cnt_structs; 
  //! The position reached in the buffer while adding objects
  loint buffer_in_mem_pos; 
  //! The reserved size (in memmory) of the buffer
  loint buffer_in_mem_size; 
  //! The name of the file
  char *file_name; 
  //! The block number whom the parse belongs to; used in order to merge different sets when parallel prosessing is in use.
  loint block_number; 
  // Used to build the file name.
  int my_id;   int taxon_in, taxon_out;
  //! The maximal size of chars to store in memory: To be set at maximum if data are processed in parallel to avoid complexities before the mergin operation occurs.
  loint MAX_BUFFER_IN_MEM_SIZE;
  bool USE_BEST_BLAST_PAIR_SCORE; // If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
 public:
/**
   @brief Inserts a new list of type T
   @remarks Assumes the internal list is empty.
**/
void insert_new_list(T *arr_new, const uint arr_new_size) {
  if(arr_new_size) {
    assert(buffer_in_mem_pos < 1); // No data inserted.
    //! Ensures the lsit is large enough:
    if(arr_new_size > buffer_in_mem_size) {
      T::enlarge(list, buffer_in_mem_size, arr_new_size, buffer_in_mem_pos);
    }
    memcpy(list, arr_new, sizeof(T)*arr_new_size);
    buffer_in_mem_pos = arr_new_size;
    buffer_in_mem_size = arr_new_size;
  }
}


  //! Prints info about the object
  void print_object_info(FILE *f) {
    assert(f);
    fprintf(f, "relations_list_object{file_cnt_structs(%d), buffer_in_mem_pos(%d), buffer_in_mem_size(%d), buffer_is_set(yes=1==%d)\n",
	    (int)file_cnt_structs, (int)buffer_in_mem_pos, (int)buffer_in_mem_size, (list != NULL));
  }
  // -----------------------GET-methods-------------------------

  //! @return the buffer of T objects.
  T *get_buffer(){return list;}
  //! @return the buffer of T objects.
  T *get_list(){return list;}
  //! @return the length of data stored in the binary file.
  loint get_file_length(){return file_cnt_structs;}
  //! @return the length of data stored in the binary file.
  loint get_file_cnt_structs(){return file_cnt_structs;}
  //! @return the number of objects found in the object list.
  loint get_buffer_in_mem_pos(){return buffer_in_mem_pos;}
  //! @return the size set for the list storint the T objects.
  loint get_buffer_in_mem_size(){return buffer_in_mem_size;}
  //! @return the file name used for the binary files.
  char *get_file_name(){return file_name;}
  //! @return the block id for this object.
  loint get_block_number(){return block_number;}
  //! @return the id for this container.
  int get_my_id(){return my_id;}
  //! @return the id for the inner taxon set for this object.
  int get_taxon_in(){return taxon_in;}
  //! @return the id for the outer taxon set for this object.
  int get_taxon_out(){return taxon_out;}
  //! return the limit for the number of objects to be stored in the object before the object starts dumping these to the binary file.
  loint GET_MAX_BUFFER_IN_MEM_SIZE(){return GET_MAX_BUFFER_IN_MEM_SIZE;}
  //! @return the absolute position of the last object in the list (the position when the (potential) binary file are also included.)
  loint get_imaginary_next_position_in_file() {
    return file_cnt_structs + buffer_in_mem_pos;
  }
  //! @return the sum of 'file_cnt_structs+buffer_in_mem_pos' 
  loint getTotalLengthOfData() {
    return (file_cnt_structs + buffer_in_mem_pos);
  }


  //! @return the sum of 'file_cnt_structs+buffer_in_mem_pos' 
  float get_sum_of_pair_distances_in_memory() {
    return get_sum_of_pair_distances_in_memory(list, buffer_in_mem_pos);
    /*
    if(list) {
      float sum = 0;
      for(uint i = 0; i < buffer_in_mem_pos; i++) {
	sum += list[i].get_distance();
      }
      return sum;
    } else return 0.0;
*/
  }

  //! @return the sum of 'file_cnt_structs+buffer_in_mem_pos' 
  static float get_sum_of_pair_distances_in_memory(T *list, uint buffer_in_mem_pos) {
    if(list) {
      float sum = 0;
      for(uint i = 0; i < buffer_in_mem_pos; i++) {
	sum += list[i].get_distance();
      }
      return sum;
    } else return 0.0;
  }

  //! @return the file content only if availability of enough memory
  T *get_file_content() { 
    T * buff = T::reserve_list(file_cnt_structs);
    return T::read_file(file_name, buff, file_cnt_structs);
  }


  //! @return the content of data written to the file
  T *get_file_content(T *&buff) {
    return T::read_file(file_name, buff, file_cnt_structs);
  }
  //! @return the index out of the object list for the given param.
  uint get_ind_out(loint index) {
    assert(index < buffer_in_mem_size);
    return list[index].get_ind_out();
  }

  //! @return the elment given the local index.
  T get_element_at_local_index(loint index) {
    if(list && (index < buffer_in_mem_size)) return list[index];
    return T();
  }
  /**
     @brief Uses a T object, thereby calling this method for another type causes error!
     @return a basket_parse object the buffer assuming it's in memory,
     @remarks Returning object used for transfering data between pipes.
  **/
  class basket_parse *get_buffer_basket();
  

  /**
     @brief Retrieves data from an intermediate (binary) file.
     @return list of these given the input.
  **/
  T* get_file_block(char *file_name, uint start_pos, uint length) {
    FILE *file = fopen(file_name, "rb");
    if (file != NULL) {
      T *buffer_ret =  NULL; //tbb::tbb_allocator<T>().allocate(length);
      try {buffer_ret = tbb::tbb_allocator<T>().allocate(length);} 
      catch (std::exception& ba) {
	if(!log_builder::catch_memory_exeception(length, __FUNCTION__, __FILE__, __LINE__)) {
	  fprintf(stderr, "!!\t An interesting error was discovered: %s."
		  "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		  "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
	}
      }
      if(buffer_ret != NULL) {
	fseek(file, sizeof(T)*start_pos, SEEK_SET);
	fread(buffer_ret, sizeof(T), length, file); // to avoid the intial char
	fclose(file);
	return buffer_ret;
      } else { fprintf(stderr, "!!\tCould not allocate %u types of class T; Aborts!\n", length);
	exit(2);
      }
    } else {
      fprintf(stderr, "!!\tCould not open %s; Aborts!\n", file_name);
      exit(2);
    }
  }

  /**
     @brief Gets the partial buffer, given the local protein index
     @return The buffer for the proteins of type T
  **/
  T *getBuffer(loint index_start, loint length_buffer) {    
    if(length_buffer>0) { // has data
      const loint end_index = index_start + length_buffer;          
      if(list && (end_index <= buffer_in_mem_size)) {
	T *buffer_ret = list + index_start;	
	return buffer_ret;
      } else {
	fprintf(stderr, "!!\tData was not found in the buffer for the given protein(end_index=%lld). Please forward the following information to the developer at oekseth@gmail.com -- index_start(%lld), length_buffer(%lld), buffer_in_mem_pos(%lld), buffer_in_mem_size(%lld), file_cnt_structs(%lld), with message generated at line %d in file %s for method %s\n", end_index, index_start, length_buffer, buffer_in_mem_pos, buffer_in_mem_size, file_cnt_structs, __LINE__, __FILE__, __FUNCTION__);
	assert(false);
      }
      //      if(list == NULL) buffer_ret=T::get_file_block(file_name, index[protein_in].get_start_pos(), length_buffer);
    }
    return NULL;
  }
  // -----------------------SET-methods-------------------------

  /**
     @brief Sets the file name for the temporary (intermediate) files
     @remarks Used in the some of the constructors.
     - The file is removed when cleaning this class: the file is therefore not found if memory is cleaned, or the program to do not abort its exectuion.
  **/
  char *set_file_name(uint taxon_in, uint taxon_out, char *FILE_BINARY_LOCATION) {
    char *id_string = T::get_class_name();
    char *f_name = new char[1270];
    assert(f_name);
    memset(f_name, '\0', 1270);
    if(FILE_BINARY_LOCATION == NULL) FILE_BINARY_LOCATION = "./";
    if(FILE_BINARY_LOCATION != NULL && strlen(FILE_BINARY_LOCATION)>2) {
      if(FILE_BINARY_LOCATION[strlen(FILE_BINARY_LOCATION)-1] != '/')
	sprintf(f_name, "%s/%d_%d.%s", FILE_BINARY_LOCATION, taxon_in, taxon_out, id_string);
      else       sprintf(f_name, "%s%d_%d.%s", FILE_BINARY_LOCATION, taxon_in, taxon_out, id_string);
    } else sprintf(f_name, "%d_%d.%s", taxon_in, taxon_out, id_string);
    return f_name;
  }
  
  //! Prints the buffer
  void print_buffer(FILE *f) {
    assert(f);
    if(buffer_in_mem_pos) {
      fprintf(f, "-\tThe buffer of length %llu, with interval [%llu, %llu], is:\n", buffer_in_mem_pos, file_cnt_structs, file_cnt_structs+buffer_in_mem_pos);
      for(uint i = 0; i < buffer_in_mem_pos; i++) {
	fprintf(f, "buffer[%u] has protein(%u)\n", i, (uint)list[i].get_ind_out());
      }
    } 
  }

  /**
     @brief Prints the buffer as an abc file
     @remarks Understanding the internal data strcture, this mehod provides a good glance.
  **/
  void print_buffer_as_abc_files(FILE *f, index_list_t *listIndex, taxa_t *listTaxa, uint taxon_in, uint taxon_out) {
    assert(f);
    assert(listIndex);
    if(buffer_in_mem_pos) {
      for(int protein_in = 0; protein_in < listIndex->get_index_used()+1; protein_in++) {
	char *key_in = NULL;
	if(listTaxa && listTaxa[taxon_in].is_a_local_key(protein_in)) key_in = listTaxa[taxon_in].getCompleteProteinName(protein_in);
	
	if(listIndex->has_data(protein_in)) {
	  for(uint out_r = listIndex->get_index_start(protein_in); out_r < listIndex->get_end_pos(protein_in); out_r++) {
	    if(listTaxa) {
	      char *key_out = listTaxa[taxon_out].getCompleteProteinName(list[out_r].get_ind_out());
	      fprintf(f, "%s\t%s\t%f\n", key_in, key_out, list[out_r].get_distance());
	      delete [] key_out; key_out = NULL;
	    } else fprintf(f, "[%u]\t\t%d\t%d\t%f\t%u\n", out_r, protein_in, list[out_r].get_ind_out(), list[out_r].get_distance(), list[out_r].get_cnt_max());
	  }
	}
	delete [] key_in; key_in = NULL;
      }
    } 
  }

  /**
     @param <arg> the object to compare this with.
     @param <print_info> if set to true prints info about in-equalities.
     @return true if the input object is equal this.
  **/
  bool is_equal(relations_list *arg, bool print_info) {
    bool is_equal = true;
    if(arg) {
      if(buffer_in_mem_pos != arg->get_buffer_in_mem_pos()) {
	if(print_info) printf("%s not set (%d != %d), ", "mem_pos", (int)buffer_in_mem_pos, (int)arg->get_buffer_in_mem_pos());
	is_equal= false;
      }
      if(list && !(arg->get_buffer())) {
	if(print_info) printf("%s not set, ", "buffer");
	is_equal= false;
      }
      if(buffer_in_mem_pos != arg->get_buffer_in_mem_pos()) {
	if(print_info) printf("%s not set, ", "buffer_in_mem_pos");
	is_equal= false;
      }
/*       if(buffer_in_mem_size != arg->get_buffer_in_mem_size()) { */
/* 	if(print_info) printf("%s not set, ", "buffer_in_mem_size"); */
/* 	is_equal= false; */
/*       } */
/*       if(file_name != arg->get_file_name()) { */
/* 	if(print_info) printf("%s not set, ", "file_name"); */
/* 	is_equal= false; */
/*       } */
      if(file_cnt_structs != arg->get_file_cnt_structs()) {
	if(print_info) printf("%s not set, ", "file_cnt_structs");
	is_equal= false;
      }

      if(block_number != arg->get_block_number()) {
	if(print_info) printf("%s not set, ", "block_number");
	is_equal= false;
      }
      // Taxon data:
      if(taxon_in != arg->get_taxon_in()) {
	if(print_info) printf("%s not set, ", "taxon_in");
	is_equal= false;
      }
      if(taxon_out != arg->get_taxon_out()) {
	if(print_info) printf("%s not set, ", "taxon_out");
	is_equal= false;
      }
    }
    return is_equal;
  }

  //! Prints the contetn of the file
  void print_file() {
    T *data = get_file_content();
    for(uint i = 0; i < file_cnt_structs; i++) data[i].print();
    delete [] data; data = NULL;
  }


  //! Dumps the buffer to a file, not resetting the data
  void dump_buffer();

  /**
     @brief Dumps the file to a buffer if the buffer is less than the memory limit set
     else upload the file to the buffer
     @param <max_buffer_size> If above dumps to a binary file.
     @return true if memory is freed.
     @remarks 
     - Launched At the start of the processing of an operation specified  operation, but before the 'inparalog' operation 
     Description of the algorithm: Tests if the total size is less than the constraints set on memory usage (i.e. the input argument).
     - yes: Tests if it has any data stored in a file (i.e. the buffer is empty):
     - yes: Opens the file, and puts the file into the buffer
     - no:  Tests if there are any data in the buffer (buffer_in_mem_pos != 0)
     - yes: Inserts the buffer in the front of the file (i.e. in
     order for the indexes to maintain their correct position pointers
     - no:  Frees the reserved memory, and sets the variables to '0' and 'NULL' 
     - no:  Writes the buffer, if any buffer, to a file
     @author Ole Kristian Ekseth (oekseth)
  */
  bool dumpBufferOrFile(uint max_buffer_size) {
    const uint this_size = file_cnt_structs + buffer_in_mem_pos;
    bool memory_is_freed = false;
    if(this_size < max_buffer_size) {  
      // Tests if it has any data stored in a file (i.e. the list is empty):
      if(file_cnt_structs > 0) { // has data in memory: opens 
	if(list == NULL) {
	  assert(buffer_in_mem_size == 0);
	  assert(buffer_in_mem_pos == 0);
	  // First reserves memory and then reads the file into it
	  list = T::reserve_list(file_cnt_structs);
	  list =  T::read_file(file_name, list, file_cnt_structs);
	  buffer_in_mem_size = file_cnt_structs; buffer_in_mem_pos = file_cnt_structs; file_cnt_structs = 0;
	} else {// Already has data: copies the data from the file and merges it with the exsisting data:
	  assert(buffer_in_mem_size > 0);
	  if(this_size > buffer_in_mem_size) { // The list is not large enough
	    T::enlarge_and_put_current_at_back(list, this_size,  file_cnt_structs, buffer_in_mem_pos, USE_BEST_BLAST_PAIR_SCORE);
	    buffer_in_mem_size = this_size;
	  } else {
	    // Reads the file iin the front of the list, therefore a move of the exsisting items are mandatory:
	    memmove(list+file_cnt_structs, list, sizeof(T)*buffer_in_mem_pos);
	  }
	  list =  T::read_file(file_name, list, file_cnt_structs);
	  buffer_in_mem_pos = this_size; file_cnt_structs = 0;
	  buffer_in_mem_pos = this_size; // TODO: verify this is correct.
	}
      } else if (buffer_in_mem_pos == 0) { // has no data in memory: frees the reserved memory
	memory_is_freed = true;
	buffer_in_mem_pos = 0, buffer_in_mem_size = 0;
      }
    } else {
      dump_buffer(); // writes the data to the file
    }
    if(memory_is_freed) {T::close(list, buffer_in_mem_pos, buffer_in_mem_size);}
    return memory_is_freed;
  }

  /**
     @brief Merges an input object with this, only updating the object 'this' internal object.
     @attention The argument-data is placed at the end of 'this' container.
     @remarks Included in the code inline tests deactivated when the macro NDEBUG is set.
     @author Ole Kristian Ekseth (oekseth).
  **/
  void merge_buffers(relations_list *&fs, uint &cnt_elements_in_all_relation_lists) {  
    assert(file_name);

    const bool write_all_relations_to_file = false; // FIXME: set the intial 'value' to false, and then make 'this' accessible from the command line.

    uint debug_cnt_elements_written_to_file = 0;

#ifndef NDEBUG
      const uint this_length_old = getTotalLengthOfData();
      uint data_length = 0; if(fs) {data_length = fs->getTotalLengthOfData();}
      const uint tot_elements = this_length_old + data_length;
      float data_sum_of_distances = 0; if(fs) {data_sum_of_distances = fs->get_sum_of_pair_distances_in_memory();}
      const float sum_before_merging = get_sum_of_pair_distances_in_memory();
      float sum_before = 0; for(loint i = 0; i < buffer_in_mem_pos; i++) {sum_before += list[i].get_distance();}
      assert(sum_before == sum_before_merging);
      //! Verfies that the argument buffers is correctly understood:
      float sum_arg_before = 0;
      if(fs){T *fs_list=fs->get_buffer();for(loint i=0;i < fs->get_buffer_in_mem_pos(); i++) {sum_arg_before += fs_list[i].get_distance();}}
      assert(data_sum_of_distances == sum_arg_before);
#endif
    if(fs && fs->get_buffer_in_mem_pos()) { // Merging requires some data to merge into:
      if(!(fs->get_file_cnt_structs())) {
	const loint new_size = buffer_in_mem_pos + fs->get_buffer_in_mem_pos();
	cnt_elements_in_all_relation_lists += new_size;
	if(
	   write_all_relations_to_file || (
					   //	   (new_size > MAX_BUFFER_IN_MEM_SIZE) && 
					   (cnt_elements_in_all_relation_lists > MAX_BUFFER_IN_MEM_SIZE) && 
					   //! We avoid writing 'too small' data-sets into files, as we regard the overhead as too high.
					   (new_size > 1024)  // todo: validate the optimality of the '1024' threshold:
					   )
	   ) { // Then put both into the binary file:
	  debug_cnt_elements_written_to_file = fs->get_buffer_in_mem_pos();
	  T *buffer_one = NULL;
	  if(list) buffer_one = list; //->get_buffer();
	  T::write_buffers(file_name, buffer_one, buffer_in_mem_pos, fs->get_buffer(), fs->get_buffer_in_mem_pos(), file_cnt_structs);
	  file_cnt_structs += (buffer_in_mem_pos + fs->get_buffer_in_mem_pos());
	  init_T_data();       buffer_in_mem_pos = 0; // Empty, as the data has been written to a file:
	  relations_list<T>::close(fs, false);
	} else {
#ifndef NDEBUG
	  T *root_list_copy = T::reserve_list(buffer_in_mem_pos);
	  assert(root_list_copy);
	  loint root_list_copy_size = buffer_in_mem_pos;
	  for(uint i = 0; i < (uint)buffer_in_mem_pos; i++) root_list_copy[i].set_data(list[i]);
#endif
	  T::enlarge(list, buffer_in_mem_size, new_size, buffer_in_mem_pos);
	  assert(list);
	  if(buffer_in_mem_pos) {
	    memcpy(list+buffer_in_mem_pos, fs->get_buffer(), sizeof(T)*fs->get_buffer_in_mem_pos());
	  } else { // Enough taking all the data from the argument, as 'this' is not set.
	    // TODO: Consider grabbing the pointer instead.
	    memcpy(list, fs->get_buffer(), sizeof(T)*fs->get_buffer_in_mem_pos());
	  }
#ifndef NDEBUG
	  //! Note: Extensive tests written August 08 2012 by oekseth due to problems concluded to originate in rounding of floats.
	  //! Verifies that the old part is still the same:
	  float sum_after_this = 0;
	  if(buffer_in_mem_pos) {
	    for(loint i = 0; i < buffer_in_mem_pos; i++) {sum_after_this += list[i].get_distance();}
	  }
	  assert(sum_after_this == sum_before);	    
	  float sum_arg_aft = 0;
	  if(fs) {
	    // Verifies that the copied part is completely equal to the orginal we got:
	    T * arg = fs->get_buffer();
	    for(loint i = 0; i < fs->get_buffer_in_mem_pos(); i++) {
	      assert(list[i+buffer_in_mem_pos].get_distance() == arg[i].get_distance());
	      sum_arg_aft += arg[i].get_distance();
	    }
	  }
	  assert(sum_arg_aft == sum_arg_before);
	  assert(sum_arg_aft == data_sum_of_distances);

	  //! Verfies that the sum of inserts corresponds to our expectations:
	  assert(fs); // TODO: Unsure it this assumption holds
	  const loint new_size_tmp = buffer_in_mem_pos + fs->get_buffer_in_mem_pos();
	  float sum_real = 0;
	  float sum_real_this = 0, sum_real_arg = 0;
	  for(loint i = 0; i < new_size_tmp; i++) {
	    sum_real += list[i].get_distance();
	    if(i < buffer_in_mem_pos)  {
	      sum_real_this += list[i].get_distance();
	      assert(sum_real == sum_real_this);
	    }
	    else {
	      sum_real_arg += list[i].get_distance();
	      log_builder::compare_floats(sum_real, (sum_real_this + sum_real_arg), i, __LINE__, __FILE__, __FUNCTION__);
	      //! Downscales the floating points to avoid error when the difference is due to rounding:
	    }
	  }
	  assert(sum_real_this == sum_after_this);
	  assert(sum_real_arg == sum_arg_aft);
	  const float sum_discreet = sum_after_this + sum_arg_aft;
	  log_builder::compare_floats(sum_real, sum_discreet, new_size_tmp, __LINE__, __FILE__, __FUNCTION__);
	  //! Verfiies that the sum before- and after merging is correct:
	  assert(sum_before == sum_before_merging);
	  const float inserted_sum = sum_before_merging + data_sum_of_distances;
	  log_builder::compare_floats(sum_real, inserted_sum, new_size_tmp, __LINE__, __FILE__, __FUNCTION__);

	  // Performs inline validation that everything goes as it should. Hadny if sofwtare crashes, then reproducing it will be easier, simple by asking the user running the sofwtware in debug-mode, ie with undef NDEBUG
	  bool all_are_equal = true;
	  for(loint i = 0; i < buffer_in_mem_pos; i++) {	    
	    if(!list[i].equals(root_list_copy[i])) {
	      fprintf(stderr,"[%u]\tcopy{index_out(%u)} != list{index_out(%u))}\n",
		      (uint)i, root_list_copy[i].get_ind_out(), list[i].get_ind_out());
	      all_are_equal = false;
	      assert(all_are_equal);
	    }
	  }
	  assert(all_are_equal);
	  T * arg = fs->get_buffer();
	  for(loint i = 0; i < fs->get_buffer_in_mem_pos(); i++) {
	    if(!list[i+buffer_in_mem_pos].equals(arg[i])) {
	      fprintf(stderr,"[%u]\targument{index_out(%u)} != list{index_out(%u))}\n",
		      (uint)i, arg[i].get_ind_out(), list[i].get_ind_out());
	      all_are_equal = false;
	      assert(all_are_equal);
	    }
	  }
	  assert(all_are_equal);
	  loint temp1 = 0;
	  T::close(root_list_copy,temp1, root_list_copy_size);
#endif
	  buffer_in_mem_size = new_size;
	  buffer_in_mem_pos += fs->get_buffer_in_mem_pos();
	}
      } else {
	fprintf(stderr, "!!\tTried merging a structure of type relations_list into 'this' having data written to a file. As the files written are not unique with regard to id's, this will not work. If this functionality is of need, contact the developer at oekseth@gmail.com. This error is generated at line %d in file %s found in method %s\n", __LINE__, __FILE__, __FUNCTION__);
	assert(false);
      }
    }

#ifndef NDEBUG
      //! Validate that the total number of elements before corresponds to the sum after:
      assert(tot_elements == getTotalLengthOfData()); // The 'root' should now consist of the sum of elements.

      
      //! Get the actual sum manually, and compare it with the calculated:
      const float expected_sum_of_distances_after_merge = get_sum_of_pair_distances_in_memory();
      float sum_after = 0; for(loint i = 0; i < buffer_in_mem_pos; i++) {sum_after += list[i].get_distance();}
      assert(sum_after == expected_sum_of_distances_after_merge);
      
      //! Verfies that the argument buffers is correctly understood:
      float sum_arg = 0;
      if(fs) {T *fs_list=fs->get_buffer(); for(loint i = 0; i < fs->get_buffer_in_mem_pos(); i++) {sum_arg += fs_list[i].get_distance();}}
      if(!(sum_arg_before == sum_arg) && (debug_cnt_elements_written_to_file == 0)) {
	fprintf(stderr, 
		"!!\t In the merging, expected the sum of distances before (%f = %.3E) to match the sum of distances after (%f = %.3E)."
		"-\t In this operation we wrote %u=%.3E  elements to memory.\n"
		"-\t If this message was not understood, please contact the develoepr at [oekseth@gmail.com]\n"
		"This warning was printed at [%s]:%s:%d\n",
		sum_arg_before, sum_arg_before, sum_arg, sum_arg, debug_cnt_elements_written_to_file, (float)debug_cnt_elements_written_to_file, __FUNCTION__, __FILE__, __LINE__);
	assert(sum_arg_before == sum_arg); // To know that the list we got has not changed
      assert(sum_arg == data_sum_of_distances);
      
      //! Verfies that the sum of inserts corresponds to our expectations:
      const float actual_sum_inserted = sum_after - sum_before;
      log_builder::compare_floats(actual_sum_inserted, sum_arg, buffer_in_mem_pos, __LINE__, __FILE__, __FUNCTION__);

      //! Verfiies that the sum before- and after merging is correct:
      const float inserted_sum = sum_before_merging + data_sum_of_distances;
      log_builder::compare_floats(actual_sum_inserted, inserted_sum, buffer_in_mem_pos, __LINE__, __FILE__, __FUNCTION__);
      log_builder::compare_floats(expected_sum_of_distances_after_merge, inserted_sum, buffer_in_mem_pos, __LINE__, __FILE__, __FUNCTION__);
      }
#endif
    }
  //! Inserts the protein inside the buffer and, if the buffersize is above the limit, writes the data to the file
  void enlarge() {
    if(buffer_in_mem_pos >= buffer_in_mem_size) { // The next element do not fit into the buffer
      if(buffer_in_mem_pos >= MAX_BUFFER_IN_MEM_SIZE) {
	dump_buffer(); // Dumps the data
      } else { // < MAX_BUFFER_IN_MEM_SIZE: increases the size reserved in memory
	if(buffer_in_mem_size != 0) {
	  loint length_new = buffer_in_mem_size * 2;
	  if(buffer_in_mem_size < 100) length_new = buffer_in_mem_size + 10; 
	  assert(length_new > buffer_in_mem_pos);
	  T::realloc_list(list, buffer_in_mem_size, length_new,  buffer_in_mem_pos);
	  buffer_in_mem_size = length_new;
	} else { // Buffer is empty
	  if(list == NULL) {// reserves and initializes empty buffer
	    assert(buffer_in_mem_size == 0);
	    assert(buffer_in_mem_pos == 0);
	    buffer_in_mem_size = 2;
	    //	    buffer_in_mem_size = 1000; // Intiatix with a value.
	    buffer_in_mem_pos = 0;
	    // Reserves memory and initiates:
	    list = T::init_list(buffer_in_mem_size, buffer_in_mem_pos);
	  } else fprintf(stderr, "!!\tNot able to reserve memory of size(%llu) at line %d in file %s\n", buffer_in_mem_size, __LINE__, __FILE__);
	}
      }
    }
  }


 public:
  //! Inserts an element into the list:
  void insert(uint index_out, float sim_score, overlap_t overlap_in, overlap_t overlap_out);

  //! Increments the data using the last position
  void increment_data(float sim_score, overlap_t overlap_in, overlap_t overlap_out);

  //! @return true if buffer is filled with data
  bool buffer_has_data_set() {return (list != NULL);}

  //! @return true if the class has relations stored.
  bool has_data() {
    if ((file_cnt_structs > 0 || buffer_in_mem_pos) > 0) {
      return true;
    } else return false;
  }
  //! @return true if key-index is found in the data-set:
  bool has_data(uint search_index, T *buffer, uint start, uint end) {
    assert(end <= buffer_in_mem_size);
    for(uint k = start; k< end; k++) {
      if(buffer[k].is_equal(search_index)) {
	return true;
      }
    }
    return false;
  }
  //! @return true if key-index is found in the data-set:
  bool has_data(uint protein_out, index_t block) {
    if(block.get_start_pos() >= file_cnt_structs || (file_cnt_structs == 0)) {// its in the buffer
      const uint start =  block.get_start_pos() - file_cnt_structs;
      const uint end   =  block.get_end_pos()   - file_cnt_structs;
      return T::has_data(protein_out, list, start, end);
    } else { 	// Opens the file, reads the data:
      const uint start =  block.get_start_pos();
      const uint length = block.get_length();
      T *buff = T::get_file_block(file_name, start, length);
      const bool ret_val = T::has_data(protein_out, buff, 0, length);//start, start+length);
      T::free_file_block(buff, length);
      return ret_val;
    }
  }


  /**
     @brief Sums the memory consumption for the object.
     @remarks For anaylsing the memory fingerprint.
     @return the total length of data in Bytes.
  **/
  loint getTotal_memoryConsumption_for_object() {
    loint sum = sizeof(relations_list<T>);
    sum += buffer_in_mem_size * sizeof(T);
    return sum;
  }

  /**
     @brief Sums the memory consumption for the object, disregarding the non-used containers.
     @remarks For anaylsing the memory fingerprint.
     @return the total length of data in Bytes.
  **/
  loint getTotal_used_memoryConsumption_for_object() {
    loint sum = sizeof(relations_list<T>);
    sum += buffer_in_mem_pos * sizeof(T);
    return sum;
  }

  /**
     @brief Gets the object size.
     @remarks Uses the deault size-init-values of the objects of type T
     @return the size of an empty object.
  **/
  static loint get_empty_object_size() {
    loint total = sizeof(relations_list<T>); // Do not allocate any objects of type 'index'
    return total;
  }

  // ------------------------------------------------------------------------------

  // -----------------------de-Alloctor-methods-------------------------

  //! Frees a list of this class.
  void free_file_block(T *&block, uint size) {
    tbb::tbb_allocator<T>().deallocate((T*)block, size);
  }

  //! De-allocates the object and if specified de-allocates the object.
  void free_memory(bool remove_file) {
    if(list) {T::close(list, buffer_in_mem_pos, buffer_in_mem_size);}//delete  [] list; list = NULL;}
    if(file_name && remove_file) {
      remove(file_name); // removes the file.
      delete [] file_name; file_name = NULL;
    } else if(file_name) {delete [] file_name, file_name = NULL;}
  }
  //! De-allocates the object and if specified de-allocates the object.
  static void close(relations_list *&obj, const bool remove_file) {
    if(obj) {obj->free_memory(remove_file); delete obj; obj = NULL;}
  }

  //! Initiates a type T list belogning to this class:
  void init_T_data() {
    assert(list);
    for(uint i = 0; i<buffer_in_mem_size; i++) list[i] = T();
  }
  /**
     @brief Init a list of these.
     @param <buffer_size> The number of elements to initialize
     @param <buffer> The data to update.
  **/
  void init(loint buffer_size, T *&buffer) {
    assert(buffer);
    for(uint i = 0; i<buffer_size; i++) buffer[i] = T();
  }
  //! Init a list of these.
  void init(uint start, uint end, T *&buffer) {
    for(;start<end; start++) buffer[start] = T();
  }
 

  // -----------------------Alloctor-methods-------------------------

  //! Reserves a list of these.
  T *reserve_list(loint buffer_size) {
    return T::init(buffer_size); //(T*)malloc(sizeof(T)*buffer_size);
  }

  //! Init a list of these.
  T *init_list(uint buffer_size) {
    return T::init(buffer_size); //buffer;
  }

  //! Deletes a list of these
  void delete_list(T *&buffer) {T::close(buffer, buffer_in_mem_pos, buffer_in_mem_size);} 

  //! @return the object of type T wihtout setting the file name, implying that it by def may not write object data to a binary file
  static relations_list<T> *init(bool USE_BEST_BLAST_PAIR_SCORE) {
    relations_list<T> *obj = NULL; // new relations_list<T>(USE_BEST_BLAST_PAIR_SCORE);
    try {obj = new relations_list<T>(USE_BEST_BLAST_PAIR_SCORE);} 
    catch (std::exception& ba) {
      fprintf(stderr, "!!\t Tried reserving a \"relations_list<%s>\" class set of length=%u, [%s]:%s:%d\n",  T::get_class_name(),1, __FUNCTION__, __FILE__, __LINE__);
      if(!log_builder::catch_memory_exeception(1, __FUNCTION__, __FILE__, __LINE__, true)) {
	fprintf(stderr, "!!\t An interesting error was discovered: %s."
		"The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		"Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
      }
    }
    if(obj) return obj;
    else {
      fprintf(stderr, "!!\tUnable to reserve memory for the relations_list<T> object (ie. Probably implies that your memory chip is too small). If questions, please contact the developer at oekseth@gmail.com, giving the following information: This message was genereated at line %d in file %s, found in method %s\n", __LINE__, __FILE__, __FUNCTION__);
      exit(2); // No point in continuing the work flox.
    }
  }

  //! @return the object of type T wihtout setting the file name, implying that it by def may not write object data to a binary file
  static relations_list<T> *init(uint taxon_in, uint taxon_out, bool USE_BEST_BLAST_PAIR_SCORE) {
    relations_list<T> *obj = new relations_list(taxon_in, taxon_out, USE_BEST_BLAST_PAIR_SCORE);
    if(obj) return obj;
    else {
      fprintf(stderr, "!!\tUnable to reserve memory for the relations_list<T> object (ie. Probably implies that your memory chip is too small). If questions, please contact the developer at oekseth@gmail.com, giving the following information: This message was genereated at line %d in file %s, found in method %s\n", __LINE__, __FILE__, __FUNCTION__);
      exit(2); // No point in continuing the work flox.
    }
  }

  /**
     @brief Builds a relation_list object with possibility of using temporary files as a storage.
     @param <taxon_in> The taxon id for the inner relation.
     @param <taxon_out> The taxon id for the outer relation.
     @param <block_number> The block id for this container.
     @param <max_memory_to_reserve> The maximal amount of memory to reserve before object starts dumping data to a binary file.
     @param <FILE_BINARY_LOCATION> The location (folder) to store the temporary binary files in: Therby if folder is not given write- and read priveligis for hte amount of data, a failiure could occure.
     @param <_USE_BEST_BLAST_PAIR_SCORE> If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
     @return an initialized relation_list object.
     @author Ole Kristian Ekseth (oekseth)
  **/
  static relations_list<T> *init(uint taxon_in, uint taxon_out,
				 int block_number, loint max_memory_to_reserve,
				 char *FILE_BINARY_LOCATION, bool _USE_BEST_BLAST_PAIR_SCORE) {
    relations_list<T> *obj =  new relations_list(taxon_in, taxon_out, block_number,
			      max_memory_to_reserve, FILE_BINARY_LOCATION, _USE_BEST_BLAST_PAIR_SCORE);
    if(obj) return obj;
    else {
      fprintf(stderr, "!!\tUnable to reserve memory for the relations_list<T> object (ie. Probably implies that your memory chip is too small). If questions, please contact the developer at oekseth@gmail.com, giving the following information: This message was genereated at line %d in file %s, found in method %s\n", __LINE__, __FILE__, __FUNCTION__);
      exit(2); // No point in continuing the work flox.
    }
  }
  //! Constructing the object wihtout setting the file name, implying that it by def may not write object data to a binary file
 relations_list(bool _USE_BEST_BLAST_PAIR_SCORE) :
  list(NULL), file_cnt_structs(0), buffer_in_mem_pos(0), buffer_in_mem_size(0),
    file_name(NULL), block_number(0), 
    my_id(-1), taxon_in(-1), taxon_out(-1), USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE)
  {
    MAX_BUFFER_IN_MEM_SIZE = UINT_MAX; // Do not know the location to store, therefore should not write data to binary files. //1000*sizeof(T);
  };
   
  //! Constructing the object wihtout setting the file name, implying that it by def may not write object data to a binary file
 relations_list(uint in, uint out, bool _USE_BEST_BLAST_PAIR_SCORE) :
  list(NULL), file_cnt_structs(0), buffer_in_mem_pos(0), buffer_in_mem_size(0),
    file_name(NULL), block_number(0), 
    my_id(-1), taxon_in(in), taxon_out(out), USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE)
  {
    MAX_BUFFER_IN_MEM_SIZE = UINT_MAX; // Do not know the location to store, therefore should not write data to binary files.//1000*sizeof(T);
  };

  /**
     @brief Constructs a relation_list object with possibility of using temporary files as a storage.
     @param <in> The taxon id for the inner relation.
     @param <out> The taxon id for the outer relation.
     @param <_block_number> The block id for this container.
     @param <mem_reserved> The maximal amount of memory to reserve before object starts dumping data to a binary file.
     @param <FILE_BINARY_LOCATION> The location (folder) to store the temporary binary files in: Therby if folder is not given write- and read priveligis for hte amount of data, a failiure could occure.
     @param <_USE_BEST_BLAST_PAIR_SCORE> If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
     @author Ole Kristian Ekseth (oekseth)
  **/
 relations_list(uint in, uint out, 
		int _block_number, loint mem_reserved, char *FILE_BINARY_LOCATION, bool _USE_BEST_BLAST_PAIR_SCORE) :
  list(NULL), file_cnt_structs(0), buffer_in_mem_pos(0), buffer_in_mem_size(0),
    file_name(NULL), block_number(_block_number), 
    my_id(-1), taxon_in(in), taxon_out(out), USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE)
  {
    MAX_BUFFER_IN_MEM_SIZE = mem_reserved;
    file_name = set_file_name(in, out, FILE_BINARY_LOCATION);
  };

  //! Tests writing of and reading for the object of interest
  void assert_buffer_write_and_read() {
#ifdef assert_code  
    uint insert_total = 100;//     uint insert_total = 10;
    // Build an new buffer: To avoid file-writing, insertes diretly not using wrapper-function:
    buffer_in_mem_pos = 0; buffer_in_mem_size = insert_total;
    init_buffer(insert_total, list); 
    for(uint index_in =0; index_in<insert_total; index_in++) {      
      // Asserts that the starting postion is correct:
      uint index_out = 0;   if(index_in == 0) index_out = 1;
      //      insert_new_rel(index_in, index_out, (float)index_in, index_in); 
      list[buffer_in_mem_pos].set_data(index_out, (float)index_in, index_in);
      buffer_in_mem_pos++;
    }
    T *data = get_file_content();
    for(uint index_in =0; index_in<insert_total; index_in++) {
      assert(list[index_in].is_equal(data[index_in]));
    }
    free_memory(true);  // free(data); data = NULL;
    // TODO: Some more tests could be considered.
#endif
  }
  /**
     @brief The main test function for this class.
     @remarks This method includes:
     - Formalized tests ..
     - Valgrind used verifying memory usage.
     - Examples of how this class may be used.
     @author Ole Kristian Ekseth (oekseth)
     @date 04.01.2012 by oekseth.
  **/
  static bool assert_class(const bool print_info) {
    const static char *class_name = "relations_list";
    if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
    uint in=0, out=0, block_number=0, mem_reserved=1;
    char *FILE_BINARY_LOCATION = ".";
    const bool USE_BEST_BLAST_PAIR_SCORE = false;
    relations_list<p_rel> *test_1 = relations_list<p_rel>::init(in, out,  block_number, mem_reserved, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);


    relations_list<p_rel>::close(test_1, true);
#endif
    if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
    return true;
  }
};

#endif
