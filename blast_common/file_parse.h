#ifndef file_parse_h
#define file_parse_h
/**
   @file
   @brief Contains the blast file representation for a bunch of taxon-pairs.
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
v * You should have received a copy of the GNU General Public License
 * along with orthAgogue
. If not, see <http://www.gnu.org/licenses/>.
 */
#include "parse.h"
#include "index_list.h"
#include "norm_t.h"
#include "taxa.h"
#include "basket_parse.h"
#include "relations_list.h"
#include "prot_list.h"
#include "mcl_format.h"
/**
   @class file_parse.
   @brief Contains the blast file representation for a bunch of taxon-pairs.
   @tparam <T> A rel class for the relations_list object.
   @ingroup blastfile_container
   @remark: The data is filles data about a taxa set for a combination of two taxa,
   representing a sparse sub-matrix, as shown in the illustrations provided in
   the documentation. The class may be used directly, or through the inferface found
   in 'list_file_parse.h'. Included are complex tests spanning several aspects.
   In view of time constraints, not all apsects are tested.
   @author Ole Kristian Ekseth (oekseth).
   @date 12.01.2011 by oekseth (initial).
   @date 08.09.2011 by oekseth (asserts).
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/ 
template<class T> class file_parse { // Relations in the file
   private:
  // public:
  index_list_t *listIndex;  // The index refering the protein id in taxa(1)-->taxa(2)
  relations_list<T> *listRelations;
  int taxon_in, taxon_out;   int my_id;
  bool USE_BEST_BLAST_PAIR_SCORE; // If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
 public:
  /**
     @brief Inserts a complete object
     @remarks
     # Primary usage is for MPI, ie, when macro variable "USE_MPI" is activated.
     # Uses a safe procedure, ie, copies the data after first allocating memory.
  **/
  void insert_object_using_1d_lists(uint taxon_in, uint taxon_out, char *FILE_BINARY_LOCATION, const bool USE_BEST_BLAST_PAIR_SCORE, T *rel_list, const uint rel_list_size, index_t *index_list, const uint index_list_size, uint MAX_BUFFER_SIZE) {
    assert(!listIndex);
    listIndex = index_list::init(index_list_size);
    listIndex->insert_new_list(index_list, index_list_size);
    assert(!listRelations);
    listRelations = relations_list<T>::init(taxon_in, taxon_out, /*block_id=*/0, MAX_BUFFER_SIZE, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);
    listRelations->insert_new_list(rel_list, rel_list_size);
  }
  //! @return the file content only if availability of enough memory
  T *get_file_content() {
    if(get_number_of_objects_in_file() == 0 && listRelations) return listRelations->get_buffer();
    else if(listRelations)    return listRelations->get_file_content();
    else return NULL;
  }
  //! @return the content of data written to the file
  T *get_file_content(T *&buff) {
    if(listRelations) return listRelations->read_file_content(buff);
    else return NULL;
  }

  //! @return the index_t list
  index_t *get_index_list() {
    if(listIndex) return listIndex->get_list();
    else return NULL;
  }
  //! Prints info about the object
  void print_relation_object_info(FILE *f) {
    if(listRelations) listRelations->print_object_info(f);
  }

  //! Returns the distance, given the pair identified by the arguments set.
  uint get_distance(uint protein_in, uint protein_out) {
    assert(listIndex);
    assert(listRelations);
    fprintf(stderr, "!!\tFunction %s in file %s not yet implemented. Belongs to a future scope. Contact oekseth@gmail.com.\n", __FUNCTION__, __FILE__);
    return 0;
  }
  //! Sets the USE_BEST_BLAST_PAIR_SCORE value.
  void set_use_best_blast_pair_score(bool _USE){USE_BEST_BLAST_PAIR_SCORE = _USE;}

  //! Returns true if only using the best blast pair score is set.
  bool get_use_best_blast_pair_score(){return USE_BEST_BLAST_PAIR_SCORE;}

  //! Sets data to empty value:
  void set_to_empty() {
    if(listIndex) {index_list::close(listIndex);}
    if(listRelations) {relations_list<T>::close(listRelations, true);}
    taxon_in = -1, taxon_out = -1;
    my_id = -1;
  }
  //! Prints the list using class index for the purpose.
  void print_list(){
    if(listIndex) listIndex->print_list();
  }
  //! Prints the buffer
  void print_buffer() { if(listRelations) listRelations->print_buffer();}
  //! Prints the content of the file
  void print_file(){  if(listRelations) listRelations->print_file();}
  //! Compares this structure with the input structure; print a warning if they not equal
  int debug_compare_structs(file_parse fileParse, int i, int k, char **arrKey_in, char **arrKey_out);
  //! @return the 'file_cnt_structs' variable.
  uint getLengthOfFile() {
    if(listRelations) return (uint)listRelations->get_file_length();
    else return 0;
  }

  
  //! @return the file_cnt_structs variable, defining the count of elements in the file
  uint get_number_of_objects_in_file() {return getLengthOfFile();}

  /**
     @brief Sums the memory consumption for the object.
     @remarks For anaylsing the memory fingerprint.
     @return the total length of data in Bytes.
  **/
  loint getTotal_memoryConsumption_for_object() {
    loint sum =       sizeof(file_parse<T>);
    if(listIndex)     sum += listIndex->getTotal_memoryConsumption_for_object();
    if(listRelations) sum += listRelations->getTotal_memoryConsumption_for_object();
    return sum;
  }

  //! @return the number of vectors for the proteins in this range
  uint getTotalRelations(uint protein_start, uint protein_end) {
    if(listIndex) return listIndex->get_number_of_elements(protein_start, protein_end);
    else return 0;
  }

  /**
     @brief Sums the memory consumption for the object, disregarding the non-used containers.
     @remarks For anaylsing the memory fingerprint.
     @return the total length of data in Bytes.
  **/
  loint getTotal_used_memoryConsumption_for_object() {
    loint sum =       sizeof(file_parse<T>);
    if(listIndex)     sum += listIndex->getTotal_memoryConsumption_for_object();
    if(listRelations) sum += listRelations->getTotal_used_memoryConsumption_for_object();
    return sum;
  }
  /**
     @return The imaginary next position in file.
     @remarks Used to avoid sending data received from others.
  **/
  uint get_index_T_end_of_this_object() {
    if(listRelations) return listRelations->get_imaginary_next_position_in_file();
    else return 0;
  }
  //! @return the sum of 'file_cnt_structs+buffer_in_mem_pos' 
  uint getTotalLengthOfData() {
    if(listRelations) return (uint)listRelations->getTotalLengthOfData();
    else return 0;
  }  

  //! @return the sum of the pair distances found in the T struct
  float get_sum_of_pair_distances_for_myrank_in_memory() {
    if(listRelations) return listRelations->get_sum_of_pair_distances_in_memory();
    else return 0.0;
  }  

  //! @return the sum of references in the listIndex
  uint getTotalLengthOfData_for_index() {
    if(listIndex) return listIndex->getTotalLengthOfData();
    else return 0;
  }  
  //! Returns the total length of the index_t list of objects:
  uint get_index_length() {
    if(listIndex) return listIndex->get_index_length();
    else return 0;
  }
  //! @return the buffer assuming its eigher in memory, 
  basket_parse *getBuffer() {
    if(listRelations)  return listRelations->get_buffer_basket();
    else return NULL;
  }

  //! @return the last index inseted in the T type buffer
  loint get_buffer_in_mem_pos() {
    if(listRelations) return listRelations->get_buffer_in_mem_pos();
    else return 0;
  }
  /**
     @brief Gets the partial buffer, given the local protein index
     @remarks Data is not allocated.
     @return The buffer for the proteins of type T
  **/
  T *getBuffer(uint protein_in)  {
    if(listRelations && listIndex) {
      return listRelations->getBuffer(listIndex->get_start_pos(protein_in), listIndex->get_length(protein_in));
    } else return NULL;
  }
  //! @return the start pos in the data for the given index; UINT_MAX if not found.
  uint get_index_start(uint protein_in) {
    if(listIndex) return listIndex->get_index_start(protein_in);
    else return UINT_MAX;
  }

  //! @return the start pos in the data for the given index:
  uint getProteinStartPos(uint protein_id) {return get_index_start(protein_id);}

  //! @return the start pos for the given index; UINT_MAX if not found.
  uint get_absolute_start_pos(uint protein_in) {
    if(listIndex) return listIndex->get_absolute_start_pos(protein_in);
    else return UINT_MAX;
  }
  //! @return the length of the proteins in the collection.
  int getLengthProteins() {
    if(listIndex) return listIndex->get_index_reserved();
    else return 0;
  }
  //! @return the number of protein relations for the index
  uint get_index_length(uint protein_in) {
    if(listIndex)     return listIndex->get_index_length(protein_in);
    else return 0;
  }

  //! @return the number of protein relations reserved for the index
  uint get_index_reserved() {
    if(listIndex)     return listIndex->get_index_reserved();
    else return 0;
  }
  //! @return the number of protein relations for the index
  uint getProteinLength(uint protein_id) {return get_index_length(protein_id);}

  //! @return the number of protein relations for the index
  uint get_total_pair_cnt() {
    if(listRelations) return listRelations->getTotalLengthOfData();
    else {
      fprintf(stderr, "!!\tNo data reserved for [%d][%d] in file_parse at line %d\n", taxon_in, taxon_out, __LINE__);
      return 0;
    }
  }

  //! Writes specific info about this object:
  void log_write_info(FILE *f, taxa_t *listTaxa, const bool DEBUG_print_pairs_in_file_parse_log_file) {
    fprintf(f, "-\ttaxon[%d][%d]\t", taxon_in, taxon_out);
    if(listTaxa) fprintf(f, "-\ttaxon[%s][%s]\t", listTaxa[taxon_in].getName(), listTaxa[taxon_out].getName());
    if(listIndex) fprintf(f, "listIndex{used=%d}\t", listIndex->get_index_used());
    //    if(listIndex) fprintf(f, "listIndex{used=%d, reserved=%d}\t", listIndex->get_index_used(), listIndex->get_index_reserved());
    if(listRelations) {
      fprintf(f, "listRelations{used=%d, in_file=%d}\t",
	      (int)listRelations->get_buffer_in_mem_pos(), (int)listRelations->get_file_cnt_structs());
/*       fprintf(f, "listRelations{used=%d, reserved=%d, in_file=%d}\t", */
/* 	      (int)listRelations->get_buffer_in_mem_pos(), (int)listRelations->get_buffer_in_mem_size(), */
/* 	      (int)listRelations->get_file_cnt_structs()); */
      fprintf(f, "\n"); 
      if(DEBUG_print_pairs_in_file_parse_log_file) {
	// The below output will be huge: therefore discard it in real-life situations
	listRelations->print_buffer_as_abc_files(f, listIndex, NULL, taxon_in, taxon_out); 
	//	listRelations->print_buffer_as_abc_files(f, listIndex, listTaxa, taxon_in, taxon_out); 
      }
    } else fprintf(f, "\n");
  }

  //! Writes the buffer to stdout, using integres as representation.
  void print_buffer_as_integer_file(uint taxon_in, uint taxon_out) {
    if(listRelations) 	listRelations->print_buffer_as_abc_files(stdout, listIndex, NULL, taxon_in, taxon_out); 
  }

  //! Writes the buffer to the file specified, using integres as representation.
  void print_buffer_as_integer_file(uint taxon_in, uint taxon_out, FILE *f) {
    if(listRelations) 	listRelations->print_buffer_as_abc_files(f, listIndex, NULL, taxon_in, taxon_out); 
  }
 
  //! The actual number of elements of the buffer
  mem_loc getLengthBuffer() {
    if(listRelations) return listRelations->getTotalLengthOfData();
    else return 0;
  }
  //! Sets the my_id variable
  void set_my_id(int id){my_id = id;}
  //! @return true if the right proteins are equal
  bool insert_protein_index(mem_loc index_in, mem_loc index_out) {
    assert(listRelations);
    const loint pos_in_rel_struct = listRelations->get_imaginary_next_position_in_file();
    return listIndex->insert(index_in, index_out, pos_in_rel_struct);
  }

  /**
     The 'main method': Inserts the data into the structure ('this')
     @return true if protein was inserted as a new protein, and not updating an exsisting one.
     @remarks Returnvalue to be used in conjunction with the method named 'getTotalLengthOfData()' found in class list_file_parse.
  **/
  bool insert_new_rel(uint index_in, uint index_out, float sim_score, overlap_t overlap_in, overlap_t overlap_out) {
    assert(listRelations);
    assert(listIndex);
    bool return_value = true;    
    const bool equal_out = insert_protein_index(index_in, index_out);
    if(!are_equal(index_in, index_out)) { // not equal
      // Tests if the previous right protein and this right protein are not equal  :
      if(equal_out == false) { 
	listRelations->insert(index_out,sim_score, overlap_in, overlap_out);
	listIndex->increment_length(index_in);
      } else {
	listRelations->increment_data(sim_score, overlap_in, overlap_out);
	return_value = false; // They were not inserted;
      }
    } else {// else are equal: ignores.
      return_value = false; // They were not inserted;
    }
    listIndex->set_index_in_prev(index_in);  
    listIndex->set_index_out_prev(index_out);
    return return_value;
  }


  /**
     @brief Produces the inparalogs, inserting them in the given container
     @return the number of inparalogs inserted
  **/
  uint produceInparalogs(taxa *listTaxa, bool MODE_PAIRWISE_OUTPUT_ABC, bool MODE_PAIRWISE_OUTPUT_MCL, mcl_format_t *container, char *protein_name_in, uint taxon, uint world_index, uint local_index/*, uint offset*/, float avg_div_factor) {
    if(listRelations && listIndex) {
      const int index_length = listIndex->get_index_length(local_index);      
      if(index_length > 0) {
	//! Formalizes the number of inparalogs inserted in this operation:
	const uint debug_cnt_inparalogs_inserted = (uint)index_length; 

	//! The procedure of building, using specific settings:
	const uint index_start_position = listIndex->get_start_pos(local_index);	
	//    if (index[local_index].get_length() > 0) {
	//	rel_t *buffer = listRelations->getBuffer(local_index/*, OPT_GET_FILE_LOCALLY, NULL, NULL, offset, offset*/);
	if(MODE_PAIRWISE_OUTPUT_ABC && MODE_PAIRWISE_OUTPUT_MCL) { // both are pairwise
	  for(uint i = 0; i<(uint)index_length; i++) {
	  //	  for(uint i = 0; i<index[local_index].get_length(); i++) {
	    T element = listRelations->get_element_at_local_index(index_start_position + i);
	    const uint index_real = listTaxa[taxon].getWorldIndex(element.ind_out);
	    char *protein_name = listTaxa[taxon].getCompleteProteinName(element.ind_out);
	    container->insert_inpa(world_index, index_real, protein_name_in, protein_name, element.distance, avg_div_factor);
	    delete [] protein_name;
	  }
	} else if(!MODE_PAIRWISE_OUTPUT_ABC && !MODE_PAIRWISE_OUTPUT_MCL){ // non are pairwise
	  container->set_header_inpa(world_index, protein_name_in);
	  //	  for(uint i = 0; i<index[local_index].get_length(); i++) {
	  for(uint i = 0; i<(uint)index_length; i++) {
	    T element = listRelations->get_element_at_local_index(index_start_position + i);
	    const uint index_real = listTaxa[taxon].getWorldIndex(element.ind_out);
	    char *protein_name = listTaxa[taxon].getCompleteProteinName(element.ind_out);
	    container->insert_inpa(index_real, protein_name, element.distance, avg_div_factor);
	    delete [] protein_name;
	  }
	} else {
	  container->set_header_inpa(world_index, protein_name_in);
	  for(uint i = 0; i<(uint)index_length; i++) {
	    T element = listRelations->get_element_at_local_index(index_start_position + i);
	    //	  for(uint i = 0; i<index[local_index].get_length(); i++) {
	    const uint index_real = listTaxa[taxon].getWorldIndex(element.ind_out);
	    char *protein_name = listTaxa[taxon].getCompleteProteinName(element.ind_out);
	    container->insert_inpa(world_index, index_real, protein_name_in, protein_name, element.distance, avg_div_factor);
	    delete [] protein_name;
	  }
	}
	return debug_cnt_inparalogs_inserted;
      }; 
    }
    return 0;
  }

  //! Tests if the previous value is equal
  bool is_equal(int in_prev, int out_prev) {
    assert(listIndex);
    return listIndex->index_are_equal(in_prev, out_prev);  
  }
  /**  @return true if the indexes corresponds to the value in the object   **/
  bool are_equal(mem_loc index_in, mem_loc index_out) {
    return ((index_in == index_out) && (taxon_in == taxon_out));
  }
  //! Compare the two classes:
  bool is_equal(file_parse data, const bool print_info) {
    bool is_equal = true;
    assert(listIndex);
    assert(listRelations);
    is_equal = listIndex->is_equal(data.listIndex, print_info);
    is_equal = listRelations->is_equal(data.listRelations, print_info);
    if(print_info && !is_equal) printf("\n");
    return is_equal;
  }
  //! @return true if class is empty
  bool is_empty() {//bool print_info) {
    if(!listRelations && !listIndex) return true;
    else return false;
  }
  //! @brief true if there is an overlap. @remarks Used to remove overlaps between different reads.
  bool gives_overlap(int in,int out, mem_loc protein_in, mem_loc protein_out) {
    if(taxon_in == in && taxon_out == out) { // Possible merge of the first element of this towards teh last element of the global list
      return is_equal(protein_in, protein_out);
    }
    return false;
  }

  //! Decreases the size of a list element given its index.
  void decreaseLengthOfListElement(int protein_index) {
    if(has_data(protein_index)) {
      listIndex->decrease_length(protein_index);
    }
  }

  //! @return true if buffer is filled with data
  bool buffer_has_data_set() {
    if(listRelations) {
      const bool has_data_in_buffer = listRelations->buffer_has_data_set();
      if(has_data_in_buffer) {
	assert(0 != listRelations->get_buffer_in_mem_size());
      }
      return has_data_in_buffer;
    }
    else return 0;
  }
  //! @return true if the class has relations stored.
  bool has_data() {
    if(listRelations) return listRelations->has_data();
    else return false;
  }
  //! @return true if the class has relations stored.
  bool has_data(uint protein_in) {
    if(listIndex)    return listIndex->has_data(protein_in);
    else return false;
  }
  //! @return true if data set for a given protein.
  bool has_data_slow(uint protein_in, uint protein_out) {
    if(has_data(protein_in)) {
      index_t block; listIndex->get_block(protein_in, block);
      assert(listRelations);
      return listRelations->has_data(protein_out, block);
    } else return false; 
  }

  /**
     @brief Grabs the data, and intialises the values to the argument into zero (empy) values.
     @param <data> The file_parse object to grab data from.
     @remarks Deletes own data, if any, before grapbing the objects data.
  **/
  void steal_data(file_parse *data) {
    //    file_parse data = data_[0];
    assert(data);
    //    if(true) {
    index_list::close(listIndex);
    relations_list<T>::close(listRelations, false);
    listIndex = data->listIndex;
    listRelations = data->listRelations;
    USE_BEST_BLAST_PAIR_SCORE = data->get_use_best_blast_pair_score();
    data->listIndex = NULL;
    data->listRelations = NULL;
    taxon_in = data->taxon_in; 
    taxon_out = data->taxon_out; 
  }

  /**
     @brief Grabs the data out of 'this' object, and intialises the values to the argument into zero (empy) values.
     @param <data> The file_parse object to put into
  **/
  void steal_data_inverse(file_parse *data) {
    //    file_parse data = data_[0];
    assert(data);
    if(!data->has_data()) {
      data->listIndex = listIndex; listIndex = NULL;
      data->listRelations = listRelations; listRelations = NULL;
      data->taxon_in = taxon_in;
      data->taxon_out = taxon_out;
      data->set_use_best_blast_pair_score(USE_BEST_BLAST_PAIR_SCORE);
    } else fprintf(stderr, "!!\tWrong usage of this method, as the input should not been set. Please refer to the manual or contact the developer at oekset@gmail.com. (This error is found in %s at line %d in method %s\n", __FILE__, __LINE__, __FUNCTION__);
  }


  //! Updates the size of the index
  void enlarge_index(int index_in) {
    assert(listIndex);
    listIndex->enlarge(index_in);
  }
  //! Writes the buffer, in both 'this' and the argument given, to the file set in 'this'
  void write_buffer_to_file(file_parse parse_b);

  //! Dumps the buffer to a file
  void dump_buffer_to_file() {
    assert(listRelations);
    listRelations->dump_buffer();
  }

  //! If data resides in files, get them into the local buffer
  void getFileIntoBuffer() {
    dumpBufferOrFile(UINT_MAX);
  }
  //! Dumps the file to a buffer if the buffer is less than the memory limit set
  void dumpBufferOrFile(uint max_buffer_size) {
    if(listRelations && listRelations->dumpBufferOrFile(max_buffer_size)) {
      index_list::close(listIndex);
    }
  }
  /**
     @brief Merges an input object with this, only updating the object 'this' internal object.
     @remarks Merges the argument into this object:
     -# Puts the reference table (of type index_list_t) into the one this object uses, ie, merging of them.
     -# Merges the containers having the pairs. @attention The argument-data is placed at the end of 'this' container.
     @author Ole Kristian Ekseth (oekseth).
  **/
  void merge_buffers(file_parse *parse_b_, taxa &taxa_obj) {
    file_parse parse_b = *parse_b_;
    assert(listIndex);
    assert(listRelations);
    if(parse_b.has_data()) {
    //    if(parse_b.getTotalLengthOfData()) { //buffer_has_data_set()) { // argument has data in memory:
    //    if(parse_b.buffer_has_data_set()) { // argument has data in memory:
#ifndef NDEBUG
      const uint tot_elements = getTotalLengthOfData() + parse_b.getTotalLengthOfData();
#endif
      // The starting position in the buffer to where the index shall start from:
      const loint offset_to_add_data_to =listRelations->get_imaginary_next_position_in_file();
      // 1. Merges the argument buffer with 'this' index (found in the input class)
      listIndex->merge_buffers(parse_b.listIndex, parse_b.getLengthProteins(), offset_to_add_data_to, taxa_obj);
      // 2. Merges the buffers containing the 'pairs' them self.
      listRelations->merge_buffers(parse_b.listRelations);
      const loint buffer_in_mem_pos = listRelations->get_buffer_in_mem_pos();
      if(buffer_in_mem_pos>0) { // If the buffer contains data:
	const loint rel_out = listRelations->get_ind_out(buffer_in_mem_pos-1);
	listIndex->set_index_out_prev(rel_out);
      }
#ifndef NDEBUG
      assert(tot_elements == getTotalLengthOfData()); // The 'root' should now consist of the sum of elements.
#endif
    }
  }
  //! Extending, and intializing, the buffer:
  void extend_buffer(mem_loc new_size);
  //!  Copies the buffer, given as input argument
  void copy_buffer(mem_loc buffer_in_mem_pos, T *T_buffer, mem_loc T_size);
  //! Inserts the protein inside the buffer and, if the buffersize is above the limit, writes the data to the file
  void enlarge_buffer();
  /**
     @brief Intitiates a squared list to be used representing a complete taxa.
     @param <taxon_length> The number of taxa to be processed.
     @param <hashProtein>  Used setting the expected length for each taxa-combination.
     @param <my_id> In order to identify each thread.
     @param <max_buff_size> If this limit is passed, data is written to intermediary file(s).
     @param <FILE_BINARY_LOCATION> The location of where to store the intermediary file(s).
     @param <_USE_BEST_BLAST_PAIR_SCORE> If set to true, instead of summing the input values for the same protein pair, only uses the best pair.
   **/
  static file_parse<T> **init(uint taxon_length, prot_list_t *hashProtein, uint my_id, long int max_buff_size, char *FILE_BINARY_LOCATION, bool _USE_BEST_BLAST_PAIR_SCORE) {
    assert(taxon_length > 0);
    if(taxon_length > 0) {
      file_parse<T> **list = new file_parse<T>*[taxon_length];
      if(list) {
	for(uint i = 0; i < taxon_length; i++) {
	  list[i] = new file_parse<T>[taxon_length]();
	  if(list[i]) {
	    for(uint out = 0; out <taxon_length; out++) {
	      list[i][out] = file_parse<T>(i, out, hashProtein[i].getLength(), my_id, max_buff_size, FILE_BINARY_LOCATION, _USE_BEST_BLAST_PAIR_SCORE); 
	      list[i][out].set_my_id(my_id);
	    }
	  } else {
	    fprintf(stderr, "!!\tUnable to reserve memory for the *file_parse-object of taxa-length %u (ie. Probably implies that your memory chip is too small). If questions, please contact the developer at oekseth@gmail.com, giving the following information: This message was genereated at line %d in file %s, found in method %s\n", taxon_length, __LINE__, __FILE__, __FUNCTION__);
	    exit(2); // No point in continuing the work flox.
	  }
	}
      } else {
	fprintf(stderr, "!!\tUnable to reserve memory for the **file_parse-object of taxa-length  %u (ie. Probably implies that your memory chip is too small). If questions, please contact the developer at oekseth@gmail.com, giving the following information: This message was genereated at line %d in file %s, found in method %s\n", taxon_length, __LINE__, __FILE__, __FUNCTION__);
	exit(2); // No point in continuing the work flox.
      }
      return list;
    } else fprintf(stderr, "!!\tTaxon length not set: Discard initializing the list in 'list_file_parse' class found at line %d in file %s. Please contact oekseth@gmail.com if this mesage is seen.\n", __LINE__, __FILE__);
    return NULL;
  }

  /**
     @brief Intitiates a squared list to be used representing a complete taxa.
     @remarks By dwefault initiated to empty (NULL) values.
     @param <taxon_length> The number of taxa to be processed.
   **/
  static file_parse<T> ***init(uint taxon_length) {
    assert(taxon_length > 0);
    if(taxon_length > 0) {
      file_parse<T> ***list = new file_parse<T>**[taxon_length];
      if(list) {
	for(uint i = 0; i < taxon_length; i++) {
	  list[i] = new file_parse<T>*[taxon_length]();
	  if(list[i]) {
	    for(uint out = 0; out <taxon_length; out++) list[i][out] = NULL;
	  } else {
	    fprintf(stderr, "!!\tUnable to reserve memory for the *file_parse-object of taxa-length %u (ie. Probably implies that your memory chip is too small). If questions, please contact the developer at oekseth@gmail.com, giving the following information: This message was genereated at line %d in file %s, found in method %s\n", taxon_length, __LINE__, __FILE__, __FUNCTION__);
	    exit(2); // No point in continuing the work flox.
	  }
	}
      } else {
	fprintf(stderr, "!!\tUnable to reserve memory for the **file_parse-object of taxa-length  %u (ie. Probably implies that your memory chip is too small). If questions, please contact the developer at oekseth@gmail.com, giving the following information: This message was genereated at line %d in file %s, found in method %s\n", taxon_length, __LINE__, __FILE__, __FUNCTION__);
	exit(2); // No point in continuing the work flox.
      }
      return list;
    } else fprintf(stderr, "!!\tTaxon length not set: Discard initializing the list in 'list_file_parse' class found at line %d in file %s. Please contact oekseth@gmail.com if this mesage is seen.\n", __LINE__, __FILE__);
    return NULL;
  }

  /**
     @brief Gets the object size.
     @remarks Uses the deault size-init-values of the objects of type index_list and relations_list to produce an extimated size.
     @return the size of an empty object.
  **/
  static loint get_empty_object_size() {
    loint total = sizeof(file_parse<T>);
    total += index_list::get_empty_object_size();
    total += relations_list<T>::get_empty_object_size();
    return total;
  }
  /**
     @brief Intitiates a squared list with minimal data.
     @param <taxon_length> The number of taxa to be processed.
     @param <_USE_BEST_BLAST_PAIR_SCORE> If set to true, instead of summing the input values for the same protein pair, only uses the best pair.
     @remarks Calling this method implies that (in most contexts) intermdiary files will not be used
   **/
  static file_parse<T> ***init(uint taxon_length, bool _USE_BEST_BLAST_PAIR_SCORE) {
    if(taxon_length > 0) {
      file_parse<T> ***list = new file_parse<T>**[taxon_length];
      for(uint i =0; i<taxon_length; i++) {
	list[i] = new file_parse<T>*[taxon_length];
	for(uint out = 0; out < taxon_length; out++) {//-1; out >-1; out--) {
	  list[i][out] = NULL; //file_parse<T>(i, out, _USE_BEST_BLAST_PAIR_SCORE); 
	}
      }
      return list;
    } else fprintf(stderr, "!!\tTaxon length not set: Discard initializing the list in 'list_file_parse' class found at line %d in file %s. Please contact oekseth@gmail.com if this mesage is seen.\n", __LINE__, __FILE__);
    return NULL;
  }

  /**
     @brief Cleans the memory for the private structures used in this;
     @param <remove_file> If set, removes the intermdediary file from disk (if existing).
  **/
  void finalize(bool remove_file) {
    index_list::close(listIndex);
    relations_list<T>::close(listRelations, remove_file);  
  }
  /**
     @brief Cleans the memory for the private structures used in this;
     @param <remove_file> If set, removes the intermdediary file from disk (if existing).
  **/
  void free_mem(bool remove_file) {finalize(remove_file);}

  //! Frees the reserved memory ontly for the non-inparalogs
  static void close_non_inparalogs(file_parse ***&list, const bool delete_file, uint taxon_start, uint taxon_end) {
    if(list) {
      for(uint i =taxon_start; i<taxon_end; i++) {
	if(list[i] !=NULL) {
	  for(uint k = taxon_start;k<taxon_end;k++) {
	    if(i!=k) { // It's not a possible inparalog-relation
	      if(list[i][k]) {
		list[i][k]->free_mem(delete_file);
		delete list[i][k]; list[i][k] = NULL; 
	      }
	    }
	  }
	}
      }
    }
  }
  //! Frees the reserved memory ontly for the inparalogs
  static void close_For_inparalogs(file_parse ***&list, const bool delete_file, uint taxon_start, uint taxon_end) {
    if(list) {
      for(uint i =taxon_start; i<taxon_end; i++) {
	if(list[i] !=NULL) {
	  if(list[i][i]) {
	    list[i][i]->free_mem(delete_file);
	    delete list[i][i]; list[i][i] = NULL; 
	  }
	}
      }
    }
  }
  //! Frees the reserved memory ontly for the non-inparalogs
  static void close(file_parse ***&list, const bool delete_file, uint taxon_start, uint taxon_end) {
    if(list) {
      for(uint i =taxon_start; i<taxon_end; i++) {
	if(list[i] !=NULL) {
	  for(uint k = taxon_start;k<taxon_end;k++) {
	    if(list[i][k]) {
	      list[i][k]->free_mem(delete_file);
	      delete list[i][k]; list[i][k] = NULL; 
	    }
	  }
	}
      }
    }
  }
  //! Deallocates the object given.
  static void close(file_parse ***&list, uint taxon_length, const bool delete_file) {
    if(list) {
      for(uint i =0; i<taxon_length; i++) {
	if(list[i] !=NULL) {
	  for(uint k = 0;k<taxon_length;k++) {
	    if(list[i][k]) {
	      list[i][k]->free_mem(delete_file);
	      delete list[i][k]; list[i][k] = NULL; 
	    }
	  }
	  delete [] list[i], list[i] = NULL;
	}
      }
      delete [] list; list = NULL;
    }
  }
  //! Deallocates the object given.
  static void close(file_parse *&obj, const bool delete_file) {
    if(obj) {
      obj->free_mem(delete_file);
      delete obj; obj = NULL;
    }
  }
  /**Initiates the variables, bu do not allocate memory.**/
 file_parse() : listIndex(NULL), listRelations(NULL), taxon_in(-1), taxon_out(-1),
   my_id(-1), USE_BEST_BLAST_PAIR_SCORE(false) {
  }
  /**Initiates the variables, bu do not allocate memory.**/
 file_parse(bool _USE_BEST_BLAST_PAIR_SCORE) : listIndex(NULL), listRelations(NULL), taxon_in(-1), taxon_out(-1),
   my_id(-1), USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE) {
  }
  //! Initializes the class, allocating memory for the objects residing in it.
   file_parse(uint in, uint out, bool _USE_BEST_BLAST_PAIR_SCORE) :
     taxon_in(in), taxon_out(out), my_id(-1), USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE) 
  {
    listIndex = index_list::init();
    listRelations = relations_list<T>::init(in, out, USE_BEST_BLAST_PAIR_SCORE);
  }

  //! Initializes the class, allocating memory for the objects residing in it.
   file_parse(int in, int out, int length_prots, int _block_number, mem_loc mem_reserved, char *FILE_BINARY_LOCATION, bool _USE_BEST_BLAST_PAIR_SCORE) :
     listIndex(NULL), listRelations(NULL), taxon_in(in), taxon_out(out), my_id(-1), USE_BEST_BLAST_PAIR_SCORE(_USE_BEST_BLAST_PAIR_SCORE)  
  {
    listIndex = index_list::init(length_prots);
    listRelations = relations_list<T>::init(in, out, _block_number, mem_reserved, FILE_BINARY_LOCATION, _USE_BEST_BLAST_PAIR_SCORE);
  }
  //! Test internal procedures.
  void assert_buffer_write_and_read() {
#ifdef assert_code  
    // Builds some data:
    uint insert_total = 100;//     uint insert_total = 10;
    for(uint index_in =0; index_in<insert_total; index_in++) {
      const uint i_out = (uint)(index_in-1);
      assert(false == insert_protein_index(index_in, i_out)); // Inserts starting positions
      // Asserts that the starting postion is correct:
      uint index_out = 0;   if(index_in == 0) index_out = 1;
      insert_new_rel(index_in, index_out, (float)index_in, index_in, index_in);    
    }    
    listRelations->dump_buffer();
    free_mem(true);
#endif
  }

  //! The main test function for this class  
  static void assert_class(const bool print_info) {
    const static char *class_name = "file_parse";
    if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code  
#ifndef NDEBUG  
    // Test 1:
    // - The basic operations:
    uint in=0, out=0, length_prots=2000, block_number=0, mem_reserved=1;
    const bool USE_BEST_BLAST_PAIR_SCORE = false;
    char *FILE_BINARY_LOCATION = ".";
    file_parse test_1 = file_parse(in, out, length_prots, block_number, mem_reserved, FILE_BINARY_LOCATION, false);
    // Validates internal file handling procedures, using intermediate files memory usage
    test_1.assert_buffer_write_and_read(); 

    test_1 = file_parse(in, out, length_prots, block_number, mem_reserved, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);
    assert(test_1.getLengthProteins() == (int)length_prots);
    assert(!test_1.has_data());
    // Build a simple list with each protein having one relation of distance and overlap:
    // -- This implies that the data checked will always reside in memory.

      for(uint index_in =0; index_in<length_prots; index_in++) {
	assert(false == test_1.insert_protein_index(index_in, (uint)(index_in-1))); // Inserts starting positions
	assert(false == test_1.is_equal(index_in, (uint)(index_in-1))); 
	// Asserts that the starting postion is correct:
	uint index_out = 0;   if(index_in == 0) index_out = 1;
	test_1.insert_new_rel(index_in, index_out, (float)index_in, index_in, index_in); 
    
	//    printf("inserted index %u\n", index_in);
	//    printf("-\tabs_start_pos(%u)\n", test_1.get_absolute_start_pos(index_in));
	//    printf("\nassert((1+index_in(%d))  =?= %d = test_1.get_absolute_start_pos(index_in))", index_in, test_1.get_absolute_start_pos(index_in));
	assert((1+index_in)  == test_1.get_absolute_start_pos(index_in));

	assert(index_in == test_1.get_index_start(index_in));
	assert(test_1.has_data(index_in));
	assert(1 == test_1.get_index_length(index_in));
	assert(test_1.has_data_slow(index_in, index_out));
      }


      for(uint index_in =0; index_in<length_prots; index_in++) {
	uint index_out = 0;   if(index_in == 0) index_out = 1;
	assert((1+index_in)  == test_1.get_absolute_start_pos(index_in));
	assert(test_1.has_data(index_in));
	assert(1 == test_1.get_index_length(index_in));
	assert(test_1.has_data_slow(index_in, index_out));
      }
      assert(test_1.has_data());
      test_1.free_mem(true);

      // ok
      //
      // Test 2: Increase of the memory
      in=0, out=0, length_prots=100, block_number=0, mem_reserved=200;
      test_1 = file_parse(in, out, length_prots, block_number, mem_reserved,FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);
      assert(test_1.getLengthProteins() == (int)length_prots);
      for(uint i =0; i<length_prots; i++) {
	assert(false == test_1.insert_protein_index(i, (uint)(i-1)));
      }
      test_1.free_mem(true);

      //int taxon_length = 3;
      const int length_in = 7, length_out = 2;
      test_1 = file_parse(in, out, length_prots, block_number, mem_reserved,FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);
      uint total_inserted = 0;
      for(int prot_in = 0; prot_in< length_in; prot_in++) {
	for(int prot_out = 0; prot_out< length_out; prot_out++) {
	  test_1.insert_new_rel(prot_in, prot_out, (float)prot_in, prot_in, prot_in);   
	  if(prot_in != prot_out) total_inserted++;
	}
      }
      //assert(test_1.buffer_in_mem_pos == total_inserted);

    if(true) {
      // TODO: Consider writing deeper tests.
      file_parse *test_2 = new file_parse(USE_BEST_BLAST_PAIR_SCORE);
      //      assert(test_2.is_empty(false));
      //      taxa temp_taxa = taxa();
      file_parse test_3 = file_parse(USE_BEST_BLAST_PAIR_SCORE);
      test_3.steal_data(test_2);
      assert(test_2->is_empty());
      test_1.free_mem(true);
      delete test_2; test_2 = NULL;
      //  if(false)   print_constants();
    }
#endif
#endif
    if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
  }



};

/**
   @brief Contains the blast file representation for a bunch of taxon-pairs.
   @ingroup blastfile_container
**/
//typedef class file_parse<p_rel> file_parse_t;
#endif
