#ifndef index_list_h
#define index_list_h
/**
   @file
   @brief Maps a protein-id to a set of other proteins (i.e. identifies the protein pairs).
   @ingroup parsing_ops
   @date 04.01.2012 by oekseth (cleanup).
   @author Ole Kristian Ekseth (oekseth).
**/ 

#include "mcl_format.h"
#include "index.h"

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
#include "../configure.h"
/**
   @class index_list
   @brief Maps a protein-id to a set of other proteins (i.e. identifies the protein pairs).
   @ingroup parsing_ops
   @remark Holds the index_list used in retrieving data for a spesific protein after read from a binary file.
   -- Each index_list correponds to a specific protein id.
   @author Ole Kristian Ekseth (oekseth).
   @date 18.03.2011 by oekseth (initial).
   @date 27.08.2011 by oekseth (asserts).
   @date 04.01.2012 by oekseth (cleanup).
**/ 
class index_list { 
 private:
  //  bool printed_warning;
  index_t *list;
  int index_reserved;
  int index_used;
  int index_in_prev, index_out_prev;
  //! @return true if the index size is inside the limits set of this structure
  bool test_index_size(int ind);
 public:
  /**
     @brief Inserts a new list of type index_t
     @remarks Assumes the internal list is empty.
  **/
  void insert_new_list(index_t *arr_new, const uint arr_new_size);
  //! @return the number of used objects for type index.
  int get_index_used(){return index_used;}
  //! @return the number of reserved objects for type index.
  int get_index_reserved(){return index_reserved;}
  //! @return the last inserted index value for the leftmost (inner) relation.
  int get_index_in_prev(){return index_in_prev;}
  //! @return the last inserted index value for the rightmost (outer) relation.
  int get_index_out_prev(){return index_out_prev;}
  //! Sets the previous index for the inner index.
  void set_index_in_prev(int ind) {index_in_prev = ind;}
  //! Sets the previous index for the outer index.
  void set_index_out_prev(int ind) {index_out_prev = ind;}
  //! @return the pointer to the index object of this class.
  index_t *get_list(){return list;}
  //! @return the number of protein relations for the index given.
  int get_index_length(int protein_in);
  //! @return the number of protein relations for index_list object.
  int get_index_length();
  //! @return the starting position in the rel object for the index given.
  uint get_start_pos(int index);
  //! @return the absolute start-position in the rel object for the index given.
  uint get_absolute_start_pos(int index);
  //! @return the starting position in the rel object for the index index, writing an error message if not found.
  uint get_index_start(uint protein_in);
  //! @return the end position in the rel object for the index given.
  uint get_end_pos(int index);
  //! Sets the starting position in the rel object for the index given.
  void set_start_pos(int index, mem_loc start_pos);
  //! Increments length for the given object specified by the index.
  void increment_length(int index, int length);
  //! Increments length for the given object specified by the index.
  void increment_length(int index);
  //! Decreases the given object id, given the param index.
  void decrease_length(int index);
  //! Sets the number of relations in the rel object for the index given.
  void set_length(int index, uint pos);
  //! @return the number of relations in the rel object for the index given.
  uint get_length(int index);
  //! @return the total number of relations in the rel object of this list.
  uint getTotalLengthOfData();
  //! Prints data for the index given.
  void print_index(int index) {list[index].print(index);}
  //! Prints data for the index given.
  void print_index(int index, char *lbl) {list[index].print(lbl);}
  //! Prints information of the list of objects given as input.
  void print_list(index_t *list, uint index_size);
  //! Prints information of the list of objects given as input.
  void print_list();
  
  /**
     @brief Merges two indexes
     @param <destination> The index object to receive data.
     @param <destination_size> The current size of the receiving object.
     @param <source> The index object giving data.
     @param <source_size> The current size of the sending object.
     @param <offset> The offset to use when inserting into the destination object.
     @author Ole Kristian Ekseth (oekseth)
  **/
  void merge_lists(index_t *&destination, uint destination_size, index_t *source, uint source_size, uint offset);  
  /**
     @brief Enlarges the current buffer. Optional parameters may be set.
     @param <buffer> The buffer to be enlarged.
     @param <new_size> The new size of the buffer.
     @param <start_pos> The start position to set the old buffer (an optional param, normally set to '0').
     @param <copy_length> The length of the buffer to copy into the structure.
     @author Ole Kristian Ekseth (oekseth)
  */
  void enlarge(index_t *&buffer, uint new_size, uint start_pos, uint copy_length);
  /**
     @brief Enlarges the current buffer. Optional parameters may be set.
     @param <new_size> The new size of the buffer.
     @author Ole Kristian Ekseth (oekseth)
  */
  void enlarge(uint new_size);
  /**
     @brief Sets the start position in the index object.
     @param <index_in> The leftmost index found in the blast file.
     @param <index_out> The outermost index found in the blast file.
     @param <pos_in_rel_struct> The position of the first relation this index will refer to.
     @returns true if the index_out was the first outer relation used for the index_in.
     @remarks Updates the index object only if this class have not seen the index_in before.
     @todo Consider changing the name of this function as it only sets the start postion.
     @author Ole Kristian Ekseth (oekseth)
  */
  bool insert(int index_in, int index_out, mem_loc pos_in_rel_struct);

  /**
     @brief Verifies that the mergin of an input index with this is performed correctly
     @remarks
     - Only active if NDEBUG is not defined, ie, performs without bread if NDEBUG is not defined.
     - Run two times: (First) Copies the root. (Last) Compares the indexes to the expectations.
     @param <listArgument> The data to update this data with.
     @param <listArgument_size> The new size of of the object.
     @param <offset_to_add_data_to> The incremental offset to use when the starting position the index object refers are are updated.
     @param <root_list_copy> The data is given before the merging, therefore copies the root
     @author Ole Kristian Ekseth (oekseth)
   **/
  void assert_merging_of_index_list(index_t *listArgument, const loint listArgument_size, loint offset_to_add_data_to, index_t *&root_list_copy);


  /**
     @brief Merges the input with this, without affecting the input.
     @param <listArgument> The data to update this data with.
     @param <listArgument_size> The new size of of the object.
     @param <offset_to_add_data_to> The incremental offset to use when the starting position the index object refers are are updated.
     @param <taxa_obj> Used for improved error-message-generation.
   @return true if the inserted buffer was a duplicate.
     @author Ole Kristian Ekseth (oekseth)
  **/
  bool merge_buffers(index_list *listArgument, const loint listArgument_size, loint offset_to_add_data_to, taxa &taxa_obj, taxa *taxa_obj_out); 
  //! @return The total number of elements the lists represents in the file fo the span given.
  uint get_number_of_elements(uint start, uint end);
  /**     @return true is data is set for the params.  **/
  bool has_data(index_t *list, uint index_in, uint size);
  /**     @return true is data is set for the params.  **/
  bool has_data(uint index_in);
  /**
     @brief sets an object of type index;
     @param <index_in> The index to get data from.
     @param <block> The object to set values into.
     @remarks Transfering data, sets the last param with it. Performs validation
   **/
  void get_block(uint index_in, index_t &block) {
    if(has_data(index_in)) list[index_in].get_data(block);
    else block = index_t();
  }
  /** @return true if the index object of this class equals the memory reference of the input.  **/
  bool pointer_equals(index_t *ind);
  /**
     @param <in_prev> The index to test for the inner case.
     @param <out_prev> The index to test for the outer case.
     @return true if thei are equal the data set in the object.
     @author Ole Kristian Ekseth (oekseth)
  **/
  bool index_are_equal(int in_prev, int out_prev);
  /**
     @param <arg> the object to compare this with.
     @param <print_info> if set to true prints info about in-equalities.
     @return true if the input object is equal this.
  **/
  bool is_equal(index_list *arg, bool print_info);

  /**
     @brief Sums the memory consumption for the object.
     @remarks For anaylsing the memory fingerprint.
     @return the total length of data in Bytes.
  **/
  loint getTotal_memoryConsumption_for_object();

  /**
     @brief Gets the object size.
     @remarks Uses the deault size-init-values of the objects of type index
     @return the size of an empty object.
  **/
  static loint get_empty_object_size();


  /**
     @brief Builds an initiated list of index objects.
     @param <index_size> The size of the list to return.
     @return the list.
  **/
  index_t *init_list(uint index_size);
  //! Return a pointer to an initiated object of type index_list.
  static index_list *init(){return new index_list();}
  //! Return a pointer to an initiated object of type index_list.
  static index_list *init(uint size){return new index_list(size, 0);}
  //! Deallocates the memory for the object given.
  void delete_list(index_t *&list);
  //! Deallocates the memory for this object.
  void free_mem();
  //! Deallocates the memory for the object given.
  static void close(index_list *&obj);
  //! The constructor.
  index_list();
  //! The constructor.
  index_list(uint _index_reserved, uint _index_used);
  /**
     @brief The main test function for this class.
     @remarks This method includes:
     - Formalized tests for setting of starting position and end position.
     - Valgrind used verifying memory usage.
     - Examples of how this class may be used.
     @author Ole Kristian Ekseth (oekseth)
     @date 04.01.2012 by oekseth.
  **/
  static void assert_class(const bool print_info);

}; 

/**
   @brief Maps a protein-id to a set of other proteins (i.e. identifies the protein pairs).
   @ingroup parsing_ops
*/
typedef class index_list index_list_t;

#endif
