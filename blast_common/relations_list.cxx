#include "relations_list.h"

//! Inserts an element into the list:
template <> void relations_list<p_rel_t>::insert(uint index_out, float sim_score, overlap_t overlap_in, overlap_t overlap_out) {
  enlarge();
  assert(buffer_in_mem_pos < buffer_in_mem_size);
  list[buffer_in_mem_pos].set_data(index_out,sim_score, overlap_in, overlap_out);
  buffer_in_mem_pos++;
}

//! Inserts an element into the list:
template <> void relations_list<rel_t>::insert(uint index_out, float sim_score, overlap_t overlap_in, overlap_t overlap_out) {
  if(false && !overlap_in && !overlap_out) ; // Hiding the variables.
  enlarge();
  assert(buffer_in_mem_pos < buffer_in_mem_size);
  list[buffer_in_mem_pos] = rel(index_out,sim_score);
  buffer_in_mem_pos++;
}



//! Dumps the buffer to a file, not resetting the data
template <>  void relations_list<p_rel_t>::dump_buffer() {
  if(list && buffer_in_mem_pos > 0) {
    //      assert(buffer_in_mem_pos);      
    p_rel::write_buffer(list, file_name, buffer_in_mem_pos, file_cnt_structs);
    p_rel::init(buffer_in_mem_pos, list);
    // Updates the position indicator:
    file_cnt_structs += buffer_in_mem_pos;
    buffer_in_mem_pos = 0;
  }
}



//! Dumps the buffer to a file, not resetting the data
template <>  void relations_list<rel_t>::dump_buffer() {
  if(list && buffer_in_mem_pos > 0) {
    //      assert(buffer_in_mem_pos);      
    rel::write_buffer(list, file_name, buffer_in_mem_pos, file_cnt_structs);
    rel::init(buffer_in_mem_pos, list);
    // Updates the position indicator:
    file_cnt_structs += buffer_in_mem_pos;
    buffer_in_mem_pos = 0;
  }
}

/**
   @brief Uses a p_rel object, thereby calling this method for another type causes error!
   @return a basket_parse object the buffer assuming it's in memory,
   @remarks Returning object used for transfering data between pipes.
**/
template <> class basket_parse *relations_list<p_rel>::get_buffer_basket() {
  //    if(list == NULL && file_cnt_structs > 0) {
  if (file_cnt_structs > 0) {
    if(list == NULL) {
      p_rel *buff = p_rel::get_file_block(file_name, 0, file_cnt_structs);
      basket_parse *basket = basket_parse::init(buff, file_cnt_structs, true);
      return basket;
    } else {
      fprintf(stderr, "!!\tList is not empty, and file has %llu amount of data in the file: Only returns the list residing in memory (of size %llu). Error at line %d in class relations_list.\n", file_cnt_structs, buffer_in_mem_pos, __LINE__);
      fprintf(stderr, "!!Only returns the data residing in memory: probably a bug. Found at line %d in file '%s'. Please report this to the developer, providing the sample input and machine architecture!!\n", __LINE__, __FILE__);
    }
  }
  // If at this point, all the data resides in memory:
  basket_parse *basket = basket_parse::init(list, buffer_in_mem_pos, false);
  return basket;
}



/**
   @brief Uses a rel object, thereby calling this method for another type causes error!
   @return a basket_parse object the buffer assuming it's in memory,
   @remarks Returning object used for transfering data between pipes.
**/
/*
template <> class basket_parse *relations_list<rel>::get_buffer_basket() {
  //    if(list == NULL && file_cnt_structs > 0) {
  if (file_cnt_structs > 0) {
    if(list == NULL) {
      rel *buff = rel::get_file_block(file_name, 0, file_cnt_structs);
      basket_parse *basket = basket_parse::init(buff, file_cnt_structs, true);
      return basket;
    } else {
      fprintf(stderr, "!!\tList is not empty, and file has %llu amount of data in the file: Only returns the list residing in memory (of size %llu). Error at line %d in class relations_list.\n", file_cnt_structs, buffer_in_mem_pos, __LINE__);
      fprintf(stderr, "!!Only returns the data residing in memory: probably a bug. Found at line %d in file '%s'. Please report this to the developer, providing the sample input and machine architecture!!\n", __LINE__, __FILE__);
    }
  }
  // If at this point, all the data resides in memory:
  basket_parse *basket = basket_parse::init(list, buffer_in_mem_pos);
  return basket;
}
*/
//! Does not do anything, except hyding complexity in the calling file_parse<T> object.
template <> void relations_list<rel>::increment_data(float sim_score, overlap_t overlap_in, overlap_t overlap_out) {
  if(false && !sim_score && !overlap_in && !overlap_out) ; // In order to hide the variables  
}


//! Increments the data using the last position
template <> void relations_list<p_rel>::increment_data(float sim_score, overlap_t overlap_in, overlap_t overlap_out) {
  loint prev_pos  = buffer_in_mem_pos;
  if(buffer_in_mem_pos > 0) prev_pos -=1;
  assert(prev_pos < buffer_in_mem_size);
  list[prev_pos].increment_data(sim_score, overlap_in, overlap_out, USE_BEST_BLAST_PAIR_SCORE);  
}
