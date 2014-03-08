#include "index_list.h"

/**
   @brief Inserts a new list of type index_t
   @remarks Assumes the internal list is empty.
**/
void index_list::insert_new_list(index_t *arr_new, const uint arr_new_size) {
  if(arr_new_size) {
    assert(index_used < 1); // No data inserted.
    if(arr_new_size > (uint)index_reserved) enlarge(arr_new_size);
    memcpy(list, arr_new, sizeof(index_t)*arr_new_size);
    index_used = index_reserved-1;
  }
}

//! @return true if the index size is inside the limits set of this structure
bool index_list::test_index_size(int ind) {
  if(list!= NULL)
    if(ind <= index_used)
      return true;
  return false;
}

//! @return the number of protein relations for the index given.
int index_list::get_index_length(int protein_in) {
  if(test_index_size(protein_in)) return list[protein_in].get_length();
  else return 0;
}

//! @return the number of protein relations for index_list object.
int index_list::get_index_length() {
  if(list) return index_used+1;
  else return 0;
}

//! @return the starting position in the rel object for the index given.
uint index_list::get_start_pos(int index) {
  return list[index].get_start_pos();
}

//! @return the absolute start-position in the rel object for the index given.
uint index_list::get_absolute_start_pos(int index) {
  const uint abs = list[index].get_absolute_start_pos();
  if(index <= index_used) return abs;
  else return UINT_MAX;
}

//! @return the starting position in the rel object for the index index, writing an error message if not found.
uint index_list::get_index_start(uint protein_in) {
  if(list != NULL && protein_in < (uint)index_reserved) {// has data stored
    return list[protein_in].get_start_pos();
  } else fprintf(stderr, "!!\tNot found index(%u) at line %d in class index_list: Contact the developers of the code if this error-message is seen.\n", protein_in, __LINE__);
  return 0;
}

//! @return the end position in the rel object for the index given.
uint index_list::get_end_pos(int index) {
  return list[index].get_end_point();
}

//! Sets the starting position in the rel object for the index given.
void index_list::set_start_pos(int index, mem_loc start_pos) {
  list[index].set_start_pos(start_pos);
}

//! Increments length for the given object specified by the index.
void index_list::increment_length(int index, int length) {
  list[index].increment_length(length);
}

//! Increments length for the given object specified by the index.
void index_list::increment_length(int index) {
  list[index].increment_length();
}

//! Decreases the given object id, given the param index.
void index_list::decrease_length(int index) {
  list[index].decrease_length();
}

//! Sets the number of relations in the rel object for the index given.
void index_list::set_length(int index, uint pos) {
  list[index].set_length(pos);
}

//! @return the number of relations in the rel object for the index given.
uint index_list::get_length(int index) {
  return list[index].get_length();
}


//! @return the total number of relations in the rel object of this list.
uint index_list::getTotalLengthOfData() {
  if(!list) return 0;    
  else {
    uint sum = 0;
    int length = index_used +1;
    if(length > index_reserved) length = index_reserved;
    for(int i = 0; i < length; i++) sum+=list[i].get_length();
    return sum;
  }
}

//! Prints information of the list of objects given as input.
void index_list::print_list(index_t *list, uint index_size) {
  for(uint i = 0; i < index_size; i++) list[i].print(i);
}

//! Prints information of the list of objects given as input.
void index_list::print_list() {
  for(int i = 0; i < index_reserved; i++) {
    if(list[i].get_length() != 0)  printf("\t[%d]\t .start = %u, .length=%u:\n", i, list[i].get_start_pos(), list[i].get_length());
  }
}

/**
   @brief Merges two indexes
   @param <destination> The index object to receive data.
   @param <destination_size> The current size of the receiving object.
   @param <source> The index object giving data.
   @param <source_size> The current size of the sending object.
   @param <offset> The offset to use when inserting into the destination object.
   @author Ole Kristian Ekseth (oekseth)
**/
void index_list::merge_lists(index_t *&destination, uint destination_size, index_t *source, uint source_size, uint offset) {
  if(destination_size >= source_size) {
    for(uint i = 0; i< source_size;i++) {
      if(source[i].get_length() > 0) {
	destination[i].set_start_pos(offset + source[i].get_start_pos());
	destination[i].set_length(source[i].get_length());
      }
    }
  } else fprintf(stderr, "!!\tAn error occured in merge(..) foudn in class 'index'\n");
}

/**
   @brief Enlarges the current buffer. Optional parameters may be set.
   @param <buffer> The buffer to be enlarged.
   @param <new_size> The new size of the buffer.
   @param <start_pos> The start position to set the old buffer (an optional param, normally set to '0').
   @param <copy_length> The length of the buffer to copy into the structure.
   @author Ole Kristian Ekseth (oekseth)
*/
void index_list::enlarge(index_t *&buffer, uint new_size, uint start_pos, uint copy_length) {
  index_t *buff_temp = new index_t[new_size]();//(index_t*)malloc(sizeof(index_t)*new_size);
  uint pos = start_pos;     if(copy_length > 1) pos+=copy_length-1;
  if(buff_temp) {  // Copys the old data into the new structure:
    for(uint i = pos; i < new_size; i++) buff_temp[i] = index_t();     // Initiates the buffer:
    // Inserts the buffer in the front of the file (i.e. in
    // order for the indexes to maintain their correct position pointers:    
    memcpy(buff_temp +start_pos, buffer, sizeof(index_t)*copy_length);
    delete [] buffer;
    buffer = buff_temp;
  } else {
    fprintf(stderr, "!!\tUnable to reserve memory for the index_t of protein-length(%u) object (ie. Probably implies that your memory chip is too small). If questions, please contact the developer at oekseth@gmail.com, giving the following information: This message was genereated at line %d in file %s, found in method %s\n",
	    new_size, __LINE__, __FILE__, __FUNCTION__);
    exit(2); // No point in continuing the work flow.
  }
}


/**
   @brief Enlarges the current buffer. Optional parameters may be set.
   @param <new_size> The new size of the buffer.
   @author Ole Kristian Ekseth (oekseth)
*/
void index_list::enlarge(uint new_size) {
  // TODO: valider parametrne bruke her!
  if((int)new_size > (int)(index_reserved-1)) {
    uint temp_size = new_size*2;
    if(new_size < 100) temp_size = new_size + 10; // TODO: Validate that this extra line decreases the memory consumption
    enlarge(list, temp_size, 0, index_reserved);
    index_reserved = temp_size;
  }
}

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
bool index_list::insert(int index_in, int index_out, mem_loc pos_in_rel_struct) {
  enlarge(index_in);
  if(index_in > index_used) {
    index_used = index_in;
  }
  if(index_in_prev != index_in){ // a protein on the left (inner) side
    list[index_in].set_start_pos(pos_in_rel_struct);
    return false;
  } else if(index_out_prev == index_out) {
    return true;
  }
  return false;
}


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
void index_list::assert_merging_of_index_list(index_t *listArgument, const loint listArgument_size, loint offset_to_add_data_to, index_t *&root_list_copy) {
#ifndef NDEBUG
  assert(list);
  if(listArgument && listArgument_size) {
    if(!root_list_copy) { // The data is given before the merging; copies the root:
      root_list_copy = init_list(index_reserved);
      assert(root_list_copy); // Has memory reserved.
      // The copying doing it the slow, but easy to verify-process:
      for(uint i = 0; i < (uint)index_reserved; i++) root_list_copy[i].set_data(list[i]);
    } else { // Verifies that the two lists if found in the new list in this object:
      bool all_are_equal = true;
      // Tests for the root being the same:
      for(uint i = 0; i < (uint)index_reserved; i++) {
	// As the data is assumed beeing disjoint, 
	if(root_list_copy[i].has_data()) {
	  // only tests for those data set that they have'nt changed:
	  if(!list[i].is_greater_or_equal(root_list_copy[i])) {
	    fprintf(stderr, "[%u]\tcopy{start(%u), length(%u)} != =< list{start(%u), length(%u)}\n",
		    i, root_list_copy[i].get_start_pos(), root_list_copy[i].get_length(),
		    list[i].get_start_pos(), list[i].get_length());
	    all_are_equal = false;
	  }
	}
      }
      assert(all_are_equal); 
      delete_list(root_list_copy);
      assert(!root_list_copy);
      // Tests for the argument given changed according ot the rules:
      for(uint i = 0; i < (uint)listArgument_size; i++) {
	// As the data is assumed beeing disjoint, 
	if(listArgument[i].has_data()) {
	  // only tests for those data set that they have'nt changed:
	  if(list[i].get_length() < listArgument[i].get_length()) { // TODO: validate this alternative!
	       //	  if(!list[i].is_greater_or_equal(listArgument[i], offset_to_add_data_to)) {
	    fprintf(stderr,"[%u]\tlistArgument{(start(%u)+offset(%llu)),length(%u)} != =< list{start(%u),length(%u)}\n",
		    i, listArgument[i].get_start_pos(), offset_to_add_data_to, listArgument[i].get_length(),
		    list[i].get_start_pos(), list[i].get_length());
	    all_are_equal = false;
	  }
	}
      }
      assert(all_are_equal);
      //	fprintf(stderr, "Tests have passe in index_list.h\n");
    }
  }
#else
  if(false) { // Simple method in order to avoid the compiler complaining:
    if(listArgument && root_list_copy) {
      assert(listArgument_size);
      assert(offset_to_add_data_to != -1);
    }
  }
#endif
}
/**
   @brief Merges the input with this, without affecting the input.
   @author Ole Kristian Ekseth (oekseth)
   @return true if the inserted buffer was a duplicate.
**/
bool index_list::merge_buffers(index_list *listArgument, const loint listArgument_size, loint offset_to_add_data_to, taxa &taxa_obj, taxa *taxa_obj_out) { 
  bool is_duplicate_of_previous = false;
  assert(listArgument != NULL);
  if(!(listArgument_size <= (loint)index_reserved)) { // ensures that the size of the internal index is correct.
    enlarge(listArgument_size);  
  }
  assert(listArgument_size <= (loint)index_reserved);
  // in order to know the starting postion of the last inserted proteins' data:
  loint last_inserted = -1; 
  // Prepares verification of the merging process:
  index_t *root_list_copy = NULL;
  assert_merging_of_index_list(listArgument->get_list(), listArgument_size, offset_to_add_data_to, root_list_copy);
  //
  // iterates over the listArgumentuments proteins:
  uint cnt_rows_with_error = 0;
  for(int i = 0; i< listArgument->get_index_reserved();i++) { 
    static bool printed_warning = false;
    static uint cnt_skipped = 0;
    if (listArgument->has_data(i)) {
      if(!has_data(i)) {
	const loint pos_in_mem_for_rel = offset_to_add_data_to + listArgument->get_index_start(i);
	set_start_pos(i, pos_in_mem_for_rel);
	set_length(i, listArgument->get_length(i));
      } else {
	bool is_an_error = true;
// 	if(listArgument->get_length(i) == list[i].get_length()) {
// 	  //! Then there might be an error in the parsing, which we intend fixing.
// 	  // TODO: validate correctness of this operation.
// 	  //s_duplicate_of_previous = true;
// 	  is_an_error = false;
// 	}
	cnt_rows_with_error++;

	if(is_an_error) {
	  if(!printed_warning) { // data already registred: enllistArgumentes the list:
	    
	    //	const uint length_out_i = listArgument->get_length(i);
	    //	increment_length(i, length_out_i);
	    // TODO: when error is fixed, move printint to "log_builder"!
	    printed_warning = true;
	    fprintf(stderr, 
		  "!!\t Seems like you used too many CPU's. Our first impression is that there is no gain in using the number of CPU's you are currently using. (There are other explanations, but this is the most commonplace.) This first impression (i.e. our latter thoughts) might be wrong. If you are puzzled about this message, it might be that me (i.e. orthAgogue) has made an error in my calculation. With the dabger of making you even more puzzled, here is my initial thoughts:\n"
		    "(1)\t Try to reduce the number of CPUs used. For an extreme case, try setting either \"-c 1\" or \"--cpu 1\" (which are synonyms) on you command line orthAgogue execution.\n"
		    "#\t If this does not make the day (i.e. the problem is persisting), then there might be an error in the merging of mapping-data for a protein (programmatically speaking when merging blocks of type \"index_t\").\n"
		  "#\t An alternative explanation regards the merging procedure applied by your Blast file: if the BlastP did a correct concatenation, then protein\"%s\" should have received all relations for a given taxa(\"%s\"-->\"%s\"), before receiving the next combination of pairs. We hope the implications of this error is low, as we have kept the best high-scoring pair in question.\n"
		    "-\t In brief, we expect the received protein block to cover all comparisons for \"%s\".\n"
		    "-\t If this error is seen (and point (1) above did not bring the nuts), please contact the developer, either through our web-page, directly at [oekseth@gmail.com], or try:\n"
		    "\t(2)\t(a)\t install the tool using the \"install_debug.bash\" script.\n"
		    "\t(2)\t(b)\t Experiment with different values for parameter \"--disk_buffer_size\".\n"
		    "\t(2)\t(c) If this does not solve the case, please report the error (to the 'issue' page of orthAgogue's homepage, or contact the developer directly at [oekseth@gmail.com]).\n"
		    "This error was produced at [%s]:%s:%d\n",
		    taxa_obj.getArrKey(i), taxa_obj.get_name(), taxa_obj_out->get_name(),
		    taxa_obj.getArrKey(i), __FUNCTION__, __FILE__, __LINE__);
	    //assert(false);
	  } 
	  cnt_skipped += (uint)listArgument->get_length(i);
#ifndef NDEBUG
	  assert(taxa_obj_out);      
	  if(cnt_rows_with_error ==1) {fprintf(stderr,  "[cnt-skipped=%u]\t Received %u pairs for protein[%u]=\"%s\", given taxa(\"%s\"-->\"%s\"), though it already had %u elements. We therefore skip the %u proteins. Message at [%s]:%s:%d\n", 
		cnt_skipped, (uint)listArgument->get_length(i), (uint)i, taxa_obj.getArrKey(i), 
		  taxa_obj.get_name(), taxa_obj_out->get_name(), 
		  (uint)list[i].get_length(), (uint)listArgument->get_length(i), 
					       __FUNCTION__, __FILE__, __LINE__);}
#endif
	}
      }
      if(get_start_pos(i) > last_inserted) { 	// Sets the index whos last visited
	last_inserted = get_start_pos(i);
	index_in_prev = i; // in order to know the id of the last inserted protein
      }
    } // end adding position data for protein 'i'
  }
  if(cnt_rows_with_error > 0) {
    fprintf(stderr, "!!\t %u out ouf %u proteins had errors, at index_list:%d\n", cnt_rows_with_error, (uint)listArgument->get_index_reserved(), __LINE__); 
  }

  if(is_duplicate_of_previous == false) {
    // Verifies the merging:
    assert_merging_of_index_list(listArgument->get_list(), listArgument_size, offset_to_add_data_to, root_list_copy);
    // Updates the index:
    if(index_used < listArgument->get_index_used()) index_used = listArgument->get_index_used();  
  }
  return is_duplicate_of_previous;
}


/**
   @return The total number of elements the lists represents in the file fo the span given.
*/
uint index_list::get_number_of_elements(uint start, uint end) {
  if(list) {
    if((index_reserved > 0) && ((int)end <= index_reserved)) { // has data set
      uint sum_outer = 0;
      for(uint protein_in = start; protein_in < end; protein_in++)
	sum_outer += list[protein_in].get_length();
      return sum_outer;
    }
  }
  return 0;
}

/**
   @return true is data is set for the params.
**/
bool index_list::has_data(index_t *list, uint index_in, uint size) {
  if(list != NULL) {// has data stored
    if(index_in < size) {
      if(list[index_in].get_length() > 0) return true;
      else return false;
    }
  }
  return false;
}

/**     @return true is data is set for the params.  **/
bool index_list::has_data(uint index_in) {
  return has_data(list, index_in, index_reserved);
}

/** @return true if the index object of this class equals the memory reference of the input.  **/
bool index_list::pointer_equals(index_t *ind) {
  return (ind == list);
}

/**
   @param <in_prev> The index to test for the inner case.
   @param <out_prev> The index to test for the outer case.
   @return true if thei are equal the data set in the object.
   @author Ole Kristian Ekseth (oekseth)
**/
bool index_list::index_are_equal(int in_prev, int out_prev) {
  if(in_prev == index_in_prev)
    if(out_prev == index_out_prev)
      return true;
  return false;
}

/**
   @param <arg> the object to compare this with.
   @param <print_info> if set to true prints info about in-equalities.
   @return true if the input object is equal this.
**/
bool index_list::is_equal(index_list *arg, bool print_info) {
  if(arg) {
    bool is_equal = true;
//     if(!arg->pointer_equals(list)) {
//       if(print_info) printf("%s not set, ", "index");
//       is_equal= false;
//     }
    if(index_used != arg->get_index_used())     {
      if(print_info) printf("%s not set, ", "index_used");
      is_equal= false;
    }
//     if(index_reserved != arg->get_index_reserved())    {
//       if(print_info) printf("%s not set, ", "index_reserved");
//       is_equal= false;
//     }
//     if(index_in_prev != arg->get_index_in_prev()) {
//       if(print_info) printf("%s not set, ", "index_in_prev");
//       is_equal= false;
//     }
//     if(index_out_prev != arg->get_index_out_prev()) {
//       if(print_info) printf("%s not set, ", "index_out_prev");
//       is_equal= false;
//     }
    return is_equal;
  } else return false;
}


/**
   @brief Sums the memory consumption for the object.
   @remarks For anaylsing the memory fingerprint.
   @return the total length of data in Bytes.
**/
loint index_list::getTotal_memoryConsumption_for_object() {
  loint sum = sizeof(index_list);
  sum += index_reserved * sizeof(index_t);
  return sum;
}


/**
   @brief Gets the object size.
   @remarks Uses the deault size-init-values of the objects of type index
   @return the size of an empty object.
**/
loint index_list::get_empty_object_size() {
  loint total = sizeof(index_list); // Do not allocate any objects of type 'index'
  return total;
}

/**
   @brief Builds an initiated list of index objects.
   @param <index_size> The size of the list to return.
   @return the list.
**/
index_t *index_list::init_list(uint index_size) {
  if(index_size) {
    index_t *list = new index_t[index_size];
    if(list) return list;
    else {
      fprintf(stderr, "!!\tUnable to reserve memory for the index_t of protein-length(%u) object (ie. Probably implies that your memory chip is too small). If questions, please contact the developer at oekseth@gmail.com, giving the following information: This message was genereated at line %d in file %s, found in method %s\n",
	      index_size, __LINE__, __FILE__, __FUNCTION__);
      exit(2); // No point in continuing the work flox.
    }

    //(index_t*)malloc(sizeof(index_t)*index_size);
    //  for(uint i = 0; i<index_size; i++) list[i] = index_t();
  } else return NULL;
}


//! Deallocates the memory for the object given.
void index_list::delete_list(index_t *&list) {
  delete [] list; list = NULL;
}

//! Deallocates the memory for this object.
void index_list::free_mem() {
  if(list){delete [] list; list = NULL;}
  index_reserved = 0;
}

//! Deallocates the memory for the object given.
void index_list::close(index_list *&obj) {
  if(obj) {
    obj->free_mem();
    delete obj; obj = NULL;
  }
}

//! The constructor.
index_list::index_list() :
  //   printed_warning(false), 
  list(NULL), index_reserved(0),index_used(0), index_in_prev(-1),
  index_out_prev(-1)
{

}

//! The constructor.
index_list::index_list(uint _index_reserved, uint _index_used) :
//   printed_warning(false), 
  index_reserved(_index_reserved), index_used(_index_used),
  index_in_prev(-1), index_out_prev(-1)
{
  if(index_reserved) list = init_list(index_reserved);
  else list = NULL;
}
void index_list::assert_class(const bool print_info) {
  const static char *class_name = "index_list";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  // TODO: Skriv noen tester....
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
