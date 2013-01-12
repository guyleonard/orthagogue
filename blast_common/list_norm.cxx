#include "list_norm.h"

/**
   @brief Sets the data directly
   @remarks To be used with carefullness, as dat is not reserved, ie, function "free_mem()" must not be called.
**/
void list_norm::set_data(norm ***_list, int _taxon_length, float _max_input_value) {
  list = _list, taxon_length = _taxon_length, max_input_value = _max_input_value;
}

//! Sets the DEBUG_NORM value
void list_norm::set_DEBUG_NORM(const bool debug) {
  DEBUG_NORM=debug;
}
/**
   @brief Prints the memory consumption for the object.
   @remarks For anaylsing the memory fingerprint.
**/
void list_norm::print_getTotal_memoryConsumption_for_object(FILE *f) {
  assert(f);
  loint numbers_reserved = 0;
  loint pointers_reserved = 0;
  if(list) {
    for(uint i=0; i<(uint)taxon_length; i++) {
      if(list[i]) {
	pointers_reserved += taxon_length;
	for(uint out=0; out<(uint)taxon_length; out++) {
	  if(list[i][out]) {
	    numbers_reserved++;
	  }
	}
      }
    }
  }
  const loint size_reserved = sizeof(list_norm)
    + sizeof(norm_t)*numbers_reserved   // The object itself.
    + sizeof(int)   *pointers_reserved  // The "inner" pointers.
    + sizeof(int)*taxon_length;         // The "main" pointer.
  const float avg_cnt = (float)numbers_reserved/(taxon_length);
  //    const float fill_factor = (float)size_used / size_reserved;
  const float size_gb = (float)size_reserved/(1024*1024*1024);
  char end_char = ' ';
  // PRints the data, but in order to avoid clutter, newline is avoided if the goal is printing the data to the log file:
  if((f==stdout) || (f == stderr)) end_char = '\n';
  fprintf(f, "- Has %lld B of memory =~ %.5f GB => averaged_cnt(%.4f)%c", size_reserved, size_gb, avg_cnt, end_char); 
}

//! Test if condidtion is passed.
void list_norm::test_condition_and_if_not_abort(bool condition, const int line, const char *file, const char *function) {
  if(!condition) {
    fprintf(stderr, "!!\tWas not able allocating data due to memory constraints. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", line, file, function);
    exit(2);
  }
}

//! Inserts the relations into the structure using local coordinates
void list_norm::insert(uint taxon_in, uint taxon_out, float arg, uint protein_in, uint protein_out, struct taxa *listTaxa){
  insert(taxon_in, taxon_out, arg, protein_in, protein_out, listTaxa, false);
}
//! Inserts the relations into the structure using local coordinates
void list_norm::insert(uint taxon_in, uint taxon_out, float arg, uint protein_in, uint protein_out, struct taxa *listTaxa, const bool prot_pair_seen_before) {
  assert(list);
  assert(list[taxon_in]);
  if(!list[taxon_in][taxon_out]) list[taxon_in][taxon_out] = new norm();
  //  if(!list[taxon_in][taxon_out]) list[taxon_in][taxon_out] = new norm(taxon_in, taxon_out, taxon_length, PRINT_NORMALIZATION_BASIS, DEBUG_NORM, max_input_value);
  test_condition_and_if_not_abort(list[taxon_in][taxon_out], __LINE__, __FILE__, __FUNCTION__);
  list[taxon_in][taxon_out]->insert(taxon_in, taxon_out, arg, protein_in, protein_out, listTaxa, prot_pair_seen_before, DEBUG_NORM);    
}

//! Returns true if data is set.
bool list_norm::has_data(int in) {
  if((in < taxon_length) && list && list[in]) return true;
  else return false;      
}

//! @return the object at the given position: If not set returns NULL
norm_t *list_norm::get_list_element(int in, int out) {
  if(in < taxon_length && out < taxon_length) {
    if(list && list[in]) return list[in][out];
  }
  return NULL;
}

//! @return the object at the given position, and removes the local reference: If not set returns NULL
norm_t *list_norm::get_list_element_and_remove_reference(int in, int out) {
  if(in < taxon_length && out < taxon_length) {
    if(list && list[in]) {
      norm_t *element = list[in][out];
      list[in][out] = NULL; // Removes the local reference, making it safe to 'delete'
      return element;      
    }
  }
  return NULL;
}

//! Prints the element given the coordinates
void list_norm::print_element(uint in, uint out) {
  assert(in < (uint)taxon_length);
  assert(out < (uint)taxon_length);
  if(list && list[in]) {
    if(list[in][out]) list[in][out]->print(in, out);
  }
}

//! Print data given the params.
void list_norm::print_basis() {
  if(list) {
    for(uint i=0; i<(uint)taxon_length; i++) {
      if(list[i]) {
	for(uint out=0; out<(uint)taxon_length; out++) {
	  if(list[i][out]) list[i][out]->print(i, out);
	}
      }
    }
  }
}

//! Sums the values for the list of objects, updating the input parameters.
void list_norm::get_object_sum(loint &cnt, float &sum, uint &cnt_zero) {
  cnt = 0, sum = 0.0, cnt_zero = 0;
  if(list) {
    for(uint i=0; i<(uint)taxon_length; i++) {
      if(list[i]) {
	for(uint out=0; out<(uint)taxon_length; out++) {
	  if(list[i][out]) {
	    cnt      += list[i][out]->get_cnt();
	    sum      += list[i][out]->get_sum(max_input_value);
	    cnt_zero += list[i][out]->get_cnt_zeros();
	  }
	}
      }
    }
  }
}

//! Sets the total sum for the object given:
void get_total_score(list_norm &obj, int taxon_length, uint &cnt, float &sum, uint &cnt_zero, const float max_input_value) {
  cnt = 0, sum = 0, cnt_zero = 0;
  for(uint in = 0; in<(uint)taxon_length; in++) {
    for(uint out = 0; out<(uint)taxon_length; out++) {
      norm_t *element = obj.get_list_element(in, out);
      if(element) {
	cnt      += (uint)element->get_cnt();
	sum      += element->get_sum(max_input_value);
	cnt_zero += (uint)element->get_cnt_zeros();
      }
    }
  }
}
//! Merges two objects
void list_norm::merge_basis(list_norm &obj, uint i, uint out) {
  assert(list);
  assert(list[i]);
  norm_t *element = obj.get_list_element(i, out);
#ifndef NDEBUG
  loint cnt = 0, cnt_arg = 0;
  float sum = 0, sum_arg = 0;
  uint cnt_zero = 0, cnt_zero_arg = 0;
  if(list[i][out]) {
    cnt      = list[i][out]->get_cnt();
    sum      = list[i][out]->get_sum(max_input_value);
    cnt_zero = list[i][out]->get_cnt_zeros();
  }
  if(element) {
    cnt_arg      = element->get_cnt();
    sum_arg      = element->get_sum(max_input_value);
    cnt_zero_arg = element->get_cnt_zeros();
  }
#endif
  if(element) {      
    if(list[i][out]) {
      list[i][out]->merge(*element);
    }  else {
      list[i][out] = obj.get_list_element_and_remove_reference(i, out);
    }
  }  // else no point in merging if it does not have any data.	

#ifndef NDEBUG
  const int cnt_sum = cnt + cnt_arg;
  const float sum_sum = sum + sum_arg;
  const uint sum_zero = cnt_zero + cnt_zero_arg;
  if(list[i][out]) {
    const int   diff_cnt = cnt_sum - list[i][out]->get_cnt();
    const float diff_sum = sum_sum - list[i][out]->get_sum(max_input_value);
    const float diff_zero = sum_zero - list[i][out]->get_cnt_zeros();
    assert(!diff_cnt);
    assert(!diff_sum);
    assert(!diff_zero);
  }
#endif
} 

/**
   @brief Merging two structures.
   @param <obj> The object who gives away data.
   @remark 
   # Do neither release any memory, or allocates any, due to assumption that the sizes of both equals.
   # Validateion executed in the sub-call.
**/
void list_norm::merge_basis(list_norm &obj) {
#ifndef NDEBUG
  if(obj.get_max_input_value()>max_input_value) max_input_value = obj.get_max_input_value();
  //! Collection data pre to the operation:
  const bool print_verbouse = false; // Set to true if problems arises after code changes in the future.
  assert(list);
  uint cnt=0; float sum=0; uint cnt_zero=0;
  uint cnt_temp; float sum_temp; uint cnt_zero_temp;
  get_total_score(obj, taxon_length, cnt_temp, sum_temp, cnt_zero_temp, max_input_value);
  cnt += cnt_temp, sum += sum_temp, cnt_zero += cnt_zero_temp;
  if(print_verbouse)  printf("-\targument has: cnt(%u), sum(%f), zeros(%u) at line %d in file %s\n", cnt_temp, sum_temp, cnt_zero, __LINE__, __FILE__); 
  get_total_score(*this, taxon_length, cnt_temp, sum_temp, cnt_zero_temp, max_input_value);
  if(print_verbouse)  printf("-\t'this' has: cnt(%u), sum(%f), zeros(%u) at line %d in file %s\n", cnt_temp, sum_temp, cnt_zero, __LINE__, __FILE__); 
  cnt += cnt_temp, sum += sum_temp, cnt_zero += cnt_zero_temp;
  if(print_verbouse)  printf("-- Before:     In total cnt(%u), sum(%f), cnt_zero(%u) in list_norm\n", cnt, sum, cnt_zero);  
#endif
  for(uint i=0; i<(uint)taxon_length; i++) {
    assert(list[i]);
    if(obj.has_data(i)) {
      for(uint out=0; out<(uint)taxon_length; out++) {
	merge_basis(obj, i, out);
      }
    }
  }
#ifndef NDEBUG
  //! The resulting test of this method:
  //! Note: For big floating-point-numbers rounding always is a problem, as rounding of decimals may enough times cause a large difference:
  get_total_score(*this, taxon_length, cnt_temp, sum_temp, cnt_zero_temp, max_input_value);
  if(!(cnt == cnt_temp) //|| (!is_not_different(sum_temp, sum)) 
     || !(cnt_zero_temp == cnt_zero)) {
    
    fprintf(stderr, "!!\tError in norm-value-merging at line %d in file %s. Please contact the developer at oekseth@gmail.com\n", __LINE__, __FILE__);
    printf("-- Before:     In total cnt(%u), sum(%f), cnt_zero(%u) in list_norm\n", cnt, sum, cnt_zero); 
    printf("-- After:      In total cnt_temp(%u), sum_temp(%f), cnt_zero_temp(%u) in list_norm\n", cnt_temp, sum_temp, cnt_zero_temp); 
    printf("-- Difference: cnt(%d), sum(%f), cnt_zero(%d) in list_norm\n", cnt_temp - cnt, sum_temp - sum, cnt_zero_temp - cnt_zero); 
    assert(cnt == cnt_temp);
    assert((int)sum_temp == (int)sum);
    assert(cnt_zero_temp == cnt_zero); 
  }
#endif
  
}


//! Method below used to avoid the difficulties with round floating-pointer-numbers.
bool list_norm::is_not_different(float s1, float s2) {
  float diff = s1-s2;
  if(diff> 2 || diff < -2) return false;
  else return true;
}

//! Initiates the list using the already set data
void list_norm::init_list() {
  assert(!list); // verifies that we do not overwrite any data.
  assert(taxon_length); // Verifies its set
  list = new norm_t**[taxon_length];
  assert(list != NULL);
  if(list) {
    for(uint i = 0; i< (uint)taxon_length; i++) {
      list[i] = new norm_t*[taxon_length];
      if(list[i]) {
	for(uint out = 0; out < (uint)taxon_length; out++) list[i][out] = NULL;
      } else {
	fprintf(stderr, "!!\tWas not able allocating data due to memory constraints. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", __LINE__, __FILE__, __FUNCTION__);
	exit(2);
      }
      assert(list[i]);
    }
  } else {
    fprintf(stderr, "!!\tWas not able allocating data due to memory constraints. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", __LINE__, __FILE__, __FUNCTION__);
    exit(2);
  }
}

/**
   @brief Merging several structures into one.
   @param <index_root> The part of the index to 'return, ie, to give a way and set to NULL, thereby not removing it when doing the 'cleaning' of the list_norm object.
   @param <n_threads> The total number of squared **norm structures.
   @param <input> The object who gives away data
   @param <max_input_value> The maximal sim score found in the blast file.
   @return The result object.
   @remark Do not deallocates memory, but transform those with value of '0' according to the rules given.
**/
void list_norm::merge_basis_zero(uint index_root, uint n_threads, list_norm **&input, float max_input_value) {
#ifndef NDEBUG
  assert(input);
  assert(list);
  uint cnt=0; float sum=0; uint cnt_zero=0;
  for(uint i = 0; i< (uint)n_threads; i++) {
    uint cnt_temp; float sum_temp; uint cnt_zero_temp;
    get_total_score(*input[i], taxon_length, cnt_temp, sum_temp, cnt_zero_temp, max_input_value);
    cnt += cnt_temp, sum += sum_temp, cnt_zero += cnt_zero_temp;
  }
#endif
  const float zero_add = max_input_value + 1;
  // Seperate foor-loop in order for the compiler to ease the roll-out of it.
  for(uint i = 0; i< (uint)n_threads; i++) {
    input[i]->changeZeroToValue(zero_add); // norm[in][out].changeZeroToValue(zero_add); // updates.
  }       
  for(uint in = 0; in<(uint)taxon_length; in++) {
    for(uint out = 0; out<(uint)taxon_length; out++) {
      for(uint i = 0; i< (uint)n_threads; i++) {
	if(i != index_root) merge_basis(*input[i], in, out);
      }
    }
  }
#ifndef NDEBUG
  uint cnt_temp=0; float sum_temp=0; uint cnt_zero_temp=0;
  get_total_score(*this, taxon_length, cnt_temp, sum_temp, cnt_zero_temp, max_input_value);
  
  assert(cnt_temp == cnt);
  assert((int)sum_temp == (int)sum);
  assert(cnt_zero_temp == cnt_zero);
#endif
}

//! Converts the '0's to a mutliplication of the inputs:
void list_norm::changeZeroToValue(float max_input_pluss_one) {
  if(list) {
    for(uint i=0; i<(uint)taxon_length; i++) {
      if(list[i]) {
	for(uint out=0; out<(uint)taxon_length; out++) {
	  if(list[i][out]) list[i][out]->changeZeroToValue(max_input_pluss_one);
	}
      }
    }
  }
}


//! Prints the data in form of the labels given in the input file:
void list_norm::print_labels(struct taxa *listTaxa) {
  if(list) {
    for(uint i=0; i<(uint)taxon_length; i++) {
      if(list[i]) {
	for(uint out=0; out<(uint)taxon_length; out++) {
	  if(list[i][out]) {
	    list[i][out]->print_label(listTaxa, i, out);
	    list[i][out]->print_label(listTaxa, out, i);
	  }
	}
      }
    }
  }
}

//! Prints the details of the normalisation procedure, given the argument list.
void list_norm::print_normalization_basis(float **arrAvgNorm) {
  assert(arrAvgNorm);
  if(DEBUG || PRINT_NORMALIZATION_BASIS) { 
    printf("\n# The basis for the normalization procedure (defined in file 'norm_t.cxx' at line %d) is:\n", __LINE__);
    printf("# The resulting arrAvgNorm (i.e. what the similarity-score(i,j) is to be divided upon) is:\n");
    if(arrAvgNorm) {
      for(uint i = 0; i< (uint)taxon_length; i++) {
	if(arrAvgNorm[i]) {
	  for(uint out = 0; out< (uint)taxon_length; out++) {
	    printf("\tarrAvgNorm[%i][%i] = %f\n", i, out, arrAvgNorm[i][out]);
	  }	  
	} else printf("(not set for [%d])\n", i);
      }
    } else printf("(not set)\n");
  }
}


//! Converts the zero-value for the given index, setting the paramters to the updated values.
void list_norm::changeZeroToValue(uint taxon_in, uint taxon_out, float max_input_pluss_one, mem_loc &cnt, float &sum) {
  assert(list);
  assert(list[taxon_in]);
  if(list[taxon_in][taxon_out]) {
    list[taxon_in][taxon_out]->changeZeroToValue(max_input_pluss_one);
    cnt = list[taxon_in][taxon_out]->get_cnt();
    sum = list[taxon_in][taxon_out]->get_sum(max_input_value);
  } else {cnt = 0, sum = 0;}
}

//!  Builds the basis for the normative proecdure
float** list_norm::build_basis() { 
  assert(list);
  assert(taxon_length);
  float **arrAvgNorm = (float**)malloc(sizeof(float*)*taxon_length);
  test_condition_and_if_not_abort(arrAvgNorm, __LINE__, __FILE__, __FUNCTION__);
  for(uint i = 0; i< (uint)taxon_length; i++) {
    arrAvgNorm[i] = (float*)malloc(sizeof(float)*taxon_length);  
    test_condition_and_if_not_abort(arrAvgNorm[i], __LINE__, __FILE__, __FUNCTION__);
  }
  for(uint i = 0; i< (uint)taxon_length; i++) {
    if(list[i] && list[i][i]) arrAvgNorm[i][i] = list[i][i]->get_basis(max_input_value); // (sum of sim-scores)/(total number of protein-pairs)
    else arrAvgNorm[i][i] = 1.0;
  }

  const float zero_add = max_input_value + 1;
  // Iterates over the leftmost (inner) taxon in the pair.
  // (As specified in the documentation, a unique number represents a taxon).
  for(uint i = 0; i< (uint)taxon_length; i++) { 
    assert(list[i]);
    // Iterates over the rightmost (outer) taxon in the pair.  (As specified in the documentation, a unique number represents a taxon).
    // The 'list' contains information about the similarity scores for a given collection of protein pairs.
    for(uint out = i+1; out< (uint)taxon_length; out++) {
      assert(list[out]);
      // Changes possible values from zero to real value:
      mem_loc cnt_first =0;  float sum_first = 0; changeZeroToValue(i, out, zero_add, cnt_first, sum_first);
      mem_loc cnt_second=0; float sum_second = 0; changeZeroToValue(out, i, zero_add, cnt_second, sum_second);
      if(cnt_first || cnt_second) { 	// To avoid division by zero.
	arrAvgNorm[out][i] // Calculates the averaged sum:
	  =  (sum_first + sum_second) // The 'sum' corresponds to all of the similarity scores summed together
	  /  (cnt_first + cnt_second) // The 'cnt' corresponds to number of proteins used as a basis for the summing.
	  ;
	arrAvgNorm[i][out] = arrAvgNorm[out][i]; // Makes them reciprocal.
      } else {
	arrAvgNorm[i][out] = 1.0; // No data set
	arrAvgNorm[out][i] = 1.0; // No data set
      }
    }
  }
  print_normalization_basis(arrAvgNorm);
  log_builder::test_memory_condition_and_if_not_abort(arrAvgNorm!=NULL, __LINE__, __FILE__, __FUNCTION__);
  //    printf("---\t\t arrAvgNorm_is_Set(%d) at line %d in file %s\n", arrAvgNorm!= NULL, __LINE__, __FILE__);
  return arrAvgNorm;
}

//! Changes this new type of matrix to the old format for backwards compability.
norm_t **list_norm::change_to_old_matrix() {
  norm **obj = norm::init_arrNorm(taxon_length, PRINT_NORMALIZATION_BASIS, DEBUG_NORM, max_input_value);
  for(uint i = 0; i< (uint)taxon_length; i++) { 
    assert(list[i]);
    for(uint out = 0; out< (uint)taxon_length; out++) {
      if(list[i][out]) obj[i][out].merge(*list[i][out]);
    }
  }
  return obj;
}

#ifdef USE_MPI 
#include "mpi_transfer_list_norm.h"
/**
   @brief Sends- and receives the list of norm_t* objects accross nodes.
**/
/**
   @brief Sends- and receives the list of norm_t* objects accross nodes.
**/
void list_norm::mpi_make_data_consistent_accross_nodes(MPI_Comm comm, int myrank, int number_of_nodes) {
  mpi_transfer_list_norm transfer = mpi_transfer_list_norm(taxon_length);
  if(myrank != 0) {
    if(taxon_length) { // Only sends the data if it's set:
      transfer.build_list(max_input_value, taxon_length, list);
      transfer.send_data_to_root(comm, myrank, number_of_nodes);
    }
    //! Fills- and receives the updated list:
    transfer.broadcast_get_and_merge_data(this, comm, myrank, number_of_nodes);
  } else { // the root:
    transfer.get_and_merge_data(this, comm, myrank, number_of_nodes);
    //! In order to include data of 'this', the object is first cleadn, and then rebuilt:
    transfer.deallocate_internal_list();
    transfer.build_list(max_input_value, taxon_length, list);
    transfer.broadcast_send_merged_data(comm, myrank, number_of_nodes);
  }
   transfer.assert_result(this, comm, myrank, number_of_nodes);
   transfer.free_mem();
}

#endif

//! Allocates- and intiates memory for the 'norm_t' structure
list_norm **list_norm::init_list(uint list_size, uint taxon_length, float max_input_value, bool PRINT_NORMALIXATION_BASIS, bool DEBUG_NORM) {
  list_norm **obj = new list_norm*[list_size];
  test_condition_and_if_not_abort(obj, __LINE__, __FILE__, __FUNCTION__);
  for(uint i =0;i < list_size;i++) {
    obj[i] = new list_norm(taxon_length, max_input_value, PRINT_NORMALIXATION_BASIS, DEBUG_NORM);
    test_condition_and_if_not_abort(obj[i], __LINE__, __FILE__, __FUNCTION__);
  }
  return obj;
}

//! Allocates- and intiates memory for the 'norm_t' structure
list_norm **list_norm::init_list(uint list_size) {
  list_norm **object = new list_norm*[list_size];
  log_builder::test_memory_condition_and_if_not_abort(object!=NULL, __LINE__, __FILE__, __FUNCTION__);
  object[0] = NULL;
  test_condition_and_if_not_abort(object, __LINE__, __FILE__, __FUNCTION__);
  for(uint i =0;i < list_size;i++) {
    object[i] = NULL;
  }
  return object;
}

//! Deallocates memory for the list of list_norm objects.
void list_norm::close(list_norm **&obj, uint list_size) {
  if(obj) {
    for(uint i =0;i < list_size;i++) {
      if(obj[i]) {
	obj[i]->free_mem();
	delete obj[i]; obj[i] = NULL;
      }
    }
    delete [] obj; obj = NULL;
  }
}

//! Deallocates memory for the list of list_norm objects.
void list_norm::close(list_norm *&obj) {
  //    printf("arrNorm(1==%d) at line %d in file %s\n", (NULL!=obj), __LINE__, __FILE__);
  if(obj) {
    obj->free_mem();
    delete obj; obj = NULL;
  }
}

//! Deallocates memory for the object given.
void list_norm::free_mem() {
  if(list) {
    if(taxon_length) {
      //	printf("closes a list of length(%d) at line %d in file %s\n", taxon_length, __LINE__, __FILE__);
      for(uint i = 0; i< (uint)taxon_length; i++) {
	if(list[i]) {
	  for(uint out = 0; out< (uint)taxon_length; out++) {
	    if(list[i][out]) delete list[i][out]; list[i][out] = NULL;
	  }
	  delete [] list[i]; list[i] = NULL;
	}
      }
    }
    delete [] list; list = NULL;
  }
}

//! The constructor.
list_norm::list_norm() : list(NULL), taxon_length(0), max_input_value(0.0), PRINT_NORMALIZATION_BASIS(false), DEBUG_NORM(false)
{;};

//! The constructor.
list_norm::list_norm(uint t_length, float m_value, bool norm_basis, bool d_norm) :
  list(NULL), taxon_length(t_length), max_input_value(m_value), PRINT_NORMALIZATION_BASIS(norm_basis), DEBUG_NORM(d_norm)
{
  if(t_length) init_list();
};

