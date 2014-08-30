#include "taxon_data.h"
//! Prints data about this object.
void taxon_data::print_data(uint index) {
  if(index < prots_reserved) printf("%u\t%s\n", index, buffer[index]);
  else printf("%u\t%s\n", index, ".nil.");
}

//! Returns the name given the param.
char *taxon_data::get_name(uint index) {
  if(index <= prots_used) return buffer[index];
  else return NULL;
}

//! Returns true if input name is the same as the taxon-name of this object.
bool taxon_data::is_equal(char *name) {
  return (0== strcmp(taxon_name, name)); 
}

//! Sets the taxon name
void taxon_data::setTaxonName(char *taxon) {
  if(taxon != NULL) {
    const loint size = strlen(taxon);
    assert(size);
    taxon_name = new char[size+1]; taxon_name[size] = '\0';
    log_builder::test_memory_condition_and_if_not_abort(taxon_name!=NULL, __LINE__, __FILE__, __FUNCTION__);
    strncpy(taxon_name, taxon, size);
  }
}

//! Sets the param with value of this object.
void taxon_data::get_taxon_name(char *&name, loint &size_arg) {
  const loint size_this = strlen(taxon_name);
  if(!name || (size_this > (size_arg-1))) {
    if(name) {delete [] name; }
    name = new char[size_this+1];
    memset(name, '\0', size_this);
    log_builder::test_memory_condition_and_if_not_abort(name!=NULL, __LINE__, __FILE__, __FUNCTION__);
    size_arg = size_this+1;
  } 
  name[size_this]='\0'; // Sets an end
  strncpy(name, taxon_name, size_this);
}


//! Returns the last inserted protein: If empty, returns NULL
char *taxon_data::getLastInsertedProtein(){
  if(prots_used >0) return buffer[prots_used-1];
  return NULL;
}


//! prints the buffer
void taxon_data::printf_buffer() {
  printf("\t\t\tThe buffer for taxon %s contains the following proteins:\n", taxon_name);
  for(uint i = 0;i<prots_used;i++) {
    printf("\t\t\t[%d]\t%s\n", i, buffer[i]);
  }
}

/**
   @brief The 'main method' of this structure: 
   @param<name> Inserts it into the the list.
   @remark If buffer is full, enalarges the list.
   @return The pointer to the name given as input argument.
   @date 01.12.2010
**/
char* taxon_data::insert_protein(char *name) {
  if(name && ((prev_protein_name == NULL) || (0!= strcmp(prev_protein_name, name)))) { // there is a new protein
    if(prots_used == prots_reserved) { // Resize
      enlarge_buffer(10);
    }
    assert(buffer);
    const loint size = strlen(name);
    if(size) {
      if(!buffer[prots_used]) {buffer[prots_used] = new char[size+1]; buffer[prots_used][size] = '\0';}
      assert(buffer[prots_used]);    // Asserts that enough memory was allocated
      strncpy(buffer[prots_used], name, size); 
      prev_protein_name = buffer[prots_used];
      prots_used++;
      return prev_protein_name;
    } else {
      fprintf(stderr, "!!\tDiscarded name(%s), as the length of it is unkown. This error was generated at line %d in file %s. Please contact oesketh@gmail.com if this message is seen:\n", name, __LINE__, __FILE__);
    }
  }
  return NULL;
}


//! Inserts the label at specific index
void taxon_data::insert_label(char *name, uint index) {
  if(name) {
    if(index > prots_used) prots_used = index;    
    if(index >= prots_reserved) { // Resize
      enlarge_buffer(0);
    }
    const loint size = strlen(name);
    if(size) {
      assert(buffer);
      if(!buffer[index]) {buffer[index] = new char[size+1]; buffer[index][size] = '\0';}
      assert(buffer[index]);
      strncpy(buffer[index], name, strlen(name));
      prev_protein_name = buffer[index];
    } else {
      fprintf(stderr, "!!\tDiscarded name(%s), as the length of it is unkown. This error was generated at line %d in file %s. Please contact oesketh@gmail.com if this message is seen:\n", name, __LINE__, __FILE__);
    }
  }

}

/**
   @brief Deletes the buffer holding the indexes
   @deprecated For backward compability,
*/
void taxon_data::delete_buffer() {
  if(buffer != NULL) {
    for(uint i = 0; i < prots_reserved; i++) {
      delete [] buffer[i], buffer[i] = NULL;
      //free(buffer[i]);
    } free(buffer); buffer = NULL; prots_reserved = 0, prots_used=0;
  }
  if(taxon_name!= NULL) {delete [] taxon_name;  taxon_name = NULL;}
}

//! Merges the argument 'taxon_data' with 'this'
void taxon_data::mergeData(taxon_data data, bool overlap) {
  uint offset = 0;
  if(overlap) offset = 1;
  else if (prots_used > 0) {
    if(compare_strings(buffer[prots_used-1], data.buffer[0])) {
      offset=1;
    }
  }
  enlarge_buffer(data.prots_used);
  assert(buffer);
#ifndef NDEBUG
  //! Asserts that pairs given are not found in this, ie, no overlapping pairs
  //  seen from the arguments perspective (in order to mkae the test simpler as
  //  it else would have been written)
  for(uint out = offset; out< data.prots_used; out++) {
    for(uint i = 0; i < prots_used; i++) {
      assert(0!=strcmp(data.buffer[out], buffer[i]));
    }
  }
#endif
  for(uint i = offset; i< data.prots_used; i++) {
    const loint size = strlen(data.buffer[i]);
    if(size) {
      if(!buffer[prots_used]) {buffer[prots_used] = new char[size+1]; buffer[prots_used][size] = '\0';}
      assert(buffer[prots_used]);    // Asserts that enough memory was allocated
      //      assert(data.buffer[i]);
      strcpy(buffer[prots_used], data.buffer[i]);
      prots_used++;
    }
  }
}

/**
   Enlarges the buffer if the data to add is larger than the free space available
   @param <min_size> The minimum size the number of elements the buffer must consist of.
**/
void taxon_data::enlarge_buffer(loint min_size)  {
  loint new_size = prots_used + min_size;  
  if(!prots_reserved) assert(!buffer);
  if(new_size && (new_size >= prots_reserved)) {
    //    fprintf(stderr, "prots_used(%d), prots_reserved(%d), new_size(%d)\tenlarge_buffer(%lld)\t(In taxon_data at lines %d)\n", (int)prots_used, (int)prots_reserved, (int)new_size, (loint)min_size, __LINE__);
    //    fprintf(stderr, "-\tin taxon_data new_size>prots_reserved\n");
    // TODO: Improve cheduling below using empirical tests measuring relative time cost versus memory consumption.
    new_size += 2;
    //    new_size *=2;
    //    if(prots_reserved < 10) new_size = prots_reserved + 2; // To make ig grow lesser in the start.
    //    else if(prots_reserved < 100) new_size = prots_reserved + 10; // To make ig grow lesser in the start.
    //    printf("new_size=%u in taxon_data.cxx\n", new_size);
    char **buffer_new = (char**)malloc(sizeof(char*)*new_size);
    assert(buffer_new);
    log_builder::test_memory_condition_and_if_not_abort(buffer_new!=NULL, __LINE__, __FILE__, __FUNCTION__);
    if(buffer) { // Buffer has data
      for(loint i = 0; i<prots_reserved; i++) {
	//	assert(buffer[i]);
	buffer_new[i] = buffer[i]; // Copies the memory references
      }
      free(buffer); buffer = NULL;
    } else assert(prots_reserved <= 0); // Buffer is empty.
    for(loint i = prots_reserved; i< new_size; i++) {
      buffer_new[i] = NULL;
    }
    assert(!buffer);
    buffer = buffer_new;
    prots_reserved = new_size;
  } 
}

//! @return true if the identifier is the same as the argument (given as input).
bool taxon_data::is_taxon(char *arr_taxon_name) {
  return compare_strings(taxon_name, arr_taxon_name);
}

//! Prints the given string as-is:
void print_string(char *start, int length) {
  char *end = start + length;
  while(start < end) {putchar(*start), start++;}
}
#ifdef USE_MPI

/**
   @brief Sets the label of proteins
   @param <protein_length> The number of proteins to use as input.
   @param <listOfProteins> The 1d-list holding the protein names.
   @param <listOfLengths>  The list holding the length of each protein.
   @return The number of chars used for the names inserted.
   @remarks Important is to use the correct offsets when calling these functions!
**/
uint taxon_data::set_labels_using_1d_object(uint protein_length, char *listOfProteins, uint *listOfLengths) {
  assert(!prots_used);
  enlarge_buffer(protein_length);
  uint listOfProteins_cnt = 0;
  for(uint i = 0; i < protein_length; i++) {
    const uint size = listOfLengths[i];
    if(!buffer[i]) {buffer[i] = new char[size+1]; buffer[i][size] = '\0';}
    assert(buffer[i]);    // Asserts that enough memory was allocated
    char *name = listOfProteins + listOfProteins_cnt;
    strncpy(buffer[i], name, size); 
    prev_protein_name = buffer[i];
    listOfProteins_cnt += size;
  }
  prots_used = protein_length;
  return listOfProteins_cnt;
}

/**
   @brief Compares in order to verify the equalness.
   @param <obj> The object to compare with this object.
   @param <verbouse> Exaplains all the differences found.
   @return true if they are equal.
**/
bool taxon_data::compare_object(taxon_data &obj, bool verbouse) {
  bool is_equal = true;
  //! Compare the names
  if(0!=strcmp(taxon_name, obj.get_taxon_name())) {
    is_equal = false;
    if(verbouse) fprintf(stderr, "!!\ttaxon_name(%s) != obj.taxon_name(%s) at line %d in file %s\n", taxon_name, obj.get_taxon_name(), __LINE__, __FILE__);
  }
  if(is_equal || verbouse) { // If it's already false, only do it if we are interested in the complete set of differences.
    if(buffer) {
      for(uint i = 0; i < prots_used; i++) {
	if(0!=strcmp(buffer[i], obj.get_name(i))) {
	  is_equal = false;
	  if(verbouse) fprintf(stderr, "!!\tprotein(%s) != obj.protein(%s) at line %d in file %s\n", buffer[i], obj.get_name(i), __LINE__, __FILE__);
	}
      }
    } else {
      if(obj.get_length()) { // Argument has data, but not this:
	is_equal = false;
	if(verbouse) fprintf(stderr, "!!\tThis object is empty, while argument has length(%d) at line %d in file %s\n", obj.get_length(), __LINE__, __FILE__);	  
      }
    }
  }
  return is_equal;
}


/**
   @brief Inserts the taxon-name into the list.
   @remarks In order to easyly be sendt with the mpi-procedure, each name is ended with an '\0' symbol.
**/
void taxon_data::insert_the_name_into_argument(char *&taxonNamesArr, uint &taxonNamesArr_size, uint &taxonNamesArr_used) {
  //! Ensures the length of the data corresponds to our needs:
  if(!taxonNamesArr) {
    assert(!taxonNamesArr_used);
    assert(!taxonNamesArr_size);
    const uint tmp_size = 40*100; // A random number implying 100 taxa in total:
    char *tmp = new char[tmp_size];
    memset(tmp, '\0', tmp_size);
    taxonNamesArr_size = tmp_size, taxonNamesArr = tmp;
  } else {
    const uint required_length = taxonNamesArr_used + strlen(taxon_name); // An estimate that most likely will be 1.5 above the actual consumption.
    if(required_length > taxonNamesArr_size) {
      assert(taxonNamesArr_size);
      const uint tmp_size = 2*required_length;
      char *tmp = new char[tmp_size];
      memset(tmp+taxonNamesArr_used, '\0', tmp_size-taxonNamesArr_used);
      memcpy(tmp, taxonNamesArr, taxonNamesArr_used);
#ifndef NDEBUG
      //! Validates the enlarging of the lists:
      for(uint i = 0; i < taxonNamesArr_used; i++) {
	assert(taxonNamesArr[i] == tmp[i]);
      }
#endif
      delete [] taxonNamesArr; taxonNamesArr_size = tmp_size, taxonNamesArr = tmp;
    }
  }
  memcpy(taxonNamesArr+taxonNamesArr_used, taxon_name, strlen(taxon_name));
#ifndef NDEBUG
  if(0!=strcmp(taxon_name, taxonNamesArr+taxonNamesArr_used)) {
    fprintf(stderr, "!!\t\ttaxon_name(%s) != updated(%s) at line %d in file %s\n", taxon_name, taxonNamesArr+taxonNamesArr_used, __LINE__, __FILE__);
    assert(0==strcmp(taxon_name, taxonNamesArr+taxonNamesArr_used));
  }
#endif
  taxonNamesArr_used += strlen(taxon_name);
  taxonNamesArr[taxonNamesArr_used] = '\0'; // Appends a trailing char:
  taxonNamesArr_used++; // In order not to overwrite the '\0'-symbol.
}
/**
   @brief Asserts that the the list is correctly transformed.
   @remarks Only active when DBEBUG macro varaible is not defined.
**/
void taxon_data::assert_1d_list(uint bufferLength_startpos, uint &start_position_in_buffer, char *bufferTmp, uint *bufferLength) {
#ifndef NDEBUG
  //! Validates the set of lists built:
  uint index_index = bufferLength_startpos;
  uint offset = start_position_in_buffer;
  for(uint i = 0; i < prots_used; i++, index_index++) {
    //	if(buffer[i])
    if(strlen(buffer[i]) != bufferLength[index_index]) {
      printf("!!\t%s\tstring-length-error:\t%d VS %u (for index %u VS %u)\n", buffer[i], (int)strlen(buffer[i]), bufferLength[index_index], i, index_index);
         assert(strlen(buffer[i]) == bufferLength[index_index]); 
    }
    if(0!=strncmp(buffer[i], bufferTmp + offset, strlen(buffer[i]))) {
      printf("!!string-transformation-error:\t%s VS ", buffer[i]);
      print_string(bufferTmp+offset, strlen(buffer[i]));
      printf("\n");      
      assert(0==strncmp(buffer[i], bufferTmp + offset, strlen(buffer[i]) ));
    }
    offset += bufferLength[index_index];
    start_position_in_buffer += bufferLength[index_index];
  }
#endif
}

//! Prints the 1d-buffer:
void print_1d_buffer(char *bufferTmp, uint *bufferLength, uint bufferLength_size) {
  if(bufferLength_size) {
    assert(bufferLength);
    char *start = NULL;
    if(bufferTmp) start = bufferTmp; // + bufferTmp_char_used;
    printf("\t- Prints the buffer (of size=%u) at line %d in file %s\n", bufferLength_size, __LINE__, __FILE__); 
    for(uint i = 0; i < bufferLength_size; i++) {
      if(bufferLength && bufferLength[i]) {
	//	  printf("[%u] (length=%u)\n", i, bufferLength[i]);	
	printf("[%u] ", i);
	//printf("[%u] (length=%u)", i, bufferLength[i]);	
	if(start) {
	  char *end = start + bufferLength[i];
	  assert(end != start);
	  while(start != end) {
	    //	      printf("'%c'", *start);
	    putchar(*start);
	    start++;
	  }
	} else printf(" (undef) ");
	printf("\n");
      }
    } 
  } else      printf("(empty)\t- Prints the buffer (of size=%u) at line %d in file %s\n", bufferLength_size, __LINE__, __FILE__); 
}
//! Transform "**buffer" into "*bufferTmp" and "int *bufferLength". Is validated internally. Used for MPI
void taxon_data::transform_protein_label_buffer_into_1d_list(char *&bufferTmp, uint &bufferTmp_char_used, uint &bufferTmp_size,uint *&bufferLength,uint &bufferLength_size, uint bufferLength_startpos){
#ifndef NDEBUG
  const uint start_position_in_buffer = bufferTmp_char_used;
#endif
  if(buffer && prots_used) {
    //! Ensures that the index-list is long enough:
    if(!bufferLength) {
      assert(!bufferLength_startpos);
      assert(!bufferLength_size);
      bufferLength = new uint[prots_used];
      log_builder::test_memory_condition_and_if_not_abort(bufferLength!=NULL, __LINE__, __FILE__, __FUNCTION__);
      memset(bufferLength, 0, sizeof(uint)*bufferLength_size);
      bufferLength_size = prots_used;
    } else {
      assert(bufferLength_size>0);
      const uint required_length = (bufferLength_startpos+prots_used);
      if(required_length > bufferLength_size) {
	//! Resizes the list:
	const uint tmp_size = 2*(bufferLength_startpos+prots_used);
	uint *tmp = new uint[tmp_size];
	log_builder::test_memory_condition_and_if_not_abort(tmp!=NULL, __LINE__, __FILE__, __FUNCTION__);
	memset(tmp + bufferLength_startpos, 0, sizeof(uint)*(tmp_size-bufferLength_startpos));
	memcpy(tmp, bufferLength, sizeof(uint)*bufferLength_startpos);
#ifndef NDEBUG
	//! Validates the enlarging of the lists:
	for(uint i = 0; i < bufferLength_startpos; i++) {
	//	for(uint i = 0; i < bufferLength_size; i++) {
	  assert(bufferLength[i] == tmp[i]);
	}
#endif
	// Update the variable:
	delete [] bufferLength; bufferLength = tmp; bufferLength_size = tmp_size;
      }
    }
    //! Ensures that the buffer-list is long enough:
    if(!bufferTmp) {
      assert(!bufferTmp_char_used);
      assert(!bufferTmp_size);
      const uint tmp_size = 2*40*prots_used; // An estimate that most likely will be 1.5 above the actual consumption.
      char *tmp = new char[tmp_size];
      memset(tmp, '\0', bufferTmp_size);
      bufferTmp = tmp; bufferTmp_size = tmp_size;
      log_builder::test_memory_condition_and_if_not_abort(bufferTmp!=NULL, __LINE__, __FILE__, __FUNCTION__);
    } else {
      const uint required_length = bufferTmp_char_used + 40*prots_used; // An estimate that most likely will be 1.5 above the actual consumption.
      if(required_length > bufferTmp_size) {
	assert(bufferTmp_size);
	const uint tmp_size = 2*required_length;
	char *tmp = new char[tmp_size];
	log_builder::test_memory_condition_and_if_not_abort(tmp!=NULL, __LINE__, __FILE__, __FUNCTION__);
	memset(tmp + bufferTmp_char_used, '\0', tmp_size-bufferTmp_char_used);
	memcpy(tmp, bufferTmp, bufferTmp_char_used);
	//! Validates the enlarging of the lists:
	for(uint i = 0; i < bufferLength_size; i++) {
	  assert(bufferTmp[i] == tmp[i]);
	}
	// Update the variable:
	delete [] bufferTmp; bufferTmp = tmp; bufferTmp_size = tmp_size;
      }
    }
    for(uint i = 0; i < prots_used; i++) {
      if(buffer[i] and buffer[i][0] != '\0') {
	const uint length = strlen(buffer[i]);
	if((bufferTmp_char_used + length) >= bufferTmp_size) {
	  // Enlarges the list:
	  const int new_size = bufferTmp_size*2;
	  char *buffer_new = new char[new_size];
	  memset(buffer_new + bufferTmp_size, '\0', bufferTmp_size); // Initates the part of the buffer not to be copied into.
	  memcpy(buffer_new, bufferTmp, bufferTmp_char_used);
	  delete [] bufferTmp;
	  bufferTmp = buffer_new; bufferTmp_size = new_size;
	}
	memcpy(bufferTmp+bufferTmp_char_used, buffer[i], length);
	bufferTmp_char_used += length; 
	bufferLength[bufferLength_startpos + i] = length; // Sets for the given protein the index.	  
      }
    }
#ifndef NDEBUG
    uint start_position_in_buffer_temp = start_position_in_buffer;
    assert_1d_list(bufferLength_startpos, start_position_in_buffer_temp, bufferTmp, bufferLength);
#endif
  }
}
#endif
/**
   @param <size> The size of this object not initialized.
   @return A list of this object
*/
taxon_data *taxon_data::init_empty(uint size) {
  return taxon_data::init(size);
  //  return (taxon_data*)malloc(sizeof(taxon_data)*size);
}

/**
   @param <size> The size of this object initilazed.
   @return A list of this object
*/
taxon_data *taxon_data::init(uint size) {
  if(size) {
    taxon_data *obj = new taxon_data[size]();
    log_builder::test_memory_condition_and_if_not_abort(obj!=NULL, __LINE__, __FILE__, __FUNCTION__);
    return obj;
  } 
  else return NULL;
}

//! The constructor.
taxon_data::taxon_data()  :
  prots_reserved(0), prots_used(0), buffer(NULL) {
  taxon_name = NULL;
  prev_protein_name = NULL;
}

//! @brief The constructor. @param <reserve_size> The length of the data to reserve memory for.
taxon_data::taxon_data(loint reserve_size) :  
  prots_reserved(0), prots_used(0), buffer(NULL)
{
  taxon_name = NULL;
  if(reserve_size) enlarge_buffer(reserve_size);
  prev_protein_name = NULL;
}

#ifdef assert_code
// The following code contains test function for this class:
//! Testing the merging of two data strctures
void taxon_data::assert_mergeData() {
  // Intializes the inner
  taxon_data list = taxon_data();  
  char *name_in[] = {"sara", "helmut", "jon", "lars"};
  //    uint name_in_size = 4;
  uint name_in_cnt = 0;
  list.insert_protein(name_in[name_in_cnt++]);
  list.insert_protein(name_in[name_in_cnt++]);

  // Intializes the outer:
  taxon_data other = taxon_data();
  char *name_out[] = {"sara", "helmut", "jon", "lars"};
  uint name_out_cnt = 0;
  other.insert_protein(name_out[name_out_cnt++]);
  other.insert_protein(name_out[name_out_cnt++]);
  list.mergeData(other, false);
  uint start = 0, end = name_in_cnt; 
  uint cnt_in = 0;     uint size_buffer = 0; char **buffer = list.getBuffer(size_buffer);
  for(uint i = start; i < end; i++) {
    assert((0 == strcmp(name_in[cnt_in++], buffer[i]))); // compares the argument witht the actual nam
  }
  start = name_in_cnt, end = size_buffer;
  uint cnt_out = 0; //if(show_positive_results) printf("length of 'name_out_cnt' = %u\n", name_out_cnt);
  for(uint i = start; i < end; i++) {
    assert((0 == strcmp(name_out[cnt_out++], buffer[i]))); // compares the argument witht the actual name
  }
  list.delete_buffer();   other.delete_buffer();
  taxon_data one = taxon_data();
  taxon_data two = taxon_data();
  char *name = "mironov";
  one.insert_protein(name);
  two.insert_protein(name);
  two.mergeData(one, true);
  size_buffer = 0; buffer = two.getBuffer(size_buffer);
  assert((0 == strcmp(name, buffer[0]))); // compares the argument witht the actual name
  assert((1 == size_buffer)); // compares the acutal length with the length computed
  one.free_mem(); two.free_mem();
}

//! Testing that the isnertin of proteins are correct:
//! - If not ewqual to the previous, insert
//! - If equal to the previous, disregard insertion, return NULL
void taxon_data::assert_insert_protein() {
  taxon_data list = taxon_data();
  //  test_cnt++; if(show_positive_results) printf("??\t Running test %u: Inserting protein\n", test_cnt);
  char *name = "hanne"; uint real_prots_used = 0;
  char *updated_name = list.insert_protein(name); real_prots_used++;
  assert((0 == strcmp(name, list.getFirstInsertedProtein()))); // compares the argument witht the first inserted name
  assert((0 == strcmp(name, list.getLastInsertedProtein()))); // compares the argument witht the last inserted name
  assert((0 == strcmp(updated_name, name))); // compares the argument witht the actual name
  assert( (real_prots_used == list.getProteinsUsed())); // verifies that the size is correct
  updated_name = list.insert_protein(name); 
  assert(updated_name == NULL); // compares the argument witht the actual name
  assert( (real_prots_used == list.getProteinsUsed())); // verifies that the size is correct
  name = "far";
  updated_name = list.insert_protein(name); real_prots_used++;
  assert((0 == strcmp(updated_name, name))); // compares the argument witht the actual name
  assert( (real_prots_used == list.getProteinsUsed())); // verifies that the size is correct
  uint size_buffer;
  char *name_new[] = {"hanne", "far"};
  char **buffer = list.getBuffer(size_buffer);
  for(uint i =0; i<size_buffer; i++)
    assert((0 == strcmp(name_new[i], buffer[i]))); // compares the argument witht the actual name
  list.free_mem();
  
}

//! Asserting that the taxon name is set correct, using the functions for inserting, and retrieval
void taxon_data::assert_setTaxonName(taxon_data &list) {
  char *taxon = "human";
  list.setTaxonName(taxon); // Inserts the name
  assert(list.is_taxon(taxon)); // The name should be inserted
}


#endif
//taxon_data::taxon_data()


//! The main test function for this class  
void taxon_data::assert_class(const bool print_info) {
  const static char *class_name = "taxon_data";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  taxon_data list = taxon_data();
  list.assert_setTaxonName(list);
  list.assert_mergeData();
  list.assert_insert_protein();
  list.free_mem();
  list.delete_buffer();
  list = taxon_data();
  list.setTaxonName("human___");
  list.setTaxonName("human___");
  list.setTaxonName("human___");
  list.setTaxonName("human___");
  list.delete_buffer();
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}


