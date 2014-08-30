#include "parse_read_blocks.h"


//! Asserts that the block ends with a new-line, and if not prints the debug-info
void parse_read_blocks::assert_newline_end(char *buffer_main, long int size) {
  if(!(buffer_main[size-1] == '\n')) {
    char *r = strrchr((buffer_main + size), '\n');
    if(r) {
      printf("(in parse_read_blocks.cxx at line %d\terror_size(%ld):\t", __LINE__, (buffer_main + size) - r);
      printf("[\n");
      while(r!= (buffer_main+size)) putchar(*r), r++;
      printf("]\n");
    }
    assert(buffer_main[size-1] == '\n');
  }
}

//! Gets the block using fseek:
char *parse_read_blocks::get_block(char *file_name, long int pos_in_file, long int size) {
  if(file_name && size) {
    char *buffer_main;
    allocate(buffer_main, size); // Uses offset of '100' to ensure the holw line comes in.
    read_file(file_name, pos_in_file, buffer_main, size);    
    return buffer_main;
  } else return NULL;
}


/**
   @brief Inserts the data inside the structure:
   @remarks Size of each block will (on average) be proportional to the size of the data used for input.
**/
void parse_read_blocks::insert(loint file_pos) {
  if(file_pos > 0 && (loint)file_pos > last_file_pos_in_mem) { // an input of '0' means that a seperator is not found
    if(size == position_in_arr) {
      const loint new_size = size*2;
      arr = (loint*)realloc(arr, sizeof(loint)*new_size);
      log_builder::test_memory_condition_and_if_not_abort(arr!=NULL, __LINE__, __FILE__, __FUNCTION__);
      memset(arr+size, 0, sizeof(loint)*size);
      size = new_size;
    }
    assert(arr);
    arr[position_in_arr] = file_pos;
    position_in_arr++;
  }
}


//! Allocates and initializes the read-buffer to be sent downward the pipe:
void parse_read_blocks::allocate(char *&buffer_main, const long int buffer_in_mem_size) {
  buffer_main = new char[buffer_in_mem_size+1];//(char*)malloc(sizeof(char)*(buffer_in_mem_size+1));
  assert(buffer_main);
  log_builder::test_memory_condition_and_if_not_abort(buffer_main!=NULL, __LINE__, __FILE__, __FUNCTION__);
  memset(buffer_main, '\0', buffer_in_mem_size+1);// Filling in the buffer
}

//! Allocates and initializes the read-buffer to be sent downward the pipe:
void parse_read_blocks::allocate(char *&buffer_main, char *&buffer_end, const long int buffer_in_mem_size) {
  buffer_main = new char[buffer_in_mem_size+1]; //(char*)malloc(sizeof(char)*(buffer_in_mem_size+1));//new char[buffer_in_mem_size + 1];
  log_builder::test_memory_condition_and_if_not_abort(buffer_main!=NULL, __LINE__, __FILE__, __FUNCTION__);
  assert(buffer_main);
  buffer_end = buffer_main + buffer_in_mem_size;
  memset(buffer_main, '\0', buffer_in_mem_size+1);// Filling in the buffer:
}

//! Reads the file:
void parse_read_blocks::read_file(char *file_name, long int pos_in_file, char *buffer_main, const long int buffer_in_mem_size) {
  assert(buffer_main);
  assert(file_name);
  assert(buffer_in_mem_size);
  FILE *f = fopen(file_name, "rb");
  if(f) {
    fprintf(stderr, "in parse_read_blocks(%d) --\tfseek called with value %llu\n", __LINE__, (loint)pos_in_file);
    fseek(f, pos_in_file, SEEK_SET);
    const long int read_length = fread(buffer_main, sizeof(char), buffer_in_mem_size, f);
    assert(read_length <= buffer_in_mem_size);
    fclose(f);
  } else {
    fprintf(stderr, "!!\tUnable to open file %s at line %d in file %s, found in method %s. Please contact oekseth@gmail.com it this message is seen at your terminal.\n", file_name, __LINE__, __FILE__, __FUNCTION__);
    assert(false);
  }
}


/**
   @brief Reads from the given file-pointer
   @param <file_input> The file pointer tor ead from.
   @param <buffer_main> The reserved char array to read the data into.
   @param <buffer_in_mem_size> The variable to get number of chars to read.
   @return true if the there are no more chars to read.
   @author Ole Kristian Ekseth (oekseth)
**/
bool parse_read_blocks::read_file(FILE *file_input, char *buffer_main, const long int buffer_in_mem_size) {
  assert(buffer_main);
  assert(file_input);
  assert(buffer_in_mem_size);
  assert(file_input);
  const long int fread_cnt = fread(buffer_main, sizeof(char), buffer_in_mem_size, file_input);
  if(fread_cnt == buffer_in_mem_size) return false;
  else return true;
}

//! Returns the section of the data using the parse-blocks as a basis:
string_section *parse_read_blocks::get_string_section_from_file(FILE *file_input) {
  long int buffer_in_mem_size = getNextReadLength();
  if(buffer_in_mem_size > 1) {
    if(!feof(file_input)) {
      char *buffer_main = NULL, *buffer_end = NULL; // Holding the grabbed string.
      allocate(buffer_main, buffer_end, buffer_in_mem_size+100); // Uses offset of '100' to ensure the holw line comes in.
      bool at_end_of_file = read_file(file_input, buffer_main, buffer_in_mem_size);
      if((buffer_main[buffer_in_mem_size-1] != '\n') && (buffer_main[buffer_in_mem_size]!='\n') && (buffer_main[buffer_in_mem_size]!='\0') ) {
	// TODO: Consider removing the lines below as the case hopefully is an unlikly one.
	fprintf(stderr, "If the following lines are found, generated from line %d in %s, please contact the developer at oekseth@gmail.com\n", __LINE__, __FILE__);
	fprintf(stderr, "An error at line %d in file parse_read_blocks.cxx:\n", __LINE__);
	if(true) fprintf(stderr, "--Found the char('%c') instead of a new line.\n", buffer_main[buffer_in_mem_size-1]);
	if(true) fprintf(stderr, "--Found the char('%c') instead of a new line.\n", buffer_main[buffer_in_mem_size]);
	if(false) {
	  printf("the block to send is:[\n");
	  for(loint i = 0; i < (loint)buffer_in_mem_size; i++) {
	    putchar(buffer_main[i]);
	  }
	  printf("\n]\n");
	}
	assert(false);
      }
      string_section_t *section = string_section::init(buffer_main, buffer_main + buffer_in_mem_size, 
						       /*prev_index=*/0, /*not important*/0,
						       /*end_of_file=*/at_end_of_file);
      return section;
    } 
  }
  return NULL;
}


//! Writing of the blocks in memory to stdout.
void parse_read_blocks::printBlocks() {
  loint total_length = 0;
  loint prev_pos = 0;
  assert(arr);
  for(loint i = 0; i < size; i++) {
    loint length = arr[i];
    if(i != 0) length-=prev_pos;
    if(arr[i] > 0) {
      printf("block %llu:\t in [%llu, %llu], with a size of %llu chars\n", i, prev_pos, arr[i], length);
      total_length += arr[i];
    }
    prev_pos =  arr[i];
  }
}

/**! Returns the length of the block, incrementing the poiner for the next read
   @Return: '0' if no more blocks to read
*/
long int parse_read_blocks::getNextReadLength() {
  assert(arr);
  if(read_postion_to_retrieve < size) {
    long int file_pos = last_file_pos_in_mem;
    if(read_postion_to_retrieve > 0) file_pos = arr[read_postion_to_retrieve-1];
    long int length = 0;
    do {
      length = arr[read_postion_to_retrieve] - file_pos; // new_pos - old_pos
      read_postion_to_retrieve++;
    } while((read_postion_to_retrieve < position_in_arr) && (length < disk_buffer_size)); 
    //printf("\t get string-block %u/%u, with size=%u, at parse_read_blocks:%d\n", (uint)read_postion_to_retrieve, (uint)size, (uint)length, __LINE__); // FIXME: remove.
    return length;
  }
  return 0;
}


/**
   @brief Generate (writes) a log file describing the details of the memory signature for this 'unit'
   @remarks 
   - Useful for analyzing- optimizing the memory allocation procedures.
   - In optimised mode does not produce the log file.
**/
void parse_read_blocks::log_generate_memory_allocation_overview(uint n_threads) {
#ifndef NDEBUG
  assert(n_threads);
#ifdef LOG_WRITE_MEMORY_ALLOCATIONS_DURING_BLASTP_PARSING_FOR_PARSE_BLOCKS
  assert(arr);
  struct tms tmsstart; clock_t clock_log;
  if((clock_log = times(&tmsstart)) == -1) // starting values
    fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
  char sep = '_';
  // The log file, and contents of it:
  if(sizeof(FUNC_LOG_FOLDER_NAME) > 0) sep = '/';
  char str[200]; memset(str, '\0', 200);
  sprintf(str, "%s%c%s.log.%u", FUNC_LOG_FOLDER_NAME, sep, "parse_read_blocks_memory", n_threads);
  FILE *f_log = fopen(str, "w");
  if(f_log) {
    fprintf(f_log, "Stored %.2fMB in memory, spreading it over %llu elements, with the specifi divsions for each of the %u threads given as, usinig and parse_read_blocks object:\n",
	    (float)arr[position_in_arr]/(float)(1024*1024), position_in_arr+1, n_threads
	    );
    if(thread_list) {
      for(uint thread = 0; thread < n_threads; thread++) { // Sets data for each thread.
	fprintf(f_log, "-\tThread[%d] is sheduled to read %.2fMB of data, starting at position %lld in the file, with a total length of %lld chars to read\n",
		thread, (float)thread_list[thread].get_total_length(arr)/((float)1024*1024),
		arr[thread_list[thread].get_index_start()], thread_list[thread].get_total_length(arr));
      } 
    } else fprintf(f_log, "(Threads have, when writing this file, not been allocated a subset of files to process.)\n");
    {
      fprintf(f_log, "#\tThe blocks are spread on a total of %lld containers:\n", position_in_arr);
      for(uint i = 0; i < position_in_arr+1; i++) {
	loint start_pos = last_file_pos_in_mem; 
	if(i != 0) start_pos = arr[i-1];
	const loint size = arr[i] - start_pos;
	fprintf(f_log, "-\tindex[%u] starts at %lld with length %lld)\n", i, arr[i], size);
      }      
    }
    // Generates the estimate on the time consumption:
    struct tms tmsend; clock_t end; long clktck = 0;
    if((end = times(&tmsend)) == -1) fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
    if((clktck =  sysconf(_SC_CLK_TCK)) < 0) fprintf(stderr, "!!\tAn error in sys-configuration at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
    const clock_t t_real = end - clock_log;
    const double log_time = t_real/(double)clktck;
    fprintf(f_log, "--\tUsed in total %10.8f seconds for this operation\n", log_time);
    fclose(f_log);
  }
#endif
#endif
}

/**
   @brief Initiates the file reading, using the data set. Thread safe.
   @param <_n_threads> The number of threads that weill be calling this function.
   @param <file_name> The name of the file to be read. (Should be the same file as block data are stored for.) 
   @remarks
   - To be called when the file is read- and the 'blocks' stored in the list named 'arr' are consistent.
   - Sets the index for the positions to retrieve in the objects of type thread_reading.
   @author Ole Kristian Ekseth (oekseth)
**/
void parse_read_blocks::init_thread_file_reading(uint _n_threads, char *file_name) {
  assert(file_name);
  assert(_n_threads);
  if(arr && (position_in_arr > 0)) {
    n_threads = _n_threads;
    const loint file_length = arr[position_in_arr-1];    
    loint avg_length = (file_length/n_threads);
    const uint MIN_SIZE = 6*1024; // If size is less than this base, only uses one thread
    // Some verifications: a ssummary is the line 'file_length/n_threads' below.
    if(avg_length < MIN_SIZE) {
      n_threads = (uint)file_length/MIN_SIZE;
      if(n_threads == 0) {
	n_threads = 1;
	if(1024 < file_length) avg_length = MIN_SIZE;
	else avg_length = file_length;
      } else avg_length = file_length/n_threads; // The total number of chars for each thread to read.
    }
    thread_list = new class thread_reading[n_threads]();
    assert(thread_list);
    uint curr_index = 0;
    uint start_pos = 0;
    uint end_pos = 0;
    for(uint i = 0; i < n_threads; i++) { // Sets data for each thread.
      loint size_block = 0; 
      if(i!= (n_threads-1)) {
	while((size_block < avg_length) && (curr_index < position_in_arr)) {
	  if(curr_index > 0) {
	    size_block += arr[curr_index] - arr[curr_index-1];
	  } else size_block = arr[0];
	  end_pos++, curr_index++;
	}
	thread_list[i] = thread_reading(disk_buffer_size, file_name, start_pos, end_pos, (loint)last_file_pos_in_mem); 
      } else {
	assert(end_pos <= position_in_arr);
	size_block = arr[position_in_arr-1];
	if(end_pos>0) size_block -= arr[end_pos-1]; // if more than one thread
	start_pos = curr_index;
	thread_list[n_threads-1] = thread_reading(disk_buffer_size, file_name, start_pos, position_in_arr, (loint)last_file_pos_in_mem);
	//end_pos = position_in_arr; // This index is never accesses,k terbhy it does not cause any access to unset memory.
	if(false)  { 	// TODO: Block code included for supplementing the documentation, illustrating the consepts of this compressed alogroithm.
	  printf("\n--\tIn parse_read_blocks at line %d:\n", __LINE__);
	  printf("\t\tThe arr for indexes above in range [%llu, %llu]\n", arr[start_pos], arr[position_in_arr]);
	  for(uint i = 0; i < position_in_arr+1; i++) {
	    loint start_pos = 0;
	    if(i) start_pos = arr[i-1];
	    else start_pos = last_file_pos_in_mem;
	    const int size = arr[i] - start_pos;
	    printf("\t\t\tarr[%u] = %llu giving size %d at start_pos(%llu)\n", i, arr[i], size, start_pos);
	    if(size > 0) {
	      FILE *file = fopen(file_name, "rb");
	      //	      if(file) {
		fseek(file, start_pos, SEEK_SET);
		char *block = new char[size+1];
		memset(block, '\0', size+1);
		const int read_size = fread(block, size, 1, file);
		assert(read_size <= size);
		printf("\n[\n");
		puts(block);
		printf("]\n");
	      delete [] block; fclose(file);
	    }
	  }
	}
      }
      //      end_pos = start_pos = end_pos;    
    }
  } else {
    n_threads = 0;
  }
}

// !Returns true if all the data in the fiels are read:
bool parse_read_blocks::all_data_read(bool print_result) {
  bool all_is_read = true;
  for(uint i = 0; i < (uint)n_threads; i++) {
    assert(print_result==false);
    bool temp = thread_list[i].file_is_read(i, print_result);
    if(!temp) all_is_read = false;
  }
  return all_is_read;
}

//! Returns the section of the data using the parse-blocks as a basis:
string_section *parse_read_blocks::get_string_section(uint thread_id) {
  if(thread_id < n_threads) {
    return thread_list[thread_id].get_string_section(thread_id, arr); 
  } else if(n_threads != 0) {
    fprintf(stderr, "!!(parse_read_blocks.cxx at line %d)\tThread id(%u) above reserved threads(%u). An error. Please contact the developer at email oekseth@gmail.com\n", __LINE__, thread_id, n_threads);
    assert(false);
    return NULL;
  } else return NULL;
}

//! Frees the memory:
void parse_read_blocks::free_mem() {
  if(arr) {free(arr); arr=NULL; size = 0, position_in_arr = 0;}
  thread_reading_t::free_mem(thread_list, n_threads), delete [] thread_list, thread_list = NULL, n_threads = 0;
}

parse_read_blocks::parse_read_blocks(uint _disk_buffer_size) : disk_buffer_size(_disk_buffer_size), size(20), position_in_arr(0), read_postion_to_retrieve(0), thread_list(NULL), n_threads(0), last_file_pos_in_mem(0) {
  arr = (loint*)malloc(sizeof(loint)*size);
  log_builder::test_memory_condition_and_if_not_abort(arr!=NULL, __LINE__, __FILE__, __FUNCTION__);
  for(loint i = 0; i< size; i++) arr[i] = 0;
}

#ifdef assert_code
//! Asserts the trehad reading
void parse_read_blocks::assert_thread_file_reading() {
  const uint _n_threads = 3;
  init_thread_file_reading(_n_threads, "all.blast");
  assert(thread_list);
  loint sum = 0;
  for(uint i = 0; i < n_threads; i++) {
    const loint tot = thread_list[i].get_total_length(arr);
    sum+= tot;    
  }
  assert(arr);
  assert(sum == arr[position_in_arr-1]);
}
#endif

void parse_read_blocks::assert_class(const bool print_info) {
  const static char *class_name = "parse_read_blocks";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  // Iinitial tests:
  uint length = 10;
  parse_read_blocks temp = parse_read_blocks(1024);
  for(uint i = 0; i < length; i++) {
    temp.insert(i+1); // Sets the psostion in the file: must be greater than one
  }
  // TODO: Below code depricated due to lbock-size set at compilation level. Consider change, or alternative remove coce below.
  /*
    for(uint i = 0; i < length; i++) {
    const uint curr = (uint)temp.getNextReadLength();  
    //    assert((uint)1 == curr);
    }*/
  temp.assert_thread_file_reading();
  temp.free_mem();
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

