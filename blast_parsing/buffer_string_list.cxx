#include "buffer_string_list.h"


//! Merges two strings:
void buffer_string_list::merge(char *&str, loint &length, char *&temp, loint temp_length) {
  if(str== NULL) {str = temp, length = temp_length;}
  else { 
    const loint res_size = temp_length + length;
    char *res = new char[res_size+1]; res[res_size] = '\0';
    assert(res);
    memcpy(res, str, length); memcpy(res + length, temp, temp_length);
    delete [] temp;
    delete [] str; str = res; length = res_size;
  }
}

/**
   @brief Generate (writes) a log file describing the details of the memory signature for strings stored in memory during parsing
   @remarks
   - Useful for analyzing- and optimizing the memory allocation procedures.
   - In optimised mode does not produce the log file.
**/
void buffer_string_list::log_generate_memory_allocation_overview(uint n_threads) {
#ifndef NDEBUG
  assert(n_threads);
#ifdef LOG_WRITE_MEMORY_ALLOCATIONS_DURING_BLASTP_STRINGS_STORED_IN_MEMORY
  assert(list);
  struct tms tmsstart; clock_t clock_log;
  if((clock_log = times(&tmsstart)) == -1) // starting values
    fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
  char sep = '_';
  // The log file, and contents of it:
  if(sizeof(FUNC_LOG_FOLDER_NAME) > 0) sep = '/';
  char str[200]; memset(str, '\0', 200);
  sprintf(str, "%s%c%s.log.%u", FUNC_LOG_FOLDER_NAME, sep, "buffer_string_list_memory", n_threads);
  FILE *f_log = fopen(str, "w");
  if(f_log) {
    loint mem_bytes_usage_first = 0, mem_bytes_usage_second = 0;
    for(loint i = 0; i < size_index; i++) {
      mem_bytes_usage_first += list[i].get_size_first_buffer();
      mem_bytes_usage_second += list[i].get_size_second_buffer();
    }
    fprintf(f_log, "Stored %.2fMB in memory, spreading it over %.2fMB on the first blocks and %.2fMB on the second blocks, using a total of %llu objects of type buffer_string.\n",
	    (float)(mem_bytes_usage_first+mem_bytes_usage_second)/(float)(1024*1024), 
	    (float)(mem_bytes_usage_first/(float)(1024*1024)), 
	    (float)(mem_bytes_usage_second/(float)(1024*1024)),
	    size_index
	    );
    fprintf(f_log, "#\tThe precise distribution of the blocks were:\n");
    for(loint i = 0; i < size_index; i++) {
      const loint mem_bytes_usage_first = list[i].get_size_first_buffer();
      const loint mem_bytes_usage_second = list[i].get_size_second_buffer();
      fprintf(f_log, "-\tAt index[%lld] stored %.2fMB in memory, spreading it over %.2fMB on the first blocks and %.2fMB on the second blocks.\n",
	      i,
	      (float)(mem_bytes_usage_first+mem_bytes_usage_second)/(float)(1024*1024), 
	      (float)(mem_bytes_usage_first/(float)(1024*1024)), 
	      (float)(mem_bytes_usage_second/(float)(1024*1024)));
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

//! Prints the properties of the item last inserted.
void buffer_string_list::print_item_last_inserted() {
  if(list) {
    printf("Object at index[%llu] has data: ", current_index);
    list[current_index].print_info();
  } else printf("(no itmes set for buffer_string_list object.)\n");
}

//! Prints the properties of the items in the list
void buffer_string_list::print_items() {
  if(list) {
    printf("has %llu elements in the list:\n", size_index);
    for(loint i = 0; i < size_index; i++) {
      printf("Object at index[%llu] has data: ", i);
      list[i].print_info();
    }
  } else printf("(no itmes set for buffer_string_list object.)\n");
}
//! Returns the string given, using the blast file stored in memory.
void buffer_string_list::get_string(loint &length, char *&str, loint &found_at_index_pos) {
  assert(list);
  found_at_index_pos = current_index;
  if(current_index < size_index) {
    if(list[current_index].has_first_block() || list[current_index].has_main_block()) {
      loint tmp_length;
      char *tmp = NULL; list[current_index].get_string(tmp_length, tmp);
      list[current_index].free_mem((uint)current_index);
      merge(str, length, tmp, tmp_length);
      current_index++;
    } else if (current_index != (size_index-1)){ // iterate until a main-block is found:
      loint temp_length = 0;
      char *temp = NULL; list[current_index].get_string(temp_length, temp); // recursive
      current_index++;
      merge(str, length, temp, temp_length);
      get_string(length, str, found_at_index_pos);	
    } else { // At the end, but the main block not found.
      // Description: A situation like this may occure at the start of the file
      // reading when filepos [100, 300] is the the memoery while the file structure is
      // requested reading filepos [100, 600].      
      if(true) {
	list[current_index].free_mem((uint)current_index); if(str) delete [] str; length = 0, str = NULL;
      } else {
	// An alternative solution, as a backup if future inputs do not follow our expectations
	loint tmp_length;
	char *tmp = NULL; list[current_index].get_string(tmp_length, tmp);
	list[current_index].free_mem((uint)current_index);
	merge(str, length, tmp, tmp_length);
	current_index++;
      }
      current_index++;
    }
  } 
}

 
//! Returns false if no mare data is of the file in memory
bool buffer_string_list::get_section(string_section *&section, loint &found_at_index_pos) {
  char *str = NULL; loint length = 0; this->get_string(length, str, found_at_index_pos);
  if(str) {
    section = new string_section(str, str + length, 0, 0, false);
    return true;
  } else {
    return false;
  }
}
//! Inserts the overflow buffer at the correct location:
void buffer_string_list::copy_into_first_buffer(uint block_cnt, char *buff, loint size) {
  if(size > 0) {
    if((loint)block_cnt <size_index) {
      list[block_cnt].copy_into_first_buffer(buff, size);
    } else {
      fprintf(stderr, "!!\tAn error in class buffer_string_list at line %d: (%u, %llu)\n", __LINE__, block_cnt, size); fflush(stderr);
      assert(false);
    }
  }
}

//! Inserts the main buffer at the correct location:
void buffer_string_list::set_second_buffer(uint block_cnt, char *buff, loint size) {
  if(size > 0) {
    if((loint)block_cnt <size_index) {
      list[block_cnt].set_second_buffer(buff, size);
    } else {
      fprintf(stderr, "!!\tAn error in class buffer_string_list at line %d: (%u, %llu)\n",   __LINE__, block_cnt, size); fflush(stderr);
      assert(false);
    }
  }
}


/**
   @brief Includes the data given into the buffer available for all the threads.
   @param <section> The object containing the chars from the blast file to store.
   @param <block_length> The remaining length of data to be merged with the block following this data in the input blast file.
   @return false if data is not inserted into memory.
   @remarks If the index block provided is above the length of the list reserved, data is not stored in this object. This implies that alternative method of getting the data is required. In the turboOrto project class thread_reading provides this option. For examples of usage, see class pipe_parse_parse.
   @author Ole Kristian Ekseth (oekseth)
   @data 05.01.2012 by oekseth (documentation and debuging).
**/
bool buffer_string_list::update(string_section_t *&section, loint block_length) {
  assert(section);
  if((loint)(section->block_cnt) < size_index) {      
    // 1. Updates the stringBuffer
    if((loint)(section->block_cnt+1) < size_index) {      
      loint rest_length = section->get_data_length() - block_length;
      copy_into_first_buffer(section->block_cnt+1, section->buffer_main_start+block_length, rest_length);
    } else {
      fprintf(stderr, "!!\t\tan error occured in %s at line %d: index[%d] above the reserved limit(%d). If in debug mode, aborts. Else the last block of the file is not included in the filtering process. Therefore recomens you contact oekseth@gmail com providing this error message and infromation about the file size used.\n", __FILE__, __LINE__, (int)section->block_cnt+1,  (int)size_index); 
      assert(false);
    }
    set_second_buffer(section->block_cnt, section->buffer_main_start, block_length);
    // 2. Removes the pointer in the section avoiding memory to be cleared
    section->buffer_main_start = NULL;
    section->buffer_main_logical_end = NULL;
    return true;
  } else {
    return false;
  }
}
//! Returns the file size.
loint buffer_string_list::get_file_size(char *file) { 
  if(file) {
    struct stat sb;
    if (stat(file, &sb) == -1) { 
      fprintf(stderr, "!!\tUnable to open file '%s' at line %d in class buffer_string_list: ", file, __LINE__);
      perror("Exits");
      exit(EXIT_FAILURE);
    }
    const loint size = (loint)sb.st_size;
    return size;
  } else {
    fprintf(stderr, "!!\tFile(%s) not defined at line %d in class buffer_string_list.cxx\n", file, __LINE__); fflush(stderr);
    assert(false);
    return 0;  
  }
}

loint buffer_string_list::get_max_size_in_buff(char *f_name) {
  const loint file_size = buffer_string_list::get_file_size(f_name);
  // TODO: Consider using different tipes of memory-usages-chemase (LOW, MEDIM, GREADY
    
  loint temp_size = (loint)(file_size * 1.05); // Estimates that approx 1/20 of the file size will be used reserving things in memory later on. 
  char *temp = NULL;
  try {temp = new char[temp_size];} 
  catch (std::exception& ba) {
    if(!log_builder::catch_memory_exeception(temp_size, __FUNCTION__, __FILE__, __LINE__)) {
      fprintf(stderr, "!!\t An interesting error was discovered: %s."
	      "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
	      "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
    }
  }
  //  char *temp = new char[temp_size];
  if(!temp) while(temp == NULL && (temp_size > 0)) {
      temp_size = (loint)(temp_size * 0.5);
      //      temp = new char[temp_size];
      try {temp = new char[temp_size];} 
      catch (std::exception& ba) {
	if(!log_builder::catch_memory_exeception(temp_size, __FUNCTION__, __FILE__, __LINE__)) {
	  fprintf(stderr, "!!\t An interesting error was discovered: %s."
		  "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		  "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
	}
      }
    } 
  if(temp != NULL) delete [] temp;  
  return temp_size;
}
//! Returns an estimated number of the block_count
loint buffer_string_list::get_estimate_of_memory_blocks_cnt(uint disk_buffer_size, char *f_name) {
  if(f_name) {
#ifdef MEMORY_CONSUMPTION_LEVEL
    //if(MEMORY_CONSUMPTION_LEVEL > 0)
    {
      const loint dbsize = (loint)disk_buffer_size;
      const loint max_file_size = buffer_string_list::get_max_size_in_buff(f_name);
      if(max_file_size > 0) {
	const loint upper_size = (loint)(1+(max_file_size*1.1));
	loint block_cnt = 3;
	const loint estimated_line_size = 150; // The worst case owerflow.
	assert(dbsize > estimated_line_size);
	if(upper_size > dbsize) {
	  block_cnt = (upper_size/(dbsize-estimated_line_size)); // Uses a factor to be on the safe side.
	}
	// Below two blocks added as an owerflow for both ordinarty- and extra is needed.
	block_cnt+=4; // To allow for the last owerflow block:
	//	printf("in buffer_string_list.cxx extimated block_cnt is %d blocks, oor more rpecicely %f blocks\n", (int)block_cnt, (float)(upper_size/dbsize));
	//	printf("in buffer_string_list.cxx this due to (%d/%d) ~ %f\n", (int)max_file_size, (int)dbsize,	       (float)(max_file_size/dbsize));
	if(block_cnt < 3) {	  
	  block_cnt = 3;	  
	}	
	return block_cnt;
      } else return 0;
    } //else return 0;
#else
    return 0;
#endif
  } else fprintf(stderr, "!!\tFile not specified (input=NULL) in class buffer_string_list at line %d. Contact author of the code if thsi error is seen.\n", __LINE__);
  return 0;
}

buffer_string_list::buffer_string_list(loint size) : current_index(0), size_index(size) {
  list = buffer_string::init(size);
  assert(list);
  //  fprintf(stderr, "in buffer_string_list at line %d initialised the 'list' with %d elements\n", __LINE__, (int)size);
}
void buffer_string_list::close(buffer_string_list &tmp) {tmp.free_mem();}
void buffer_string_list::close(buffer_string_list *&tmp) {
  if(tmp) {tmp->free_mem(), delete tmp, tmp= NULL;}
}
void buffer_string_list::free_mem() {
  if(list) buffer_string::close(list, size_index);
  current_index = 0;
}

  
void buffer_string_list::assert_class(const bool print_info) {
  const static char *class_name = "buffer_string_list";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  const loint size = 10;
  {
    buffer_string_list tmp = buffer_string_list(size);
    // Tests when all the blocks contains both first- and second buffer:
    char *a = new char[10]; memset(a, 'a', 10);
    for(loint i = 0; i < size; i++) {
      tmp.copy_into_first_buffer((uint)i, a, 10);
      char *b = new char[20]; memset(b, 'b', 20);
      tmp.set_second_buffer((uint)i, b, 20);
    }
    char *corr = new char[30]; memset(corr, 'a', 10);  memset(corr+10, 'b', 20);
    for(loint i = 0; i < size; i++) {
      loint len = 0;
      loint temp = 0;
      char *res = NULL; tmp.get_string(len, res, temp);
      //      assert(0 == strncmp(res, corr, 30));
      //assert(len == 30);
      delete [] res;
    }
    delete [] a, delete [] corr;    
    tmp.free_mem();
  }
  // Tests the merging (1):
  {
    buffer_string_list tmp = buffer_string_list(size);
    tmp.copy_into_first_buffer((uint)0, "a", 1);
    char *b = new char[20]; memset(b, 'b', 20);
    tmp.copy_into_first_buffer((uint)1, "a", 1);
    tmp.set_second_buffer((uint)1, b, 20);
    loint len = 0;
    char *corr = new char[22]; memset(corr, 'a', 2);  memset(corr+2, 'b', 20);
      loint temp = 0;
      char *res = NULL; tmp.get_string(len, res, temp);
    // Her slettes '0'??
    //assert(len == 22);
    //    assert(0 == strncmp(res, corr, 22));
    tmp.free_mem();
    delete [] corr;
    delete [] res;
  }
  {
    // Tests the merging (2):
    buffer_string_list tmp = buffer_string_list(size);
    for(loint i = 0; i < 3; i++) {
      tmp.copy_into_first_buffer((uint)i, "a", 1);
    }
    char *b = new char[20]; memset(b, 'b', 20);
    tmp.copy_into_first_buffer((uint)3, "a", 1);
    tmp.set_second_buffer((uint)3, b, 20);
    loint len = 0;
    char *corr = new char[24]; memset(corr, 'a', 4);  memset(corr+4, 'b', 20);
      loint temp = 0;
      char *res = NULL; tmp.get_string(len, res, temp);
    //    assert(len == 24);
    //    assert(0 == strncmp(res, corr, 24));
    tmp.free_mem();
    delete [] corr;
    delete [] res;
  }
  // delete [] corr;
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
  
