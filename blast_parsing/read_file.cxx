#include "read_file.h"

//! Allocates and initializes the read-buffer to be sent downward the pipe:
void read_file::allocate(char *&buffer_main, char *&buffer_end, const long int buffer_in_mem_size) {
  buffer_main = new char[buffer_in_mem_size+1];//new char[buffer_in_mem_size + 1];
  log_builder::test_memory_condition_and_if_not_abort(buffer_main!=NULL, __LINE__, __FILE__, __FUNCTION__);
  //  buffer_main = (char*)malloc(sizeof(char)*(buffer_in_mem_size+1));//new char[buffer_in_mem_size + 1];
  buffer_end = buffer_main + buffer_in_mem_size;
  memset(buffer_main, '\0', buffer_in_mem_size+1);// Filling in the buffer:
}

/**
   @brief Reads the data
   @param <str>    The char buffer to put things into.
   @param <size>   The bytes each of the elements occupie eg, for char, this implies '1'
   @param <nmemb>  The maximum number of elements to read.
   @param <stream> The stream to read data from.
   @return the number of elements read.
**/
static size_t f_read(char *str, size_t size, size_t nmemb, FILE* stream) {
  const size_t disk_b = BLOCK_FILE_READ_SIZE;
  assert(str);
  if(size > disk_b) {
    loint block_cnt = size/((size_t)BLOCK_FILE_READ_SIZE);
    size_t read_cnt = 0; // the position in the string
    for(uint i = 0; i < block_cnt-1; i++) {
      read_cnt += fread(str+read_cnt, BLOCK_FILE_READ_SIZE, nmemb, stream);
    }
    read_cnt += fread(str+read_cnt, (size - read_cnt), nmemb, stream);
    return read_cnt;
  } else return fread(str, size, nmemb, stream);
}
//! Reads the file:
void read_file::read_from_file(char *buffer_main, const long int buffer_in_mem_size) {
  if(file_input) {f_read(buffer_main, sizeof(char), buffer_in_mem_size, file_input);}
#ifdef USE_MPI
  else {
    assert(mpi_obj);
    mpi_obj->read_chars(buffer_main, buffer_in_mem_size);
  }
#endif
}

//! Prints the data in the buffer:
void read_file::print_buffer(char *buffer_main, char *buffer_end) {
  char *word = buffer_main;
  printf("\nThe block read is:\n|");
  while(word < buffer_end) putchar(*word++);
  printf("\n|\n");
}

bool read_file::read_from_file(char *buffer_main, long int buffer_in_mem_size, long int buffer_rest_size, long int &size_current_buffer) {
  const long int chrs_to_read = buffer_in_mem_size - buffer_rest_size;    
  long int fread_cnt =0;
  if(mpi_obj) assert(!file_input);

  if(file_input) {
    fread_cnt = (long int)f_read(buffer_main+buffer_rest_size, sizeof(char), chrs_to_read, file_input);
  } 
#ifdef USE_MPI
  else {
    assert(mpi_obj);
    fread_cnt = (long int)mpi_obj->read_chars(buffer_main+buffer_rest_size, chrs_to_read);
  }
#endif
  size_current_buffer = buffer_rest_size + fread_cnt;
  assert(size_current_buffer <= buffer_in_mem_size);
  bool end_of_file = false;       // Marks if its the end of the file
  if(fread_cnt < chrs_to_read) end_of_file = true;
  return end_of_file;
}

//! Finding the last new_line char:
char *read_file::get_logical_end(char *buffer_main,char *buffer_main_end, const bool end_of_file, long int &buffer_rest_size) {
  char *buffer_main_logical_end = NULL;
  if(!end_of_file) {
    buffer_main_logical_end = strrchr(buffer_main, '\n');
    if(buffer_main_logical_end != NULL) {
      buffer_rest_size = buffer_main_end - buffer_main_logical_end-1; 
      //      buffer_rest_size = buffer_main_end - buffer_main_logical_end -1; // to remove the newline
    } else buffer_rest_size = 0;
  } else {
    buffer_rest_size = 0;
    buffer_main_logical_end = buffer_main_end; // At the end of the file sets it to the end
  }
  return buffer_main_logical_end;
}

//! Copies the owerflowing data from the previous buffer to this:
void read_file::copy_overflow(char *buffer_main, char *buffer_rest, long int buffer_rest_size) {
  if (buffer_rest_size != 0) { // there are chars to insert from the previous retrival.
    memcpy(buffer_main, buffer_rest, sizeof(char)*buffer_rest_size); // Adds char from the previous
  }
}


//! Updates the buffer with data to be included on the next run:
void read_file::update_overflow_buffer(char *buffer_main_start, char *buffer_rest, char *buffer_rest_start, char *buffer_main_logical_end, long int buffer_rest_size, long int buffer_rest_total_size, long int size_current_buffer) {
  buffer_rest_start = buffer_rest;
  memset(buffer_rest_start, '\0', buffer_rest_size+1);
  if (buffer_rest_size == 0) {
    buffer_main_logical_end = buffer_main_start + size_current_buffer;
  } else { // Has chars outside the boundary
    memset(buffer_rest_start, '\0', buffer_rest_total_size); 
    buffer_main_logical_end+=1;
    memcpy(buffer_rest_start, buffer_main_logical_end, sizeof(char) * buffer_rest_size);
  }
}  

/**
   @brief Sets the position of the file pointer
   @return true if operation was possible
**/
bool read_file::set_file_pointer_position(loint file_position) {
  if(-1 != fseek(file_input, file_position, SEEK_SET)) return true; // positions the file
  else return false;
}

/**
   @brief Remove content of overflow buffer.
   @remarks Useful when jumpring around in a file, only retrieving the content, ie, not merging it     
**/
void read_file::set_overflow_buffer_to_emtpy() {
  buffer_rest_size = 0;
//   if(buffer_rest_size && buffer_rest) {
//     memset(buffer_rest, '\0', buffer_rest_size);
//   } 
}
//! Verifies the correctness of the block
void read_file::verify_correctness_of_block(int myrank, char *start, char *buffer_end, char seperator, const int line_caller, const char *file_caller) {
  assert(start);
  assert(buffer_end);
#ifndef NDEBUG
  const bool verbouse = false; // Handy if errors- or changes for the future happends.
  if(verbouse) printf("Verifies that all the data sent is correct, given data from line %d in file %s\n", line_caller, file_caller);
  while(start < buffer_end - 5) {
    char *line_start = start;
    char *end = strchr(start, '\n');
    if(end) {
      if(verbouse) printf("has the following chars: |\n");
      uint sep_cnt = 0;
      while(start != end) {
	if(*start == seperator) sep_cnt++;
	if(verbouse) putchar(*start);
	start++;
      } 
      if(verbouse)       printf("\n|\n");
      if(sep_cnt!=2) {
	fprintf(stderr, "!!\t Identified a difficult row in the blastp file; suggest you remove it, or fix the problems with the fields of it. To help identify details, we will now include details which should be forwarded to [oekseth@gmail.com]:\n[myrank=%d]\tPlease verify the setting used launching this software. (Merging of blocks in blast_p failed: Found #(sperators)=%u, at %s:%d, originally called from %s:%d\n", myrank, (uint)sep_cnt, __FILE__, __LINE__, file_caller, line_caller);
	log_builder::throw_warning(software_error, __LINE__, __FILE__, __FUNCTION__, "handling of blastp-row failed");
	//! Print the last line: important during debugging, ie..e understanding of the problem. Last tie it helped solving a problem was April 26, 2013.
	char *line_end = strchr(line_start, '\n');
	assert(line_end);
	printf("[%d]\tThe line (length=%d) causing the trouble (sep_cnt = %d and seperator=%c) is:\n|\n", myrank, (int)(line_end-line_start), sep_cnt, seperator);
	while(line_start != line_end) {putchar(*line_start); line_start++;}
	printf("\n|\n");
	fprintf(stderr, "\nSummary:\n-\t The error is either due to (a) bugs in the BLASTp-file you provided or (b) that the software did not handle it (ie, your input-file). The BLAST program has a tendency generating format-errors in its result. It's easy to chech if the error is due to your blast-input-file: if the above printed line is found in your BLAST-input, the BLAST-input is wrong. Else the problem is with this software. Please contact the developer at [oekseth@gmail.com] if you did not understand this error-message. The error-message was generated at [%s]:%s:%d.\n", __FUNCTION__, __FILE__, __LINE__);
	//	assert(false); 
      }
      start = end +1;
    } else start = buffer_end;
  }
#endif
}
/**
   Reads the file and dump a block of chars to the next element in the pipe
   @return An object holding the chars read.
**/
struct string_section *read_file::read_file_blocks(char seperator, uint &block_cnt, const bool perform_tests) {
  const long int buffer_in_mem_size = disk_buffer_size; 
  // Holding the grabbed strign
  char *buffer_main = NULL; 
  try {buffer_main = new char[buffer_in_mem_size+1];}  
  catch (std::exception& ba) {
    if(!log_builder::catch_memory_exeception(buffer_in_mem_size+1, __FUNCTION__, __FILE__, __LINE__)) {
      fprintf(stderr, "!!\t An interesting error was discovered: %s."
	      "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
	      "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
    } 
  }
  //  char *buffer_main = new char[buffer_in_mem_size+1];//new char[buffer_in_mem_size + 1];
  log_builder::test_memory_condition_and_if_not_abort(buffer_main!=NULL, __LINE__, __FILE__, __FUNCTION__);
  char *buffer_main_start = buffer_main; // holds the start position
  char *buffer_main_end = buffer_main + buffer_in_mem_size; // holds the start position
  char *buffer_rest_start = buffer_rest;
  bool file_is_not_empty = true;
  //! Assures that there are more data to read.
  if(file_input) {
    if(feof(file_input)) file_is_not_empty = false;
  } else {
#ifdef USE_MPI
    if(mpi_obj) {if(!(mpi_obj->has_more_data_in_file_to_read())) {file_is_not_empty = false;}}
    else 
#endif
      file_is_not_empty = false;
  }
  if(file_is_not_empty) { // Assuring not the end of the file (27.10.2011 by oekseth)
    memset(buffer_main, '\0', buffer_in_mem_size+1);// Filling in the buffer:
    copy_overflow(buffer_main, buffer_rest_start, buffer_rest_size); // Copies the owerflowing data from the previus
    //    copy_overflow(buffer_main, buffer_rest_start, buffer_rest_size); // Copies the owerflowing data from the previus
    long int size_current_buffer = 0;
    const bool end_of_file = read_from_file(buffer_main, buffer_in_mem_size, buffer_rest_size, size_current_buffer);
    buffer_main_end = buffer_main + size_current_buffer; 
    if(perform_tests)    {      
      verify_correctness_of_block(myrank, buffer_main_start, buffer_main_end, seperator, __LINE__, __FILE__);
    }

    buffer_rest_size = 0; // resets;    
    // Finding the last new_line char:
    char *buffer_main_logical_end = get_logical_end(buffer_main, buffer_main_end, end_of_file, buffer_rest_size);
    if(end_of_file || buffer_main_logical_end != NULL) {
      // Updates the size and the data of the part to be copied to the next run:      
      long int size_this_buff = (buffer_main_logical_end - buffer_main_start); 
      if(*buffer_main_logical_end == '\n') {
	//	printf("er lik newline\n");
	size_this_buff++;
      }
      update_overflow_buffer(buffer_main_start, buffer_rest, buffer_rest_start, buffer_main_logical_end, buffer_rest_size, buffer_rest_total_size, size_current_buffer);
      
      if(perform_tests)    verify_correctness_of_block(myrank, buffer_main_start, buffer_main_logical_end, seperator, __LINE__, __FILE__);

      //! Sends the buffer (Starting on a line start and ending on a newline) to the buffer
      string_section *section = string_section::init(buffer_main_start,
				  buffer_main_logical_end,
				  chars_sendt_to_parsing, block_cnt, end_of_file);  
      //      section->print_segment(buffer_main_start, buffer_main_logical_end);
      //      section->print_data_block();
      block_cnt++; // Increments for the next iteration
      chars_sendt_to_parsing += size_this_buff; // 
      return section;
    } 
  }
  // At this point, the file is at its end:
  delete [] buffer_main;
  return NULL;
} 
 
//! Constructor:
read_file::read_file(uint _disk_buffer_size, char *f_name) :
  disk_buffer_size(_disk_buffer_size), myrank(0), mpi_obj(NULL) {
  if(f_name != NULL) {
    file_input = fopen(f_name, "rb"); // assumes that reading in binary mode goes faster
    log_builder::test_file_reading_and_if_not_abort(file_input, f_name, __LINE__, __FILE__, __FUNCTION__);
  } else {
    char string[100]; memset(string, '\0', 100); sprintf(string, "File not set for name '%s'", f_name);
    log_builder::throw_error_and_abort(string, __LINE__, __FILE__, __FUNCTION__);
  }
  buffer_rest_total_size = disk_buffer_size;
  BUFFER_REST = new char[buffer_rest_total_size];
  memset(BUFFER_REST, '\0', buffer_rest_total_size);
  buffer_rest = &BUFFER_REST[0];
  buffer_rest_size = 0;
  chars_sendt_to_parsing = 0;

}

#ifdef USE_MPI
//! The constructor for MPI-procedure-initiation.
read_file::read_file(uint _disk_buffer_size, MPI_Comm comm, int _myrank, char *f_name, loint reading_file_start_position, loint reading_file_end_position, loint reading_file_length) :
  disk_buffer_size(_disk_buffer_size),
  myrank(_myrank),
  file_input(NULL)
{
  mpi_obj = new mpi_read(comm, myrank, f_name, reading_file_start_position, reading_file_end_position, reading_file_length, disk_buffer_size);
  assert(mpi_obj);
  buffer_rest_total_size = disk_buffer_size;
  BUFFER_REST = new char[buffer_rest_total_size];
  memset(BUFFER_REST, '\0', buffer_rest_total_size);
  buffer_rest = &BUFFER_REST[0];
  buffer_rest_size = 0;
  chars_sendt_to_parsing = 0;

}
#endif
//! To be called after the operation is finalized
void read_file::free_mem() {
  if(file_input) fclose(file_input);
#ifdef USE_MPI
  mpi_read::close(mpi_obj);  
#endif
  if(BUFFER_REST) {delete [] BUFFER_REST; BUFFER_REST = NULL;}
}
void read_file::assert_private_parts() {
#ifdef assert_code
  // Tests overlap filtering.
  char *buffer_main = NULL, *buffer_end = NULL; // Holding the grabbed string.
  const long int buffer_in_mem_size = 450;
  allocate(buffer_main, buffer_end, buffer_in_mem_size);

  free(buffer_main); // Sould normally be freed in a pipe-class further down the slope.
#endif
}
void read_file::assert_class(const bool print_info) {
  const static char *class_name = "read_file";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  // Iinitial tests:
  char *test_file = "test_file.txt";
  uint size_line = 100, cnt_newlines = 4,  cnt_block = 2;
  char *test_block = new char[size_line*cnt_newlines+1];
  const char seperator = '|';
  for(uint i = 0; i < (size_line*cnt_newlines); i++) test_block[i] = 'a';  
  for(uint i =1; i < cnt_newlines; i++) {
    test_block[((i-1)*size_line) + 4] = seperator;
    test_block[((i-1)*size_line) + 8] = seperator;
    test_block[i*size_line-1] = '\n';
  }
  test_block[size_line*cnt_newlines-1] = '\n';

  test_block[size_line*cnt_newlines] = '\0';
  FILE *file = fopen(test_file, "w");
  uint cnt_chars=0;
  if(file != NULL) {
    for(uint i = 0; i < cnt_block; i++) {
      fputs(test_block, file);
      cnt_chars +=strlen(test_block);
    }
    fclose(file);
  }  
  read_file temp(1024, test_file);
  uint block_cnt = 0;
  //  printf("found at line %d in read_file.cxx %d elements.\n", __LINE__, (uint)temp.read_file_blocks(block_cnt));
  string_section *section = temp.read_file_blocks(seperator, block_cnt, false);
  //  printf("writes %d blocks, with a total size of %u in and cnt_chars=%u in read_file.cxx\n", cnt_block, size_line*cnt_newlines+1, cnt_chars);
  //  printf("block_cnt=%u and cnt_block = %u in read_file.cxx\n", block_cnt, cnt_block);
  //  printf("found at line %d in read_file.cxx with length(string)=%d.\n", __LINE__, (int)section->get_data_length());
  assert(cnt_chars == section->get_data_length()); // 4 lines excluding the newline
  //  assert(399 == section->get_data_length()); // 4 lines excluding the newline
  //section->finalize();
  string_section::close(section);

  section = temp.read_file_blocks(seperator, block_cnt, false);

  //  printf("found at line %d in read_file.cxx %d elements.\n", __LINE__, (int)section->get_data_length());
  // // assert(399 == section->get_data_length()); // 4 lines excluding the newline
  // // section->finalize();
  // // free(section);

  // // section = temp.read_file_blocks(block_cnt);

  // // if(section != NULL) {
  // //   assert(0 == section->get_data_length()); // 4 lines excluding the newline
  // //   section->finalize();
  // //   free(section);
  // // }


  
  delete [] test_block;
  temp.free_mem();
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
