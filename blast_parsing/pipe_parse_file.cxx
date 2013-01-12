#include "pipe_parse_file.h"

void print_seg(char *start_p, char *end) {
  //  printf("The string:\n|\n");
  char *start = start_p;
  while(start < end) putchar(*start), start++;
  printf("\n");
  //  printf("\n|\n");
}


//! Returns the section of the data using the parse-blocks as a basis:
string_section *pipe_parse_file::get_section() {
  FILE *file = fileRead->file_input;
  return parseBlocks->get_string_section_from_file(file);
}
 
//! Sets the file to the start_pos of the file
void pipe_parse_file::rewind_file(parse_read_blocks *_parseBlocks) {
  assert(fileRead);
  assert(_parseBlocks);
  fileRead->rewind_file();
  block_cnt = 0; // resets the block count
  parseBlocks = _parseBlocks;
  FIRST_RUN = false;
}

/**
   @return An object consisting of a char-block with extra meta-data.
   @remark Using method found in object of type read_file.
**/
struct string_section *pipe_parse_file::read_file_blocks() {
  return fileRead->read_file_blocks(seperator, block_cnt, true); 
}
 
int cnt = 0;
/*overrride*/void* pipe_parse_file::operator()(void* item) {
  if(false) printf("i=%ld", (long int)item); // Hiding an unsed item from teh compiler
  //  printf("starts\tat line %d in file %s\n", __LINE__, __FILE__);
  const clock_t tstart = times(NULL);
  string_section_t *section = NULL;
  if(FIRST_RUN)  {
    section = read_file_blocks();
    if(section) {
      chars_in_file_read_first_run += (loint)section->get_data_length();
      if(!section->end_of_file) chars_in_file_read_first_run++;
    }
    const clock_t tend = times(NULL); 
    log->append_measurement(read_first, 0, tstart, tend); 
  } else {
    section = get_section();
    const clock_t tend = times(NULL); 
    log->append_measurement(read_second, 0, tstart, tend);
  }
  //  debug_sum_time_usage(false, (tend-tstart).seconds());
  //  printf("returns\tat line %d in file %s\n", __LINE__, __FILE__);
  return section;
}

//! Constructor:
pipe_parse_file::pipe_parse_file(char *f_name, log_builder_t *_log) : tbb::filter(  /*serial=*/ true), seperator('_'),  disk_buffer_size(1024*1024*10), log(_log), block_cnt(0), parseBlocks(NULL), FIRST_RUN(true), FILE_INPUT_NAME(f_name), chars_in_file_read_first_run(0), myrank(0) {
  fileRead = read_file::init(1024*1024*10, f_name);
}


//! Constructor:
pipe_parse_file::pipe_parse_file(char _seperator, uint _disk_buffer_size, char *f_name, log_builder_t *_log) : tbb::filter(  /*serial=*/ true), seperator(_seperator), disk_buffer_size(_disk_buffer_size), log(_log), block_cnt(0), parseBlocks(NULL), FIRST_RUN(true), FILE_INPUT_NAME(f_name), chars_in_file_read_first_run(0), myrank(0) {
  fileRead = read_file::init(disk_buffer_size, f_name);
}

#ifdef USE_MPI
//! Constructor
pipe_parse_file::pipe_parse_file(char _seperator, uint _disk_buffer_size, MPI_Comm comm, int _myrank, char *f_name, log_builder_t *_log, loint reading_file_start_position, loint reading_file_end_position, loint reading_file_length) : tbb::filter(  /*serial=*/ true), seperator(_seperator), disk_buffer_size(_disk_buffer_size), log(_log), block_cnt(0), parseBlocks(NULL), FIRST_RUN(true), FILE_INPUT_NAME(f_name), chars_in_file_read_first_run(0), myrank(_myrank) {
    if(true) {
      fileRead = read_file::init(disk_buffer_size, comm, myrank, f_name, reading_file_start_position, reading_file_end_position, reading_file_length);
    } else {
      if(reading_file_length) fileRead = read_file::init(disk_buffer_size, comm, myrank, f_name, reading_file_start_position, reading_file_end_position, reading_file_length);
      else   fileRead = read_file::init(disk_buffer_size, f_name);
    }
}
#endif

//! To be called after the operation is finalized
void pipe_parse_file::close_file(int n_threads) {
  if(fileRead) {read_file::close(fileRead);}
#ifndef NDEBUG // Verifies that all the data in the file have been read:
  if(fileRead) {
    const loint real_file_size = buffer_string_list::get_file_size(FILE_INPUT_NAME);
    FILE *f_log = log_builder::get_log_file_pointer(n_threads, "pipe_parse_file", __FILE__, __LINE__, true);
    if(f_log) {
      fprintf(f_log, "-\tIn total %lld chars (of %lld possible) were read for file %s, giving a difference of %lld chars.\n", chars_in_file_read_first_run, real_file_size,
	      FILE_INPUT_NAME, (loint)((loint)chars_in_file_read_first_run - (loint)real_file_size)); // D
      fclose(f_log);
    }
    if(real_file_size != chars_in_file_read_first_run) {
      fprintf(stderr, "!!\treal_file_size(%d) != chars_in_file_read_first_run(%d) with a difference(%d). (This error is found at line %d in file %s. If seen, please contact the developer at oekseth@gmail.com.)\n",
	      (int)real_file_size, (int)chars_in_file_read_first_run, (int)(real_file_size-  chars_in_file_read_first_run),
	      __LINE__, __FILE__);
      assert(real_file_size == chars_in_file_read_first_run);
    }   
  }
#endif
}



void pipe_parse_file::assert_class(const bool print_info) {
  const static char *class_name = "pipe_parse_file";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  // Iinitial tests:
  char *test_file = "test_file.txt";
  uint size_line = 100, cnt_newlines = 4,  cnt_block = 2;
  char *test_block = new char[size_line*cnt_newlines+1];
  for(uint i = 0; i < (size_line*cnt_newlines); i++) test_block[i] = 'a';
  for(uint i =1; i < cnt_newlines; i++) {
    test_block[i*size_line-1] = '\n';
  }
  test_block[size_line*cnt_newlines-1] = '\n';

  test_block[size_line*cnt_newlines] = '\0';
  FILE *file = fopen(test_file, "w");
  uint cnt_chars = 0;
  if(file != NULL) {
    for(uint i = 0; i < cnt_block; i++) {
      fputs(test_block, file);
      cnt_chars +=strlen(test_block);
    }
    fclose(file);
  }

  log_builder_t *log = log_builder::init();//(NULL, '_', 0, 0, 0);
  //  log_builder_t *log = new log_builder[1];//(NULL, '_', 0, 0, 0);
  pipe_parse_file temp(test_file, log);
  
  string_section *section = temp.read_file_blocks();
  assert(cnt_chars == section->get_data_length()); // 4 lines excluding the newline
  //  assert(399 == section->get_data_length()); // 4 lines excluding the newline
  string_section::close(section);
  // // section->finalize();
  // free(section);

  // section = temp.read_file_blocks();

  // assert(399 == section->get_data_length()); // 4 lines excluding the newline
  // section->finalize();
  // free(section);

  // section = temp.read_file_blocks();

  // if(section != NULL) {
  //   assert(0 == section->get_data_length()); // 4 lines excluding the newline
  //   section->finalize();
  //   free(section);
  // }

  // if(log) {log->free_mem(), delete log, log = NULL;}
  // delete [] test_block;
  delete [] test_block; test_block = NULL;
  temp.close_file(1);
  log_builder::close(log);
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
  //  if(false)   print_constants();
}

