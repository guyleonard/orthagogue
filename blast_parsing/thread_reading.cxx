#include "thread_reading.h"


void thread_reading::print_block() {
  printf("-\t block_index_start(%u), block_index_end(%u), current_index(%u)\n", block_index_start, block_index_end, current_index);
}
//! Returns the total length of the data to be read by this object.
loint thread_reading::get_total_length(loint *arr) {
  return get_block_length(arr, block_index_start, block_index_end-1);
}

//! Returns true if the file is read
bool thread_reading::file_is_read(uint my_id, bool print_data_) {
  if(current_index == block_index_end) return true;
  else {
    //      assert(!print_data_);
    if(print_data_) {printf("in thread_reading.cxx at line %d, my_id[%u] the data: ", __LINE__, my_id); print_block();}
    return false;
  }
}

//! Returns the length of the block given: assumes the indexes ar in the range of the arr given.
loint thread_reading::get_block_length(loint *arr, uint start_index, uint end_index) {
  const loint file_pos = get_file_pos(arr, start_index);
  return arr[end_index] - file_pos; // new_pos - old_pos
}

//! Returns the starting position in the file for the given index:
loint thread_reading::get_file_pos(loint *arr, uint index) {
  loint file_pos = 0;
  if(index > 0) file_pos = arr[index-1];
  return file_pos;
}

//! Allocates and initializes the read-buffer to be sent downward the pipe:
void thread_reading::allocate(char *&buffer_main, char *&buffer_end, const loint buffer_in_mem_size, char *caller_file, int caller_line) {
  try {buffer_main = new char[buffer_in_mem_size+1];} 
  catch (std::exception& ba) {
    if(!log_builder::catch_memory_exeception(buffer_in_mem_size+1, __FUNCTION__, __FILE__, __LINE__)) {
      fprintf(stderr, 
	      "!!\t An interesting error was discovered: %s.\n"
	      "-\t This function was called by %s:%d\n"
	      "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
	      "Error generated at [%s]:%s:%d\n",  ba.what(), 
	      caller_file, caller_line,
	      __FUNCTION__, __FILE__, __LINE__);
    } else {
      fprintf(stderr, 
	      "!!\t An interesting error was discovered: %s.\n"
	      "-\t To fix this error, either (if possible) increase the accessible chunk of memory (available for orthAgogue), or decrese the disk buffer size (\"-dbs\" parameter.)\n"
	      "-\t This function was called by %s:%d\n"
	      "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
	      "Error generated at [%s]:%s:%d\n",  ba.what(), 
	      caller_file, caller_line,
	      __FUNCTION__, __FILE__, __LINE__);
    }

  }
  assert(buffer_main);
  //  buffer_main = new char[buffer_in_mem_size+1];//new char[buffer_in_mem_size + 1];
  buffer_end = buffer_main + buffer_in_mem_size;
  memset(buffer_main, '\0', buffer_in_mem_size+1);// Filling in the buffer:
}


/**! Returns the length of the block, incrementing the poiner for the next read
   @Return: '0' if no more blocks to read
*/
loint thread_reading::getNextReadLength(loint *arr) {
  if(current_index < block_index_end) {
    const loint file_pos = get_file_pos(arr, current_index);
    loint length = 0;
    do { // Makes the block as large as possible, but not exceeding the maximal limit
      length = arr[current_index] - file_pos; // new_pos - old_pos
      current_index++;
    } while((current_index < block_index_end) && (length < (loint)disk_buffer_size)); 

    return length;
  }
  return 0;
}

/**
   @brief Reads the file using the position given as input.
   @param <file_position_array>  The array where the position data is located.
   @param <buffer_main> The reserved char array to read the data into.
   @param <buffer_in_mem_size> The variable to get the length of the chars inserted into the char parameter.
   @return true if the there are no more chars to read.
   @author Ole Kristian Ekseth (oekseth)
**/
bool thread_reading::read_file(loint *file_position_array, char *buffer_main, loint &buffer_in_mem_size) {
  if(!file_input) {
    file_input = fopen(file_name, "rb");
    if(file_input) {
      // Below only neccessary to move the pointer if its not the start of the file:
      loint pos_in_file = global_file_reading_start;
      if(block_index_start > 0) pos_in_file = (loint)file_position_array[block_index_start-1];      
      fseek(file_input, pos_in_file, SEEK_SET); // positions the file
      buffer_in_mem_size = buffer_in_mem_size - pos_in_file; // Updates the length using the position offset.
    } else fprintf(stderr, "!!\tNot possible open file %s at line %d in class thread_reading\n", file_name, __LINE__);
  }
  const size_t chrs_read = fread(buffer_main, sizeof(char), buffer_in_mem_size, file_input);
  if(chrs_read != (size_t)buffer_in_mem_size) {
    memset(buffer_main+chrs_read, '\0', (size_t)buffer_in_mem_size-chrs_read);
    buffer_in_mem_size = chrs_read; // Updates the size
    return true;
  } else return false; 

}

  /**
     @brief Reads a block from the input file using positions in input.
     @param <my_id> The identifier specifying the if of the file pointer to use.
     @param <file_position_array>  The array where the position data is located.
     @return the section of the data using the parse-blocks as a basis:
     @remarks Is thread-safe. Steps taken are:
     - Get number of chars to read (length), summing up several 'lengths" if they in total constitute a sum of bytes less that the variable disk_buffer_size.
     - On the first call moves the file pointer to the point were this object starts the read. (For an example of this usage, see class pipe_parse_merge in combination with class buffer_string_list.
     @author Ole Kristian Ekseth (oekseth)
  **/
string_section *thread_reading::get_string_section(loint *file_position_array) {  
  loint buffer_in_mem_size = getNextReadLength(file_position_array);
  if(buffer_in_mem_size > 1) {
    if((file_input == NULL) || !feof(file_input)) {
      char *buffer_main = NULL, *buffer_end = NULL; // Holding the grabbed string.
      allocate(buffer_main, buffer_end, buffer_in_mem_size+100, __FILE__, __LINE__); // Uses offset of '100' to ensure the holw line comes in.
      memset(buffer_main, '\0', buffer_in_mem_size+100); // TODO: consider removing this- and the 100-term above, and see if things workds as they should!
      bool at_end_of_file = read_file(file_position_array, buffer_main, buffer_in_mem_size);
      string_section_t *section = new string_section_t(buffer_main, buffer_main + buffer_in_mem_size, 
						       /*prev_index=*/0, 
						       /*not important*/0, at_end_of_file);   
      return section;
    } 
  }
  return NULL;
}

string_section *thread_reading::get_string_section(uint my_id, loint *file_position_array) {
  if(false) printf("my_id=%u\n", my_id); // For use in debugging if needed.
  return get_string_section(file_position_array);
}
//! De-allocates memory for the collection of objects, given as input.
void thread_reading::free_mem(thread_reading *list, uint size) {
  for(uint i = 0; i < size; i++) list[i].free_mem();
}

//! The constructor.
thread_reading::thread_reading(uint _disk_buffer_size, char *_file_name, uint _block_index_start, uint _block_index_end, loint _global) :
  disk_buffer_size(_disk_buffer_size), file_name(_file_name), file_input(NULL), block_index_start(_block_index_start), block_index_end(_block_index_end), current_index(block_index_start), global_file_reading_start(_global)
{}

//! The constructor.
thread_reading::thread_reading() :   disk_buffer_size(1024*1024), file_name(NULL), file_input(NULL), block_index_start(0), block_index_end(0), current_index(0), global_file_reading_start(0)
{}

void thread_reading::assert_class(const bool print_info) {
  const static char *class_name = "thread_reading";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  // Iinitial tests:
  //  uint length = 10; 
  thread_reading temp = thread_reading();
  // TODO: Tests of this class constitures of those usage-cases found in other classes. Consider wirting seperate tests for this class.
  /*
    for(uint i = 0; i < length; i++) {
    temp.insert(i+1); // Sets the psostion in the file: must be greater than one
    }
    for(uint i = 0; i < length; i++) {
    const uint curr = (uint)temp.getNextReadLength();  
    assert((uint)1 == curr);
    }
  */
  temp.free_mem();
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

