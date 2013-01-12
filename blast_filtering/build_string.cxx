#include "build_string.h"

//! Prints the buffer stored in this object.
void build_string::print() {
  printf("\n|");
  fputs(logical_start, stdout);
  printf("\n|");
}

//! Writes the string to the file, and clears the memory
void build_string::write_data_to_file(FILE  *out_file, const bool PIPE) {
  if(size_blast_all > 0 && logical_start != NULL) {
    //    printf("Writes %u data to row:\n", size_blast_all);
    FILE *file_out;
    if(PIPE) file_out = stdout;
    else file_out = out_file;
    fputs(logical_start, file_out);
    tbb::tbb_allocator<char>().deallocate(logical_start, size_blast_all);
    size_blast_all = 0;
    logical_start = NULL;
    logical_end   = NULL;
  }
}

// Sets the header
void build_string::setHeader(uint world_index_in) {
  uint strl_len = 1;
  if(world_index_in > 0)
    strl_len = (int)log10(world_index_in)+1;
  resize_if_to_large(strl_len);
  sprintf(logical_end, "%u ", world_index_in);
  logical_end += strl_len + 1;
}


//! Frees the memory allocated for this object
void build_string::free_mem() {
  tbb::tbb_allocator<char>().deallocate(logical_start, size_blast_all);
  size_blast_all = 0;
  logical_start = NULL;
  logical_end   = NULL;
}


//! Called when no data shall be printed to the stream
void build_string::avoidPrinting() {logical_end = NULL;}


//! Return true if data shall be printed
bool build_string::hasData() {return (logical_end!=NULL);}

//! Pointer to beginning of the sequence
char* build_string::begin() {return (char*)(logical_start);}

//! Appends the chars to the buffer
void build_string::append(char* first, char* last) {
  memcpy(logical_end, first, last-first);
  logical_end += (last-first); // Adding the length
}

/**
   @Brief adds the content of the input to the string located in this object.
   @param <start> Pointer to the first char in the sequence of chars.
   @param <end> Pointer to the last char in the sequence of chars.
*/
void build_string::copy(char *start, char *end) {
  const mem_loc length_argument = resize_if_to_large(end-start); 
  memcpy(logical_end, start, length_argument);
  logical_end += (length_argument);
}

//! adds the content of the input to the string located in this object.
void build_string::copy_struct(build_string row) {
  memcpy(logical_end, row.logical_start, row.size());
  logical_end += row.size();
}

/**! Increases the size of the buffer if it does not fit
   @ProgramPoint: At insertion of a new part of the char array
   @arg 'length_argument': the length of the input
   @returns: the size of the argument to be copied
*/ 
mem_loc build_string::resize_if_to_large(mem_loc length_argument) {
  const mem_loc length_curr = logical_end-logical_start; // the length of the current array used
  const mem_loc length_new = length_argument + length_curr; 
  if(length_new >= size_blast_all) {
    const mem_loc length_new_temp = length_new + 400;
    char *temp = (char*)tbb::tbb_allocator<char>().allocate(sizeof(char)*length_new_temp);
    std::memset(temp, ' ', length_new_temp);// the string with empty spaces
    memcpy(temp, logical_start, size_blast_all);
    tbb::tbb_allocator<char>().deallocate(logical_start, size_blast_all);
    logical_start = temp;
    logical_end = temp + length_curr;
    size_blast_all = length_new_temp;
  }
  return length_argument;
}




//! Allocates a list of 'this' class, but do not initialize it.
build_string *build_string::allocate_class(uint size) {
  return (build_string*)malloc(sizeof(build_string)*size);
}
//! Appends the line-end characters specified for the 'blast' file format
void build_string::end_blast_line() {
  if(size_blast_all > 0) {
    *logical_end = '$',  logical_end++;
    *logical_end = '\n', logical_end++;
  }
}

//! Appends the line end ('\0' character)
void build_string::finalize() {
  if(size_blast_all > 0) {
    *logical_end = '\0';
    logical_end++;
  }
}


build_string::build_string(uint _size_blast_all) :
  last_index(0), size_blast_all(_size_blast_all)
{
  if(size_blast_all > 0) {
    logical_end = (char*)tbb::tbb_allocator<char>().allocate(sizeof(char)*size_blast_all);
    std::memset(logical_end, ' ', size_blast_all);// the string with empty spaces
    logical_start = logical_end;
  } else {
    logical_end = NULL;
    logical_start = NULL;
  }
}

#ifdef assert_code

void build_string::assert_private_parts() {
  ;
}
void build_string::assert_class(const bool print_info) {
  const static char *class_name = "build_string";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
  build_string temp(2);
  temp.free_mem();
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
#endif
