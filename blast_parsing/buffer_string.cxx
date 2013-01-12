#include "buffer_string.h"

//! Prints info about the object.
void buffer_string::print_info() {
  printf("class 'buffer_string':\t");
  if(first_buffer || size_first_buffer) printf("first_buffer(%llu)\t", size_first_buffer);
  if(second_buffer || size_second_buffer) printf("second_buffer(%llu)\t", size_second_buffer);
  printf("\n");
}

/**
   @param <buff> The data to copy into the first buffer.
   @param <size> The length of the data to copy into the first buffer.
**/
void buffer_string::copy_into_first_buffer(char *buff, loint size) {
  //  if(true) {
  assert(size_first_buffer == 0);
  first_buffer = new char[size+1]; 
  memcpy(first_buffer, buff, size);
  first_buffer[size] = '\n';
  size_first_buffer = size+1;
  //     assert(size_first_buffer == 0);
  //     first_buffer = new char[size]; 
  //     memcpy(first_buffer, buff, size);
  //     size_first_buffer = size;
  //   }
}

/**
   @brief Sets the second buffer with the params given.
   @remark Uses the given pointer, without any copying.
**/
void buffer_string::set_second_buffer(char *&buff, loint size) {
  second_buffer = buff, size_second_buffer = size;    
}


void buffer_string::print_segment(char *m_start, loint length) {
  char *start =m_start; 
  //printf("[\n");
  loint cnt = 0;
  while(cnt != length) putchar(*start++), cnt++;
  printf("\n");
  //printf("\n]\n");
}

void buffer_string::print_first_buffer() {if(first_buffer)print_segment(first_buffer, size_first_buffer);}
void buffer_string::print_second_buffer(){if(second_buffer)print_segment(second_buffer,size_second_buffer);}

//! Sets the input param as the merged result of the two internal buffers.
void buffer_string::get_string(loint &length, char *&str) {
  if(first_buffer == NULL) {length =size_second_buffer; str = second_buffer; size_second_buffer = 0;}
  else if(second_buffer == NULL) {length =size_first_buffer; str = first_buffer; size_first_buffer = 0;}
  else { // Merge:
    length = size_first_buffer + size_second_buffer+1;
    str = new char[length+1]; str[length] = '\0';
    memcpy(str, first_buffer, size_first_buffer);
    str[size_first_buffer] = '\n'; 
    memcpy(str+size_first_buffer+1, second_buffer, size_second_buffer);
    // TODO: A '\n' resides in them: verify that this is how it should bee!
    //    printf("new_string(%s) = '%s' + '%s'\n", str, first_buffer, second_buffer);
  }
}

/**
   @brief De-allocates the list of objects given as input, defined
   @param <tmp> The list of objects to de-allocate.
   @param <size_index> The number of objects provided to de-allocate.
**/
void buffer_string::close(buffer_string *&tmp, loint &size_index) {
  if(tmp) {
    for(uint i = 0; i < size_index; i++) tmp[i].free_mem(i);
    delete [] tmp; tmp  = NULL, size_index = 0;
  }
}

void buffer_string::free_mem(uint block_id) {
  if(false && ((NULL!=first_buffer) || (NULL!=second_buffer))) printf("frees block_id(%d)\n", block_id);
  free_mem();
}

void buffer_string::free_mem() {
  if(first_buffer && (size_first_buffer > 0)) {
    //    static uint deb_cnt = 0; printf("%d\t\tat line %d in buffer_string frees the first_buffer of length(%llu)\n", deb_cnt++, __LINE__, size_first_buffer);
    delete [] first_buffer, first_buffer = NULL, size_first_buffer = 0;
  } else {
    //    static uint deb_cnt = 0; printf("(not set)\t%d\t\tat line %d in buffer_string frees the first_buffer of length(%llu)\n", deb_cnt++, __LINE__, size_first_buffer);
  }
  if(second_buffer && (size_second_buffer > 0)) {
    //    static uint deb_cnt = 0; printf("%d\t\tat line %d in buffer_string frees the second_buffer of length(%llu)\n", deb_cnt++, __LINE__, size_second_buffer);
    delete [] second_buffer, second_buffer = NULL, size_second_buffer = 0;
  } else {
    //    static uint deb_cnt = 0; printf("(not set)\t%d\t\tat line %d in buffer_string frees the second_buffer of length(%llu)\n", deb_cnt++, __LINE__, size_second_buffer);
  }
}


void buffer_string::assert_class(const bool print_info) {
  const static char *class_name = "buffer_string";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  {
    buffer_string temp_0 = buffer_string();
    assert(!temp_0.has_data());
    char *str = new char[11]; memset(str, 'a', 10); str[10] = '\0';
    temp_0.copy_into_first_buffer(str, 10);
    assert(temp_0.has_data());
    loint len = 0;
    char *res = temp_0.get_string(len);
    assert(0==strncmp(str, res, 10));
    assert(len == (10+1)); len = 0;
    assert(str != NULL); // Data is only copied, and not removed
    temp_0.free_mem(0);
    if(str) delete [] str, str = NULL;
    if(res) delete [] res, res = NULL;
  }
    
  if(true) {
    buffer_string temp_0 = buffer_string();
    //    assert(!temp_0.has_data());
    char *str = new char[11]; memset(str, 'a', 10); str[10] = '\0';
    char *tst = new char[11]; memset(tst, 'a', 10); tst[10] = '\0';
    temp_0.set_second_buffer(str, 10);
    assert(temp_0.has_data());
    loint len = 0; char *res = temp_0.get_string(len);
    assert(0==strncmp(tst, res, 10));
    //    printf("length(%llu)\n", len);
    assert(len == 10); len = 0;
    temp_0.free_mem(0);
    //    assert(str == NULL); // Data is freed
    //    assert(!temp_0.has_data());
    if(tst) delete [] tst, tst = NULL;
    if(res) delete [] res, res = NULL;

  }
  {      
    buffer_string temp_0 = buffer_string();
    // Tests merging:
    char *a = new char[10]; memset(a, 'a', 10);
    char *b = new char[20]; memset(b, 'b', 20);
    temp_0.copy_into_first_buffer(a, 10);
    temp_0.set_second_buffer(b, 20);
    assert(temp_0.has_data());
    loint len = 0; char *res = temp_0.get_string(len);
    // assert(len == 30); len = 0;
    // assert(0 == strncmp(a, res, 10));
    assert(a != NULL);
    delete [] a;
    //    assert(0 == strncmp(b, res+10, 20));
    temp_0.free_mem(0);
    assert(!temp_0.has_data());
    if(res)     delete [] res; res = NULL;
      
  }
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
