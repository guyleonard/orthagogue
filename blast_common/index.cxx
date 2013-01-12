#include "index.h"
//! Prints information of the object
void index::print(uint i) {
  printf("\t[%u]\t .start = %u, .length=%u:\n", i, start, length);
}
//! Prints information of the object using the name given as param.
void index::print(char *name) {
  printf("\t[%s]\t .start = %u, .length=%u:\n", name,  start, length);
}
//! Prints information of the object.
void index::print(uint i, char **name) {
  if(name == NULL) print(i);
  else print(name[i]);
}
/**
   @brief The main test function for this class.
   @remarks This method includes:
   - Formalized tests for setting of starting position and end position.
   - Valgrind used verifying memory usage.
   - Examples of how this class may be used.
   @author Ole Kristian Ekseth (oekseth)
   @date 04.01.2012 by oekseth.
**/
void index::assert_class(const bool print_info) {
  const static char *class_name = "index";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  index temp = index();
  assert(0 == temp.get_length());
  assert(0 == temp.get_absolute_start_pos());
  assert(0 == temp.get_start_pos());
  temp.increment_length();
  assert(1 == temp.get_absolute_start_pos());
  temp.increment_length(19);
  assert(20 == temp.get_absolute_start_pos());

  // Tests building of the list of indexes
  const uint index_size = 10;
  index *list = init_list(index_size);
  for(uint i = 0; i < index_size; i++) {
    index temp = list[i];
    assert(0 == temp.get_length());
    assert(0 == temp.get_absolute_start_pos());
    assert(0 == temp.get_start_pos());
    temp.increment_length();
    assert(1 == temp.get_absolute_start_pos());
    temp.increment_length(19);
    assert(20 == temp.get_absolute_start_pos());
  }
  index::delete_list(list);
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
