#include "bucket_norm.h"
//! Deallocates the memory
void bucket_norm::free_mem(const uint taxon_length) {
  list_norm::close(arrNorm);
}

//! Initializes the class
bucket_norm *bucket_norm::init(list_norm_t *_arrNorm) {
  bucket_norm *b =  (bucket_norm*)malloc(sizeof(bucket_norm));
  b->arrNorm = _arrNorm;
  return b;
}

//! Constructs the class given the input params.
bucket_norm::bucket_norm(list_norm_t *_arrNorm) : arrNorm(_arrNorm) {}

//! Asserts the class given with test functions.
void bucket_norm::assert_class(const bool print_info) {
  const static char *class_name = "bucket_norm";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
#endif
  // Has nothing to test.
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
