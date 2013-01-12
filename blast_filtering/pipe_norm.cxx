#include "pipe_norm.h"
/**! Updates the global list of the normalization values
    @ProgramPoint: in 'void operator'
*/
void pipe_norm::mergeArrNorm(list_norm_t *norm) {
  if(arrNorm) arrNorm->merge_basis(*norm);
  else fprintf(stderr, "!!\tObject of type list_norm was not set. Error found at line %d in file %s for method %s\n", __LINE__, __FILE__, __FUNCTION__);
}

/*override*/void * pipe_norm::operator()(void* item) {
  const clock_t tstart = clock_t();
  bucket_norm *bucket = static_cast<bucket_norm*>(item);
  if(bucket) {
    if(!USE_EVERYREL_AS_ARRNORM_BASIS) mergeArrNorm(bucket->arrNorm);
    bucket->free_mem(taxon_length);
    free(bucket);
  }
  const clock_t tend= clock_t(); log->append_measurement(filter_co_orthologs, 2, tstart, tend);
  return NULL;
}

pipe_norm::pipe_norm(uint _taxon_length, const bool _use_everyrel_as_arrnorm_basis,
	    const pipe_t type, log_builder_t *_log)	      
    :
  tbb::filter(true), // serial
    log(_log),
    taxon_length(_taxon_length),
    USE_EVERYREL_AS_ARRNORM_BASIS(_use_everyrel_as_arrnorm_basis),
    PIPE_TYPE(type),
    arrNorm(NULL) 
      {     }

#ifdef assert_code
void pipe_norm::assert_private_parts() {
  ;
}
#endif

void pipe_norm::assert_class(const bool print_info) {
  const static char *class_name = "pipe_norm";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  log_builder_t *log = new log_builder[1];//(NULL, '_', 0, 0, 0);
  pipe_norm temp(2, false, ORTH, log);
  //  temp.free_mem();
  //  if(false)   print_constants();
  if(log) {log->free_mem(), delete log, log = NULL;}
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
