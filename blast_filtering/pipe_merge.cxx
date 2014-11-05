#include "pipe_merge.h" 
/**
   @brief When the list_file_struct object is returned, it's the end-of-life for pipe_merge, and a log file is generated.
   @return the processed *list_file_struct object.
**/
list_file_struct_t *pipe_merge::getFileStruct(int n_threads, bool MODE_PAIRWISE_OUTPUT_ABC, bool MODE_PAIRWISE_OUTPUT_MCL, char *FILE_BINARY_LOCATION) {
#ifdef USE_MPI
  if(!fileStruct) {
    //! If not set, allocates an empty object:
    fileStruct = list_file_struct_t::init_class(taxon_length, listTaxa, MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, FILE_BINARY_LOCATION);
  }
#else
  if(!fileStruct) {
    fprintf(stderr, "!!\tThe object holding all the not-filtered-away pairs is empty, implying that no data remains after the filtering process.\n"
	    "-\t If this message seems strange, please install the software in debug-mode: detailed summaries are found in the \"report_orthAgogue\" folder.\n"
	    "-\t If you expected true results from running orthAgogue, please contact the developer at [oekseth@gmail.com].\n"
	    "This message was generated at line %d in file %s, located in method %s.\n", __LINE__, __FILE__, __FUNCTION__);
    assert(fileStruct);
  }
#endif
#ifndef NDEBUG  
  const loint size_buffer = fileStruct->getTotalLengthOfData();
  if(PIPE_TYPE == ORTH) {
    FILE *f_log = log_builder::get_log_file_pointer(n_threads, "pipe_process_merge", __FILE__, __LINE__, true);
    if(f_log) {      
      fprintf(f_log, "--\t The ortholog operation resulted in %lld pairs.\n", size_buffer);
      fclose(f_log);
    }
  } else {
    FILE *f_log = log_builder::get_log_file_pointer(n_threads, "pipe_process_merge", __FILE__, __LINE__, false);
    if(f_log) {      
      fprintf(f_log, "--\t The inparalog operation resulted in %lld pairs (started with %lld pairs).\n", size_buffer-debug_size_fileStruct_at_init, debug_size_fileStruct_at_init);
      fclose(f_log);
    }
  }  
#endif
  return fileStruct;
}

/**! Updates the global list of the normalization values
   @ProgramPoint: in 'void operator'
*/
void pipe_merge::mergeArrNorm(list_norm_t *norm) {
  if(arrNorm) {arrNorm->merge_basis(*norm);}
  else {
    fprintf(stderr, "!!\t'list_norm' object not set at line %d in file %s for method %s. If error not understood, send an email to oekseth@gmail.com with infomration about the observed behaviour.\n", __LINE__, __FILE__, __FUNCTION__);
  }
}

/*override*/void * pipe_merge::operator()(void* item) {
  const clock_t tstart = times(NULL);
  bucket_pipe_binary *bucket = static_cast<bucket_pipe_binary*>(item);
  //printf("(merge) at [%s]:%s:%d\n",  __FUNCTION__, __FILE__, __LINE__);
  if(bucket != NULL) {
    //printf("# (merge) at [%s]:%s:%d\n",  __FUNCTION__, __FILE__, __LINE__);
    if(bucket->isNotEmpty()) {
      if(fileStruct) {
	const loint cnt_this_before = fileStruct->getTotalLengthOfData();
	const loint cnt_arg_before = bucket->structData->getTotalLengthOfData();
	fileStruct->merge(bucket->structData);
	bucket->free_structdata();
	const loint cnt_this_after = fileStruct->getTotalLengthOfData();
	//printf("(merge)\t has %u elements in the structre after the insertion-buncj, at [%s]:%s:%d\n", (uint)cnt_this_after, __FUNCTION__, __FILE__, __LINE__);
	if(!(cnt_this_after == (cnt_this_before + cnt_arg_before))) {
	  fprintf(stderr, "!!\tHad in this %lld elements and %lld in argument before parsing, while %lld pairs after, giving a difference of %lld pairs. (This error is found at line %d in %s and PIPE_TYPE==ORTH(%d). Please contact oekseth@gmail.com if this message is seen.)\n",
		  cnt_this_before, cnt_arg_before, cnt_this_after, ((cnt_this_before + cnt_arg_before)-cnt_this_after), __LINE__, __FILE__, (PIPE_TYPE == ORTH));
	  assert(cnt_this_after == (cnt_this_before + cnt_arg_before));
	}
      } else {
	fileStruct = bucket->structData;
	//printf("# (merge) cnt=%u, at [%s]:%s:%d\n",  (uint)fileStruct->getTotalLengthOfData(),  __FUNCTION__, __FILE__, __LINE__);
      }
      if(arrNorm != NULL) {
	mergeArrNorm(bucket->arrNorm);
	list_norm::close(bucket->arrNorm);
      } else {
	arrNorm = bucket->arrNorm;
      }
    } else bucket->free_mem(taxon_length); 
    free(bucket); bucket = NULL;
  }
  logid_t id = filter_orthologs;
  if(PIPE_TYPE == INPA) id = filter_inparalogs;
  const clock_t tend = times(NULL);
  log->append_measurement(id, 2, tstart, tend);
  return NULL;
}

pipe_merge::pipe_merge(uint _taxon_length,
		       const bool _use_everyrel_as_arrnorm_basis,
		       const pipe_t type, log_builder_t *_log, list_file_struct_t *&_listStructData,
		       taxa_t *_listTaxa
		       )
  :
  tbb::filter(true), // serial
  listTaxa(_listTaxa), log(_log),
  taxon_length(_taxon_length),
  USE_EVERYREL_AS_ARRNORM_BASIS(_use_everyrel_as_arrnorm_basis),
  PIPE_TYPE(type),
  arrNorm(NULL),
  fileStruct(_listStructData), debug_size_fileStruct_at_init(0)
		      //,  listStructData(_listStructData)
{
  if(PIPE_TYPE==ORTH) {
    if(listTaxa != NULL) {
      for(uint i = 0; i< taxon_length; i++) listTaxa[i].initArrInpaLimit(/*init only if empty=*/true); // allocates memory for the limits to be created for the inparalogs
    } else {fprintf(stderr, "!!\tlistTaxa in class pipe_merge at line %d not set.\n", __LINE__);}
  }
  if(fileStruct) debug_size_fileStruct_at_init = fileStruct->getTotalLengthOfData();
}    
#ifdef assert_code
void pipe_merge::assert_private_parts() {
  ;
}
#endif
void pipe_merge::assert_class(const bool print_info) {
  const static char *class_name = "pipe_merge";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  log_builder_t *log = new log_builder[1];//(NULL, '_', 0, 0, 0);
  list_file_struct_t *listStructData = NULL;
  const uint taxon_length = 2;
  //  int *listTaxa = new int;
  taxa_t *listTaxa = new taxa[taxon_length]();
  pipe_merge temp(taxon_length, false, ORTH, log, listStructData, listTaxa);
  //  temp.free_mem();
  //  if(false)   print_constants();
  delete [] listTaxa; listTaxa = NULL;
  if(log) {log->free_mem(), delete log, log = NULL;}
#endif

  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

