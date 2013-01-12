#include "pipe_bucket.h"
/**
   @brief Produces buckets of numbered items to parse. Goal is to make the threads in the rest of the pipe effective:
   @todo Should be considered removed, replaced by preestimated blocks used by each
   of the consquative threads.
   @date 21.12.2010 by oekseth
*/
struct taxon_pair* pipe_bucket::init_taxon_p(const bool inpa_ops, uint taxon_start, uint taxon_end,
					     uint taxon_length, uint n_threads)// , int n_threads)
{
  struct taxon_pair *listPair; // Holds the blocks to be used for the parsing
  const clock_t tstart = clock_t();
  listPair =  taxon_pair::init_taxon_pair(taxon_start, taxon_end,
					  taxon_length, n_threads,
					  inpa_ops, // uses the orthologs as a basis if this variable is set to false
					  listParseData, listStructData,
					  listTaxa
					  );
  const clock_t tend = clock_t();
  log->append_measurement(filter_orthologs_init, 0, tstart, tend);
  return listPair;  
}  

void* pipe_bucket::operator() (void* item) {
  if(false) printf("i=%ld", (long int)item); // Hiding an unsed item from the compiler
  const clock_t tstart = clock_t();
  if(listPair) {
    taxon_pair *buck = taxon_pair::get_taxon_pairs(listPair, list_pair_pos);      
    const clock_t tend = clock_t();
    logid_t id = filter_orthologs;
    if(is_inpa) id = filter_inparalogs;
    log->append_measurement(id, 0, tstart, tend);
    return buck;
  } else {
    // No data were set for this taxon.
    return NULL;
  }
}

//! The constructor.
pipe_bucket::pipe_bucket(const bool inpa_ops, uint taxon_start,
			 uint taxon_end,
			 uint taxon_length, log_builder_t *_log,
			 list_file_parse_t *&_listParseData,
			 list_file_struct_t *&_listStructData,
			 taxa_t *_listTaxa, uint n_threads
			 ) : 
  tbb::filter(true),   listParseData(_listParseData),   listStructData(_listStructData), listTaxa(_listTaxa), log(_log), is_inpa(inpa_ops) // serial
{ 
  list_pair_pos = 0;
  listPair = init_taxon_p(inpa_ops, taxon_start, taxon_end, taxon_length, n_threads);
}

#ifdef assert_code
void pipe_bucket::assert_private_parts() {;}
#endif
void pipe_bucket::assert_class(const bool print_info) {
  const static char *class_name = "pipe_bucket";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  log_builder_t *log = new log_builder[1];//(NULL, '_', 0, 0, 0);
  list_file_parse_t *listParseData = NULL;
  list_file_struct_t *listStructData = NULL;
  taxa_t *listTaxa = NULL;
  const uint n_threads = 1;
  pipe_bucket temp(false, 0, 2, 3, log, listParseData, listStructData, listTaxa, n_threads);
  temp.free_mem();
  if(log) {log->free_mem(), delete log, log = NULL;}
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

