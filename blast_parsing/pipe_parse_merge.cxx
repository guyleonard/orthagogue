#include "pipe_parse_merge.h"

list_file_parse<p_rel> *pipe_parse_merge::getListParseData(){
  if(parse_data_this != NULL)  return parse_data_this;
  else {
    log_builder::throw_error_and_abort("The list holding the set of protein comparisons were not set, probably due to erronous settings and/or input file", __LINE__, __FILE__, __FUNCTION__);
  }
  return NULL;
}

taxon_list_t *pipe_parse_merge::getListTaxon() {
  if(listTaxon != NULL) return listTaxon;
  else   log_builder::throw_error_and_abort("The list holding the taxa were not set, probably due to errounous settings and/or input file", __LINE__, __FILE__, __FUNCTION__);
  return NULL;
}

/**! Gets the overlap values into the taxa structure, for later use. Frees the container.
   @Changed: 30.12.2010 by oekseth
*/
void pipe_parse_merge::getInitTaxaArrOverlap(taxa_t *&listTaxa) {
#ifndef NDEBUG
  if(hashProtein && arrOverlap && taxon_length) {
    uint sum_after = 0;
    uint cnt_proteins = 0;
    for(int taxon =0; taxon< taxon_length; taxon++) {
      for(int prot = 0; prot <hashProtein[taxon].getLength(); prot++) {
	cnt_proteins++;
	assert(arrOverlap[taxon]);
	sum_after += arrOverlap[taxon][prot];
      }
    }

    assert(sum_after == debug_sum_inserted_overlap_values);  
  }
#endif
  for (int taxon = 0;taxon<taxon_length;taxon++) {
    listTaxa[taxon].set_unsafe_overlap_list(arrOverlap[taxon]);
  }
  free(arrOverlap); arrOverlap=NULL;
}

/**! Called from the outside before the second pipeline is run */
void pipe_parse_merge::initSecondRead(taxa_t *&_listTaxa, int _taxon_length, prot_list_t *&_hashProtein // The list of the hashed proteins in the collection
				      ) {
  next_block = 0;
  hashProtein = _hashProtein;
  listTaxa = _listTaxa;
  taxon_length = _taxon_length;
  FIRST_READ = false;
  arrOverlap = (overlap_t**)malloc(sizeof(overlap_t*)*taxon_length);
  for(int i = 0; i<taxon_length;i++) {
    arrOverlap[i] = (overlap_t*)malloc(sizeof(overlap_t)*hashProtein[i].getLength());
    memset(arrOverlap[i], 0, sizeof(overlap_t)*hashProtein[i].getLength());
  }
  if(!MAX_PARSE_BUFFER_SIZE) MAX_PARSE_BUFFER_SIZE = 0;
  parse_data_this = new list_file_parse<p_rel>(hashProtein, MAX_PARSE_BUFFER_SIZE, taxon_length, listTaxa, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);
  // TODO: include below when the new testfunction in list_file_parse is written.
  //  parse_data_this->set_max_buffer_size_based_on_file_properties(listTaxa, taxon_length);
  log_builder::test_memory_condition_and_if_not_abort((parse_data_this!= NULL), __LINE__, __FILE__, __FUNCTION__);
  if(listTaxon!= NULL) {
    taxon_list::close(listTaxon); 
  }
}


//! Cleans the memory for an outer overlap array
void pipe_parse_merge::free_overlap(overlap_t **arr) {
  if(arr != NULL) {
    for(int i = 0; i<taxon_length;i++) {
      free(arr[i]);
    } free(arr); arr = NULL;
  }
}

//! frees the memory for arrOverlap and the 'taxon_list'
void pipe_parse_merge::finalize_class() {
#ifndef NDEBUG
  if(hashProtein && arrOverlap && taxon_length) {
    uint sum_after = 0;
    for(int taxon =0; taxon< taxon_length; taxon++) {
      for(int prot = 0; prot <hashProtein[taxon].getLength(); prot++) {
	assert(arrOverlap[taxon]);
	sum_after += arrOverlap[taxon][prot];
      }
    }
    assert(sum_after == debug_sum_inserted_overlap_values);  
  }
#endif
  free(arrOverlap); arrOverlap = NULL;
  // TODO: When possibli validate that this is correct, and that not only thepointer should be removed (deleted)
  taxon_list::close(listTaxon);//  free(listTaxon); listTaxon = NULL;
}

void pipe_parse_merge::finalize_parse_blocks() {
  if(parseBlocks) {	parseBlocks->free_mem(), delete parseBlocks, parseBlocks = NULL;}
}


//! Builds a consistent list of the overlap scores
// @overlap_t **overlap: The input to be merges with the global list of overlaps
void pipe_parse_merge::merge_overlap(overlap_t **overl) {
  uint sum_inserted = 0;
  uint sum_before = 0;
  uint sum_total = 0;
  for(int taxon =0; taxon< taxon_length; taxon++) {
    for(int prot = 0; prot <hashProtein[taxon].getLength(); prot++) {
      sum_before += arrOverlap[taxon][prot];
      if(overl[taxon][prot] != 0) {
	//	arrOverlap[taxon][prot] = overl[taxon][prot];
	arrOverlap[taxon][prot] += overl[taxon][prot];
	sum_inserted+= overl[taxon][prot];
      }
      sum_total += arrOverlap[taxon][prot];
    }
  }
#ifndef NDEBUG
  const uint sum_expected = sum_before + sum_inserted;
  assert(sum_total == sum_expected);
  debug_sum_inserted_overlap_values += sum_inserted;
#endif
  free_overlap(overl);
}


/*overrride*/void* pipe_parse_merge::operator()(void* item) {
  const clock_t tstart = times(NULL);
  if(FIRST_READ) {
    parse_send_first_t *parse = static_cast<parse_send_first_t*>(item);
    if(parse != NULL) {
      if(parse->data_resides_in_mem) parseBlocks->update_last_file_pos_im_mem(parse->last_block_pos); 
      
      blocks_end.push(parse[0].last_block_pos);
      if(next_block == (int)parse[0].block_number) {
	next_block++;
	parseBlocks->insert(parse->last_block_pos);
	if(listTaxon!= NULL) {

	  listTaxon->merge(parse->getTaxonList());
	  parse->free_taxon_list_mem();
	  //	  free(parse->getTaxonList());	
	} else listTaxon = parse->getTaxonList();

	delete parse, parse=NULL;
      } else {
	lstFirstParse.push_back(parse[0]);
	free(parse);
	bool changed = true;

	// Iterates over the blocks
	while(changed) {
	  changed = false; // If this is chnaged, a new iteration through the block is needed
	  queue<parse_send_first_t> del_ele;
	  for(list<parse_send_first_t>::iterator it = lstFirstParse.begin(); it != lstFirstParse.end();it++) {
	    if(next_block == (int)it->block_number) {
	      parseBlocks->insert(parse[0].last_block_pos);
	      if(listTaxon!= NULL) {
		listTaxon[0].merge(it->getTaxonList());      
		free(parse->getTaxonList());	 
	      } else listTaxon = parse->getTaxonList();
	      del_ele.push(*it); // The data is extracted, deletes the data
	      next_block++;
	      changed = true;
	    }
	  }
	  while(!del_ele.empty()) { // Those data who are extracted is to be deleted
	    free(del_ele.front().t_list);
	    del_ele.pop();
	  }
	}
      }
    } 
  } else { // The SECOND read:
    parse_send_t *parse = static_cast<parse_send_t*>(item);
    if(parse && parse->has_data()) {
      if(parse->max_sim_value > max_sim_score) max_sim_score = parse->max_sim_value;      
      merge_overlap(parse->arrOverlap);
      struct protein_relation protrel  = parse->prot;
      if(parse->parseData) {
	assert(parse_data_this);
	list_file_parse<p_rel> *parse_block_c = parse->parseData;
#ifndef NDEBUG
	const uint tot_elements = parse_data_this->getTotalLengthOfData() + parse_block_c->getTotalLengthOfData();;
#endif
	//	uint cnt_elements_in_all_relation_lists = 0; // FIXME: make  this a globally access bile object.
	parse_data_this->merge_data(parse_block_c, protrel, cnt_elements_in_all_relation_lists);
#ifndef NDEBUG
	assert(!parse_block_c); // ARgument should be deleted at this point.
	assert(tot_elements == (parse_data_this->getTotalLengthOfData())); // Asserts that we-we got the correct result.
	//	const uint main_cnt = parse_data_this->getTotalLengthOfData();
#endif
      } 
    } else {
      next_block++;
    }
    free(parse);
  }
  logid id = read_first;
  if(!FIRST_READ) id = read_second;
  const clock_t tend = times(NULL); log->append_measurement(id, 2, tstart, tend);
  return NULL;
}

pipe_parse_merge::pipe_parse_merge(uint disk_buffer_size, int _taxon_length, log_builder_t *_log, lint _MAX_PARSE_BUFFER_SIZE, taxa_t *_listTaxa, char *_binary_loc, uint _CPU_TOT) :
  tbb::filter(  /*seriel=*/ true), USE_BEST_BLAST_PAIR_SCORE(false), listTaxa(_listTaxa), FILE_BINARY_LOCATION(_binary_loc), 
  MAX_PARSE_BUFFER_SIZE(_MAX_PARSE_BUFFER_SIZE), log(_log)
  , taxon_length(_taxon_length), parse_data_this(NULL), listTaxon(NULL), next_block(0),
  FIRST_READ(true), CPU_TOT(_CPU_TOT), cnt_elements_in_all_relation_lists(0), max_sim_score(0.0), debug_sum_inserted_overlap_values(0), parseBlocks(NULL),
  arrOverlap(NULL)  {
  parseBlocks = new parse_read_blocks(disk_buffer_size);
}
 
#ifdef assert_code
void pipe_parse_merge::assert_private_parts() {
  ;
}
#endif

void pipe_parse_merge::assert_class(const bool print_info) {
  const static char *class_name = "pipe_parse_merge";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  log_builder_t *log = log_builder::init();
  const lint MAX_PARSE_BUFFER_SIZE = 10000;
  taxa_t *listTaxa = NULL;
  char *FILE_BINARY_LOCATION = NULL;
  pipe_parse_merge temp(1024, 2, log, MAX_PARSE_BUFFER_SIZE, listTaxa, FILE_BINARY_LOCATION, 1);
  
  temp.finalize_parse_blocks();
  temp.finalize_class(); 
  //  if(false)   print_constants();
  log_builder::close(log);
#endif
  // Iinitial tests:
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

