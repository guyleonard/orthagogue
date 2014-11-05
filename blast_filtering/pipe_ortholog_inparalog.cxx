#include "pipe_ortholog_inparalog.h"


/**
   Returns the similarity score: if zero reciprocal vector is found, returns '0.0'
   @param  'buff_in_out': The buffer containg the vector data
   @param  'protein_out': The relative index of the outer protein
   @param  'protein_in': The relative index of the inner protein
   @param  'o_v_right': The overlap value to be set:
   @return  the right-left similarity score: if relation not found, returns '0.0'
*/
float pipe_ortholog_inparalog::getOuterSimilarity(basket_parse *basket,
						  uint taxon_out, uint taxon_in,
						  uint protein_out, uint protein_in,
						  overlap_t &overlap_right_left, overlap_t &recip_right_left_right)
{
  assert(basket);
  assert(listParseData);
  // Reads the in_out relations:
  if(listParseData->has_data(taxon_out, taxon_in, protein_out)){
    // There exists some relations between 'out' and 'inn'
    // The settings for this inner (left) protein:
    const uint length_out_in = listParseData->getProteinLength(taxon_out, taxon_in, protein_out);
    const uint start_out_in = listParseData->getProteinIndexStart(taxon_out, taxon_in, protein_out);
    for(uint ref_index_in = 0; ref_index_in<length_out_in; ref_index_in++) {
      if(basket->has_data(start_out_in, ref_index_in)) { // Tests if data is set for this protein:
	if(basket->get_index_out(start_out_in, ref_index_in) == protein_in) { // A to-way match is found
	  overlap_right_left= (overlap_t)basket->get_overlap_in(start_out_in, ref_index_in); // the overlap
	  recip_right_left_right= (overlap_t)basket->get_overlap_out(start_out_in, ref_index_in); //overlap
	  const float sim_score = basket->getSimScore(start_out_in, ref_index_in, max_input_value);
	  return sim_score;
	} 
      }
    }
  }
  return 0.0; // the reciprocal relation was not found
}


/**! Takes the average, and inserts it (in accordanse with the settings
   @ProgramPoint: If the 'aboveOverlapLimit' holds true, this function is called
   @Return: The averged value (to use as score for the proteins
   // TODO: Change made by oekseth 04. mars 2011 adding the 'INARAPLOG_OPERATION' lause, stating that inparalog operations are not to be tested when ont an Inparalog operation: needs to be verified, but no haste :)
   */
float pipe_ortholog_inparalog::insert_relation(uint my_id, uint taxon_in, uint taxon_out,
					       uint protein_in, uint protein_out,
					       float sim_in_out, float sim_out_in)
{
  const float sim_avg = (sim_in_out + sim_out_in)/2;
  if(!INPARALOG_OPERATION || (listTaxa[taxon_in].aboveInparalogLimit(/*taxon_in, */protein_in, sim_avg, INPARALOG_OPERATION)
			      && listTaxa[taxon_out].aboveInparalogLimit(/*taxon_out,*/ protein_out, sim_avg, INPARALOG_OPERATION))) {
    if(sim_avg > MIN_SIMILARITY_LIMIT) { // Store elements
      //printf("\t inserts relation, at [%s]:%s:%d\n", __FUNCTION__, __FILE__, __LINE__);
      l_fileStruct[my_id]->insert_relation(taxon_in, taxon_out, protein_in, protein_out, sim_avg);
      if(INPARALOG_OPERATION && !USE_EVERYREL_AS_ARRNORM_BASIS) { // Inserts the normative basis into the list of basis
	assert(l_arrNorm[my_id]);
	l_arrNorm[my_id]->insert(taxon_in, taxon_out, sim_avg, protein_in, protein_out, listTaxa);
      }
      return sim_avg;
    } else {
      if(DEBUG_PRINT_DISCARDED_PAIRS) printf("(discarded)\tpair (%s, %s) with averaged_similarity(%f) < similarity_limit(%f)\n", listTaxa[taxon_in].getArrKey(protein_in), listTaxa[taxon_out].getArrKey(protein_out), sim_avg, MIN_SIMILARITY_LIMIT);
    }
  } else {
    if(DEBUG_PRINT_DISCARDED_PAIRS) {
      if(listTaxa[taxon_in].arrInpaLimit[protein_in] >= sim_avg) 
	printf("(discarded)\tpair (%s, %s) with averaged_similarity(%f) < the leftmost_ortholog_limit(%f)\n", listTaxa[taxon_in].getArrKey(protein_in), listTaxa[taxon_out].getArrKey(protein_out), sim_avg, listTaxa[taxon_in].arrInpaLimit[protein_in]);
      else  {
	if(listTaxa[taxon_out].arrInpaLimit[protein_out] >= sim_avg) 
	  printf("(discarded)\tpair (%s, %s) with averaged_similarity(%f) < the rightmost_ortholog_limit(%f)\n", listTaxa[taxon_in].getArrKey(protein_in), listTaxa[taxon_out].getArrKey(protein_out), sim_avg, listTaxa[taxon_out].arrInpaLimit[protein_out]);
	else 	printf("(discarded)\tpair (%s, %s) with averaged_similarity(%f) ...due to other unknown reason (the data for the left- and right ortholog limits are: %f, %f)\n", listTaxa[taxon_in].getArrKey(protein_in), listTaxa[taxon_out].getArrKey(protein_out), sim_avg, listTaxa[taxon_in].arrInpaLimit[protein_in], listTaxa[taxon_out].arrInpaLimit[protein_out]);
      }
    }
  }
  return 0.0;
}

/**! Find the relations for this protein towards other proteins belonging to an outer taxon
   @ProgramPoint: in the 'for loop ' in the 'void* operator'
   @param 'arrOrtho' Holding the data to insert into
*/
void pipe_ortholog_inparalog::getProteinRelations(const uint my_id, basket_parse *basket_inn, basket_parse *basket_out, uint taxon_out, uint taxon_in, const uint real_protein_in/*, uint **arrOrthom*/) {
  // Base assumptions:
  assert(listTaxa);
  assert(listParseData);
  assert(basket_inn);
  assert(basket_out);
  // The settings for this inner (left) protein:
  const uint protein_in = listTaxa[taxon_in].getLocalIndex((uint)real_protein_in);
  // Containers for the maximum score:
  float this_sim_max = 0.0;
  stack<uint> max_sim_key;
  const bool debug = false;
  // Reads the in_out relations:
  if(listParseData->has_data(taxon_in, taxon_out, protein_in)) {// The binary_buffer is not empty
    // There exists some relations between 'out' and 'inn'
    const uint length_in_out = listParseData->getProteinLength(taxon_in, taxon_out, protein_in);
    const uint start_in_out = listParseData->getProteinIndexStart(taxon_in, taxon_out, protein_in);
    for(uint ref_index_out = 0; ref_index_out<length_in_out; ref_index_out++) {		    
      // 'protein_out' is the name of the id found in the list, ie, the (in-out) pair we are interested in:
      // This implies that only if it's found, we are interested in it:
      if(basket_inn->has_data(start_in_out, ref_index_out)) {
	const uint protein_out = basket_inn->get_index_out(start_in_out, ref_index_out); 
	const uint real_protein_out = listTaxa[taxon_out].getWorldIndex(protein_out);
	overlap_t overlap_right_left = 0; // The over vector from the right-to-left direction
	overlap_t recip_right_left_right = 0; // The over vector from the right-to-left direction
	const float sim_out_in = getOuterSimilarity(basket_out, taxon_out, taxon_in,
						    protein_out, protein_in,
						    overlap_right_left, recip_right_left_right);
	if(sim_out_in > 0) { // The relation exists: thereby do not ignore it
	  const overlap_t overlap_right_right = listTaxa[taxon_out].getOverlap(protein_out);
	  const overlap_t    overlap_left_left = listTaxa[taxon_in].getOverlap(protein_in);//arrOverlap[protein_in];
	  const overlap_t overlap_left_right = basket_inn->get_overlap_in(start_in_out, ref_index_out);
	  const overlap_t recip_left_right_left = basket_inn->get_overlap_out(start_in_out, ref_index_out);
	  bool is_above_overlap_limit = false;
	  if(use_improved_overlap_algo) is_above_overlap_limit = aboveOverlapLimit_improved(overlap_left_left, overlap_left_right, overlap_right_right, overlap_right_left, recip_left_right_left, recip_right_left_right, AMINO_LIMIT, PRINT_OVERLAP_VALUES_ABOVE);
	  else is_above_overlap_limit = aboveOverlapLimit(overlap_left_left, overlap_left_right, overlap_right_right, overlap_right_left, AMINO_LIMIT, PRINT_OVERLAP_VALUES_ABOVE);

	  if(is_above_overlap_limit)
	    { // Above the overlap limit set by the user
	      const float sim_in_out = basket_inn->getSimScore(start_in_out, ref_index_out, max_input_value); 
	      const float sim_avg = insert_relation(my_id, taxon_in, taxon_out,
						    protein_in, protein_out,
						    sim_in_out, sim_out_in);

	      if(!INPARALOG_OPERATION) { // Inserts a possible orhtolog pair:
		if(debug) {printf("(possible-ortholog-above-overlap-limit)\tpair (%s, %s), where sim_avg=%.3f must be above zero for being of interestat %s:%d\n",
				  listTaxa[taxon_in].getArrKey(protein_in), listTaxa[taxon_out].getArrKey(protein_out), sim_avg, __FILE__, __LINE__);}
		if(sim_avg > 0) {
		  if(sim_avg > this_sim_max) { 
		    while(!max_sim_key.empty()) {
		      const uint index_discarded = max_sim_key.top();
		      max_sim_key.pop();
		      if(DEBUG_PRINT_DISCARDED_PAIRS) {
			printf("(possible-ortholog-pair-discarded)\tpair (%s, %s) with averaged_similarity(%f) < the_better_limit(%f) defined for (%s, %s)\n",
			       listTaxa[taxon_in].getArrKey(protein_in), listTaxa[taxon_out].getProteinNameWorld(index_discarded), this_sim_max, sim_avg, listTaxa[taxon_in].getArrKey(protein_in), listTaxa[taxon_out].getProteinNameWorld(real_protein_out));
			
		      }
		    }
		    this_sim_max = sim_avg, max_sim_key.push(real_protein_out);
		  } else if(sim_avg == this_sim_max) {
		    max_sim_key.push(real_protein_out);
		  }
		}
	      }
	      if(DEBUG)  printf("\n");
	    } else {
	    if(DEBUG_PRINT_DISCARDED_PAIRS) {
	      assert(listTaxa);

	      if(false)	      printf("(discarded)\tpair (%u, %u) with overlap_values( [%d, %d, %d] and [%d, %d, %d] ) < overlap_limit(%d)\n", 					 real_protein_in, real_protein_out, overlap_left_left, overlap_left_right, recip_left_right_left, overlap_right_right, overlap_right_left, recip_right_left_right, AMINO_LIMIT);
	      else	      printf("(discarded)\tpair (%s, %s) with overlap_values( [%d, %d, %d] and [%d, %d, %d] ) < overlap_limit(%d)\n", 					 listTaxa[taxon_in].getArrKey(protein_in), listTaxa[taxon_out].getArrKey(protein_out), overlap_left_left, overlap_left_right, recip_left_right_left, overlap_right_right, overlap_right_left, recip_right_left_right, AMINO_LIMIT);
	    }
	  }
	}
      }
    }
    if(!INPARALOG_OPERATION) { // Inserts a possible orhtolog pair:
      listTaxa[taxon_in].insertInpaLimit(protein_in, this_sim_max); 	// Updates the limit list for the inparalogs:
      const uint size = max_sim_key.size();
      if(size > 0) {
	while(!max_sim_key.empty()) {
	  const uint this_protein_max = max_sim_key.top();
	  max_sim_key.pop(); 
	  if(debug) {printf("(possible-ortholog)\tpair (%s, %s), with sim_max=%.3f, at %s:%d\n",
			    listTaxa[taxon_in].getProteinNameWorld(real_protein_in), listTaxa[taxon_out].getProteinNameWorld(this_protein_max), this_sim_max, __FILE__, __LINE__);}
	  listOrtho.insertGlobalElement(real_protein_in, this_protein_max, this_sim_max);
	}
      }
    }
  }
}

/**! Produses a temporary list of arrOrtho elements
   @Description: Used in order to enhance the locality, and thereby do a speedup of the writing to memory
*/
/*
uint **createTempList(uint protein_start, uint protein_end) {
  const uint length_block = protein_end - protein_start;
  uint **arrOrtho = (uint**)malloc(sizeof(uint*)*length_block);
  for(uint i = 0; i< length_block; i++) arrOrtho[i] = NULL;
  return arrOrtho;
}
*/
/**! Updates 'listTaxa' with the arrOrtho, and finlizes the data
   @ProgramPoint: At the end of the operation in 'void operator'
*/
/*
void finalizeTempList(uint **arrOrtho, uint taxon_in, uint protein_start, uint protein_end, taxa_t *listTaxa) {
  int pos = 0;
  for(uint protein_in = protein_start; protein_in < protein_end; protein_in++) {
    listTaxa[taxon_in].arrOrtho[protein_in] = arrOrtho[pos], pos++;
  } free(arrOrtho);
}
*/
/*override*/void * pipe_ortholog_inparalog::operator()(void *item) {
  int my_id = -1;
  const clock_t tstart = times(NULL);
  taxon_pair_t *bucket = static_cast<taxon_pair*>(item);
  if(bucket) {
    {
      slock_t::scoped_lock updlock(lock_set_ps_id);
      for(uint i = 0; i< n_threads; i++) {
	if(in_use[i]==false) {in_use[i]= true, my_id=i,i =n_threads;}
      }
    }
    if(my_id > -1) {
      loint this_debug_cnt_possible_orthologs = 0;
      loint this_debug_cnt_possible_orthologs_inserted = 0;
      loint this_debug_cnt_possible_inparalogs = 0;
      loint this_debug_cnt_possible_inparalogs_inserted = 0;
      loint this_debug_proteins_evaluated = 0;
      uint taxon_in = 0, protein_end=0, protein_start = 0; 
      uint bucket_cnt = 0;
      uint debug_total_cnt[2] = {0,0};
      float debug_total_distance[2] = {0,0};
      uint debug_total_zeros[2] = {0,0};
      if(!USE_EVERYREL_AS_ARRNORM_BASIS) {
	l_arrNorm[my_id] = new list_norm(taxon_length, max_input_value, PRINT_NORMALIXATION_BASIS, DEBUG_NORM);
      }
      l_fileStruct[my_id] = list_file_parse<rel>::init_class(taxon_length, listTaxa, MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, FILE_BINARY_LOCATION);

      do {
	bucket[bucket_cnt++].getVariables(taxon_in, protein_start, protein_end); // The dataset to work ok    
	if(protein_end != 0) {
	  // One item is added to the list when there is one protein left for reading, ie, updates to the correct settings:
	  //	  if((uint)listTaxa[taxon_in].total_cnt == (protein_end+1)) protein_end = protein_end+1; 
	  //	  uint pairs_added = 0;
	  if(!INPARALOG_OPERATION) {
	    for(uint taxon_out = 0; taxon_out < taxon_length; taxon_out++) {
	      if(taxon_in != taxon_out) {
		basket_parse *basket_inn = listParseData->getBufferFromFile(taxon_in, taxon_out);
		if(basket_inn && basket_inn->has_data()) {
#ifndef NDEBUG
		  basket_inn->get_properties_and_increment(debug_total_cnt[0], debug_total_distance[0], debug_total_zeros[0]); 
#endif
		  this_debug_cnt_possible_orthologs += basket_inn->get_buffer_size();
		  basket_parse *basket_out = listParseData->getBufferFromFile(taxon_out, taxon_in);
		  if(basket_out && basket_out->has_data()) {
#ifndef NDEBUG
		    basket_out->get_properties_and_increment(debug_total_cnt[1], debug_total_distance[1], debug_total_zeros[1]); 
#endif
		    this_debug_proteins_evaluated += (protein_end - protein_start);
		    for(uint protein_in = protein_start; protein_in < protein_end; protein_in++) {
		      getProteinRelations(my_id, basket_inn, basket_out, taxon_out, taxon_in, listTaxa[taxon_in].getWorldIndex(protein_in));
		    }
#ifndef NDEBUG
		    this_debug_cnt_possible_orthologs_inserted += l_fileStruct[my_id]->getSizeofBuffer(taxon_in, taxon_out, protein_start, protein_end);
		    //		    pairs_added = (uint)l_fileStruct[my_id]->getSizeofBuffer(taxon_in, taxon_out, protein_start, protein_end);;
#endif
		  }
		  basket_parse::free_mem(basket_out);
		  basket_parse::free_mem(basket_inn);
		}
	      }
	    }
	  } else {
	    // The Inparalog operation: Iterating over the inner protein by using their 'world indexes':
	    basket_parse *basket_inn = listParseData->getBufferFromFile(taxon_in, taxon_in);
	    if(basket_inn && basket_inn->has_data()) {
	      this_debug_cnt_possible_inparalogs += basket_inn->get_buffer_size();
	      for(uint protein_in = protein_start; protein_in < protein_end; protein_in++) {
		getProteinRelations(my_id, basket_inn, basket_inn, taxon_in, taxon_in, listTaxa[taxon_in].getWorldIndex(protein_in));
	      }
#ifndef NDEBUG
	      this_debug_cnt_possible_inparalogs_inserted += l_fileStruct[my_id]->getSizeofBuffer(taxon_in, taxon_in, protein_start, protein_end);
	      //	      pairs_added = (uint)l_fileStruct[my_id]->getSizeofBuffer(taxon_in, taxon_in, protein_start, protein_end);
#endif
	      basket_parse::free_mem(basket_inn);
	    }
	  }
	}
      } while(protein_end != 0);
      
      
      taxon_pair::close(bucket);
      {
	// At the end:
	bucket_pipe_binary* bucket = (bucket_pipe_binary*)malloc(sizeof(bucket_pipe_binary));
	*bucket = bucket_pipe_binary(l_arrNorm[my_id], l_fileStruct[my_id]);
	l_arrNorm[my_id] = NULL; // As a safeguard.

	{
	  slock_t::scoped_lock updlock(lock_set_ps_id);
	  in_use[my_id] = false;
#ifndef NDEBUG
	  if(this_debug_cnt_possible_orthologs) {
	    debug_cnt_possible_orthologs += this_debug_cnt_possible_orthologs;
	    debug_cnt_possible_orthologs_inserted += this_debug_cnt_possible_orthologs_inserted;
	    assert(this_debug_cnt_possible_orthologs_inserted <= this_debug_cnt_possible_orthologs);
	  } else {
	    debug_cnt_possible_inparalogs += this_debug_cnt_possible_inparalogs;
	    debug_cnt_possible_inparalogs_inserted += this_debug_cnt_possible_inparalogs_inserted;
	    assert(this_debug_cnt_possible_inparalogs_inserted <= this_debug_cnt_possible_inparalogs);
	  }
#endif
	}
	const clock_t tend = times(NULL);
	logid_t id = filter_orthologs;
	if(INPARALOG_OPERATION) id = filter_inparalogs;
	log->append_measurement(id, 1, tstart, tend);
	//printf("(return-pipe)\t bucket-cnt=%u, at pipe_ortholog_inparalog:%d\n", bucket->structData->getTotalLengthOfData(), __LINE__);
	return bucket; 	
      }
    }
  }
  //  return NULL;
  bucket_pipe_binary* _bucket = (bucket_pipe_binary*)malloc(sizeof(bucket_pipe_binary));
  *_bucket = bucket_pipe_binary();
  return _bucket;
}


//! Frees the data bounded by this class
void pipe_ortholog_inparalog::free_data() { 
#ifndef NDEBUG
  if(!INPARALOG_OPERATION) {
    FILE *f_log = log_builder::get_log_file_pointer(n_threads, "pipe_process_build", __FILE__, __LINE__, true);
    if(f_log) {      
      fprintf(f_log, "-- \t Given %lld possible orthologs as input, ", debug_cnt_possible_orthologs);
      fprintf(f_log, "of these inserted %lld possible orthologs into new structure, to be used when building co-orthologs, giving a compression-factor of %f\n", debug_cnt_possible_orthologs_inserted, (float)((float)debug_cnt_possible_orthologs_inserted/debug_cnt_possible_orthologs));
      fclose(f_log);
    }
  } else {
    FILE *f_log = log_builder::get_log_file_pointer(n_threads, "pipe_process_build", __FILE__, __LINE__, false);
    if(f_log) {      
      fprintf(f_log, "-- \t Given %lld possible inparalogs as input, ", debug_cnt_possible_inparalogs);
      fprintf(f_log, "of these %lld classified inparalogs, giving a compression-factor of %f\n", debug_cnt_possible_inparalogs_inserted, (float)((float)debug_cnt_possible_inparalogs_inserted/debug_cnt_possible_inparalogs));
      fclose(f_log);
    }
  }  
#endif
  if(l_fileStruct) {
    for(uint i = 0; i < (uint)n_threads; i++) {
      if(l_fileStruct[i] != NULL) {
	list_file_parse<rel>::close(l_fileStruct[i], false);
      }
      assert(l_fileStruct[i] == NULL);
    }
    delete [] l_fileStruct; l_fileStruct = NULL;
  }
  list_norm::close(l_arrNorm, n_threads);
  free(in_use); in_use = NULL; 
}

//! De-allocates both the internal temporary objects of type list_file_struct for all of the threads, in addition to de-allocating the memory reserved for this object.
void pipe_ortholog_inparalog::free_additional_blocks() {
  list_file_parse<rel>::close(l_fileStruct, n_threads, false);
  free_data();
}

pipe_ortholog_inparalog::pipe_ortholog_inparalog(uint _taxon_length,
						 const bool _inparalog_operation, // set to true if its the inparalogs who shall be treated
						 const bool _use_everyrel_as_arrnorm_basis, const uint _n_threads, log_builder_t *_log,
						 id_simil_list &_listOrtho, list_file_parse_t *_listParseData, short int _AMINO_LIMIT,
						 float _max_input_value, float _MIN_SIMILARITY_LIMIT, bool _use_improved_overlap_algo,
						 bool _DEBUG_PRINT_DISCARDED_PAIRS, bool _PRINT_OVERLAP_VALUES_ABOVE,
						 bool _PRINT_NORMALIXATION_BASIS, bool _DEBUG_NORM, taxa_t *_listTaxa,
						 bool _MODE_PAIRWISE_OUTPUT_ABC, bool MODE_PAIRWISE_OUTPUT_MCL, char *_FILE_BINARY_LOCATION
						 )
  :
  tbb::filter(false), // true=serial, false=paralell
  PRINT_NORMALIXATION_BASIS(_PRINT_NORMALIXATION_BASIS), DEBUG_NORM(_DEBUG_NORM),
  debug_cnt_possible_orthologs(0), debug_cnt_possible_orthologs_inserted(0),
  debug_cnt_possible_inparalogs(0), debug_cnt_possible_inparalogs_inserted(0),
  listTaxa(_listTaxa),
  MODE_PAIRWISE_OUTPUT_ABC(_MODE_PAIRWISE_OUTPUT_ABC), MODE_PAIRWISE_OUTPUT_MCL(MODE_PAIRWISE_OUTPUT_MCL),
  FILE_BINARY_LOCATION(_FILE_BINARY_LOCATION),  PRINT_OVERLAP_VALUES_ABOVE(_PRINT_OVERLAP_VALUES_ABOVE), 
  use_improved_overlap_algo(_use_improved_overlap_algo), DEBUG_PRINT_DISCARDED_PAIRS(_DEBUG_PRINT_DISCARDED_PAIRS), 
  listOrtho(_listOrtho), listParseData(_listParseData), AMINO_LIMIT(_AMINO_LIMIT),
  max_input_value(_max_input_value), MIN_SIMILARITY_LIMIT(_MIN_SIMILARITY_LIMIT), 
  log(_log),
  taxon_length(_taxon_length), INPARALOG_OPERATION(_inparalog_operation),
  USE_EVERYREL_AS_ARRNORM_BASIS(_use_everyrel_as_arrnorm_basis),
  n_threads(_n_threads),
  list_pair_pos(0)
{

  l_fileStruct = new list_file_parse<rel>*[n_threads]; //*)malloc(sizeof(list_file_parse*)*size); // Holds the refernces to the binary files created under the parsing
  memset(l_fileStruct, 0, sizeof(list_file_parse<rel>*)*n_threads);
  //  l_fileStruct = list_file_parse<rel>::allocate_single_list(n_threads, 0, taxon_length, listTaxa, MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, FILE_BINARY_LOCATION);
  l_arrNorm = list_norm::init_list(n_threads); 
  in_use = (bool*)malloc(sizeof(bool)*n_threads);
  for(uint i = 0; i<n_threads; i++) {in_use[i]=0;}
}



#ifdef assert_code
void pipe_ortholog_inparalog::assert_private_parts() {
  // Tests overlap filtering.
  overlap_t overlap_left_left, overlap_right_right, overlap_protein_left, overlap_protein_right;
  short int AMINO_LIMIT = 1;
  bool PRINT_OVERLAP_VALUES_ABOVE = false;
  const uint old_amion_limit = AMINO_LIMIT;
  overlap_left_left = 2, overlap_right_right=2;
  AMINO_LIMIT = 0;
  overlap_protein_right = 2, overlap_protein_left = 2;
  assert(aboveOverlapLimit(overlap_left_left, overlap_right_right, overlap_protein_left, overlap_protein_right, AMINO_LIMIT, PRINT_OVERLAP_VALUES_ABOVE));
  AMINO_LIMIT = 300;
  overlap_protein_right = 2, overlap_protein_left = 2;
  assert(!aboveOverlapLimit(overlap_left_left, overlap_right_right, overlap_protein_left, overlap_protein_right, AMINO_LIMIT, PRINT_OVERLAP_VALUES_ABOVE));
  AMINO_LIMIT = old_amion_limit;
}
#endif
void pipe_ortholog_inparalog::assert_class(const bool print_info) {
  const static char *class_name = "pipe_ortholog_inparalog";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  // Iinitial tests:
  uint taxon_length = 2, n_threads = 1;
  bool inparalog_operation = true, use_everyrel = false;
  log_builder_t *log = new log_builder[1];//(NULL, '_', 0, 0, 0);
  id_simil_list listOrtho;  
  list_file_parse_t *listParseData = NULL; 
  short int AMINO_LIMIT = 2;
  float max_input_value = 10; // hold the highets input value read
  float MIN_SIMILARITY_LIMIT = 0;
  bool use_improved_overlap_algo = true;
  bool DEBUG_PRINT_DISCARDED_PAIRS = true;
  bool PRINT_OVERLAP_VALUES_ABOVE = false;
  bool PRINT_NORMALIXATION_BASIS = false; bool DEBUG_NORM = false; taxa_t *listTaxa = NULL;
  bool MODE_PAIRWISE_OUTPUT_ABC = true; bool MODE_PAIRWISE_OUTPUT_MCL = true; char *FILE_BINARY_LOCATION = NULL;
  pipe_ortholog_inparalog temp(taxon_length, inparalog_operation, use_everyrel, n_threads, log,
			       listOrtho, listParseData, AMINO_LIMIT, max_input_value, MIN_SIMILARITY_LIMIT,
			       use_improved_overlap_algo, DEBUG_PRINT_DISCARDED_PAIRS, PRINT_OVERLAP_VALUES_ABOVE,
			       PRINT_NORMALIXATION_BASIS, DEBUG_NORM, listTaxa, MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, FILE_BINARY_LOCATION
			       );
  temp.assert_private_parts();
  temp.free_additional_blocks();
  if(log) {log->free_mem(), delete log, log = NULL;}
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

