#include "pipe_struct.h"  
 
//! Clears the memory allocated for this thread
void pipe_struct::finalize_memory(const uint taxon_length) {
  if(PIPE_TYPE == DUMP)  {
#ifndef NDEBUG
    FILE *f_log = log_builder::get_log_file_pointer(n_threads, "pipe_dump", __FILE__, __LINE__, true);
    if(f_log) {      
      fprintf(f_log, "--\t The dumping operation resulted in the strings having the following number of chars:\n");
      debug_dump_result.print_result(f_log);
      fclose(f_log);
    }
    assert(lst_elements_evaluated);
    //! Sums the result for this node
    meta_pipe_struct_t result_this = meta_pipe_struct();
    for(uint i = 0; i < (uint)n_threads; i++) {
      result_this.cnt_names   += lst_elements_evaluated[i].cnt_names;
      result_this.cnt_inpa    += lst_elements_evaluated[i].cnt_inpa;
      result_this.cnt_ortho   += lst_elements_evaluated[i].cnt_ortho;
      result_this.cnt_co_orth += lst_elements_evaluated[i].cnt_co_orth;
    }
    //! Sums the number of relations in total
    meta_pipe_struct_t result_complete = meta_pipe_struct();
    if(number_of_nodes > 1) { // Then we are in 'mpi-mode':
#ifdef USE_MPI
      //! Sum the values of the inserted data:
      MPI_Allreduce(&result_this.cnt_names,   &result_complete.cnt_names, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&result_this.cnt_inpa,    &result_complete.cnt_inpa, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&result_this.cnt_ortho,   &result_complete.cnt_ortho, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&result_this.cnt_co_orth, &result_complete.cnt_co_orth, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
      assert(false); // Should not have come to this location if mpi is not used.
#endif
    } else { // We have all the interesting data in this node, ie, communication is not needed.
      result_complete.cnt_names   = result_this.cnt_names;
      result_complete.cnt_inpa    = result_this.cnt_inpa;
      result_complete.cnt_ortho   = result_this.cnt_ortho;
      result_complete.cnt_co_orth = result_this.cnt_co_orth;
    }

    //! Tests for the internal sums:
    if((result_this.cnt_inpa > 0) || (result_this.cnt_ortho > 0) || (result_this.cnt_co_orth > 0) ) {
      if(result_this.cnt_names != elements_expectation.cnt_names) {
	fprintf(stderr, "!!\t There might be an error in the orthAgogue-filtering:\n"
		"-\t The operation resulted in %u names, though there were %u names;\n"	      
		"-\t With in total %u inparalogs, %u orthologs, and %u co-orthologs;\n"
		"-\t The message was generated at at [%s]:%s:%d\n"
		"If this message is not understood, please forward it to the developer at oekseth@gmail.com",
		result_this.cnt_names, elements_expectation.cnt_names, 
		result_this.cnt_inpa, result_this.cnt_ortho, result_this.cnt_co_orth,
		__FUNCTION__, __FILE__, __LINE__);
      }
      assert(result_this.cnt_names == elements_expectation.cnt_names);
    } // else no names are wirtten, as no result were created.
    assert(result_this.cnt_inpa  == elements_expectation.cnt_inpa);
    assert(result_this.cnt_ortho == elements_expectation.cnt_ortho);
    assert(result_this.cnt_co_orth <= elements_expectation.cnt_co_orth);

    //! Tests for the total sums:
    if((result_this.cnt_inpa > 0) || (result_this.cnt_ortho > 0) || (result_this.cnt_co_orth > 0) ) {
      assert(result_complete.cnt_names == (uint)listTaxa[taxon_length-1].rel_end);
    }
    assert(result_complete.cnt_ortho == listOrtho.get_total_number_of_pairs());

    //! De-allocates the container:
    if(lst_elements_evaluated) {delete [] lst_elements_evaluated; lst_elements_evaluated = NULL;}
#endif
  } else {
    if((PIPE_TYPE == INPA_INPA) || (PIPE_TYPE == INPA_ORTH)) {
      if(!USE_EVERYREL_AS_ARRNORM_BASIS) list_norm::close(l_arrNorm, n_threads);
      assert(!l_arrNorm);
      //    free(l_arrNorm);
    }
  }
  free(listPair);
  delete [] mclData;
  delete [] in_use;
}


////! Builds the array for the normalization values doing the averaging
/*void pipe_struct::build_normArr(const uint taxon_length, norm_t **arrNorm) {  
  //  arrAvgNorm = norm_t::build_basis(taxon_length, arrNorm, PRINT_NORMALIXATION_BASIS, DEBUG_NORM, max_input_value);
  arrAvgNorm = norm_t::build_basis(taxon_length, arrNorm, PRINT_NORMALIXATION_BASIS, max_input_value);
}*/

//! Produces a row for the protein given
void pipe_struct::produce_row(uint my_id, uint protein_in, uint taxon_in) {
  assert(mclData);
  assert(my_id < n_threads);
  const uint world_index_in = listTaxa[taxon_in].getWorldIndex(protein_in);
  const uint size_outer = (uint)stackRel[world_index_in].unsafe_size(); 
  char *protein_name_in = listTaxa[taxon_in].getCompleteProteinName(protein_in);

  mclData[my_id]->set_name_index(protein_name_in, world_index_in); // sets the name of the protein (a seperate file)

  const uint debug_cnt_inparalogs_inserted = listStructData->produceInparalogs(mclData[my_id], protein_name_in, taxon_in, world_index_in, protein_in, arrAvgNorm[taxon_in][taxon_in]);
  uint debug_cnt_co_orthologs_inserted = 0;
  if(size_outer > 0) {
    mclData[my_id]->set_header_ortho_inpa(protein_name_in, world_index_in);
    o_rel_t out_pair;
    while(stackRel[world_index_in].try_pop(out_pair)) {
      const int world_out = out_pair.ind_out;
      const uint taxon_out = listTaxa[0].getTaxonId(world_out, listTaxa, taxon_length);
      const uint protein_out = listTaxa[taxon_out].getLocalIndex(world_out);
      char *protein_name_out = listTaxa[taxon_out].getCompleteProteinName(protein_out);
      const uint inserted_cnt = mclData[my_id]->insert_ortho_inpa(world_index_in, world_out, protein_name_in, protein_name_out, out_pair.distance, arrAvgNorm[taxon_in][taxon_out]);

#ifndef NDEBUG
      debug_cnt_co_orthologs_inserted += inserted_cnt;
      if(mclData[my_id]->has_sub_string(orth_inpa, protein_name_in) == false) {
	printf("(failed)\t\t inserts %s->%s not inserted, at pipe_struct:%d\n", protein_name_in, protein_name_out, __LINE__); 
	mclData[my_id]->print(orth_inpa); 
	assert(false);
      }
      if(mclData[my_id]->has_sub_string(orth_inpa, protein_name_out) == false) {
	printf("(failed)\t\t inserts %s->%s not inserted, at pipe_struct:%d\n", protein_name_in, protein_name_out, __LINE__); 
	mclData[my_id]->print(orth_inpa); 
	assert(false);
      }
#endif
      delete [] protein_name_out;
    }
  }
  //! The orthologs, named 'orthologs.*'
  const uint world_in = listTaxa[taxon_in].getWorldIndex(protein_in);
  uint size_o = 0; // size to be set by function call below:
  rel_t *ortho_p =  listOrtho.getGlobalRow(world_in, size_o);
  if(ortho_p != NULL) { // The protein has orthologs
    if(size_o > 0) {
#ifndef NDEBUG
      //! Get the start position of the string contaning the orthologs
      //      char *orthologs_start_position = mclData[my_id]->get_ortholog_current_eof(/*type_integer=*/true);
      //      uint estimated_chars_added = mclData[my_id]->get_size_of_header(protein_name_in, world_index_in);
#endif
      mclData[my_id]->set_header_pair_ortho(protein_name_in, world_index_in);
      // Add the relations:
      for(uint o = 0; o < size_o; o++) {
	const uint taxon_out   =  taxa::getTaxonId(ortho_p[o].ind_out, listTaxa, taxon_length);
	const int world_out = ortho_p[o].ind_out; const float sim_score = ortho_p[o].distance;
	assert(taxon_out < (uint)taxon_length);
	const uint protein_out = listTaxa[taxon_out].getLocalIndex(world_out);
	char *protein_name_out = listTaxa[taxon_out].getCompleteProteinName(protein_out);
	assert(protein_name_out); // todo: will this walsywa hold?
	mclData[my_id]->insert_pair_ortho(world_index_in, world_out, protein_name_in, protein_name_out, sim_score, arrAvgNorm[taxon_in][taxon_out]);
	delete [] protein_name_out;
      }
    }
  }
#ifndef NDEBUG
  //! Procedure counting elements, as seend handled sperately for each thread id:
  assert(lst_elements_evaluated);
  //! Increments the number of 'root values' found:  
  lst_elements_evaluated[my_id].cnt_names++;
  //! Inserts the number of inparalogs (later to be checked towards our expectations):
  lst_elements_evaluated[my_id].cnt_inpa += debug_cnt_inparalogs_inserted;  
  //! Inserts the number of co-orthologs (later to be checked towards our expectations):
  lst_elements_evaluated[my_id].cnt_co_orth += debug_cnt_co_orthologs_inserted;
  //! Inserts the number of orthologs (later to be checked towards our expectations):
  lst_elements_evaluated[my_id].cnt_ortho += size_o;
#endif

  delete [] protein_name_in;
  mclData[my_id]->set_line_end();
}

/**! Inserts the similarity score into the stack
   @Changed: 18.04.2011 by oekseth, added possible debug outprint when inserted
*/
void pipe_struct::insert_into_stack_rel(uint my_id, uint taxon_in, uint taxon_out, uint world_in, uint world_out, float sim_score) {
  assert(!(listOrtho.hasRelation(world_in, world_out))); // Validates that an ortholog is not inserted.
  stackRel[world_in].push(o_rel(world_out, sim_score));
  stackRel[world_out].push(o_rel(world_in, sim_score));  
  if(!USE_EVERYREL_AS_ARRNORM_BASIS) { // Not already set, therefor sets it
    const uint protein_in = listTaxa[taxon_in].getLocalIndex(world_in);
    const uint protein_out = listTaxa[taxon_out].getLocalIndex(world_out);
    if(l_arrNorm[my_id]) {
      l_arrNorm[my_id]->insert(taxon_in,  taxon_out, sim_score, protein_in,  protein_out, listTaxa);
      l_arrNorm[my_id]->insert(taxon_out, taxon_in,  sim_score, protein_out, protein_in,  listTaxa);
    }
  }
}
/**! Insert the orthologs for the for the pair "protein_in <==> protein_out"
   @ProgramPoint: if 'PIPE_TYPE==INPA_ORTH' the called from the 'void operator'
   @Description: The stages as:
   1. Iterate through the inparalogs of 'protein_in'
   2. Test if 'protein_out' has any inparalogs bound to the above inparalog
   @Comment: The difference from the 'inparalog' operation, is that this function
   only get th the list of out-paralogs for the orholog_pair, i.e.
   ->getBuffer(taxon_out, taxon_in, protein_out).
*/
uint debug_cnt = 0; 
void pipe_struct::insertOrthologs(uint my_id, rel_t *buff_in_in, uint taxon_in, uint taxon_out, uint world_in, uint world_out, uint protein_in, uint protein_out) {
  rel_t *buff_out_in = listStructData->getBuffer(taxon_out, taxon_in, protein_out);
  const uint out_k_length = listStructData->getLength(taxon_out, taxon_in, protein_out);
  const uint inp_k_length = listStructData->getLength(taxon_in, taxon_in, protein_in);
  for(uint inp_k = 0; inp_k < inp_k_length; inp_k++) { // Iterating thorugh the inparalogs
    //    bool relation_found = false;
    for(uint out_k = 0; out_k < out_k_length; out_k++) { // Iterating thorugh the inparalogs
      if(buff_out_in[out_k].ind_out == buff_in_in[inp_k].ind_out) { // (the outer orthopairprotein)->inparalog-match or a match towrads the ortho-pair
	const uint world_in_inpa = listTaxa[taxon_in].getWorldIndex(buff_in_in[inp_k].ind_out);
	if(!(listOrtho.hasRelation(world_in_inpa, world_out))) { // Ensures that an ortholog is not inserted.
	  insert_into_stack_rel(my_id, taxon_in, taxon_out, world_in_inpa, world_out, buff_out_in[out_k].distance);
	}
	//	relation_found = true;
      } else if (buff_out_in[out_k].ind_out == protein_in) { // the ortholog pair, ie, do nothing as its already added.
	;
      } else if(DEBUG_PRINT_DISCARDED_PAIRS) {
	const uint world_in_inpa = listTaxa[taxon_in].getWorldIndex(buff_in_in[inp_k].ind_out);
 	printf("(co-orthologs-discarded-due-to-non-reciprocallity)\tpair (%s, %s)\n", 
 	       listTaxa[taxon_in].getProteinNameWorld(world_in_inpa), 
 	       listTaxa[taxon_out].getProteinNameWorld(world_out));
      }
    }
  }
  //  listStructData->free_getBuffer(buff_out_in, taxon_out, taxon_in, protein_out);
}

//! Usage: Used in the function below
struct relation {
  uint world_in_inpa;
  uint world_out_inpa;
  uint index_out;
  float distance;
  relation() :       world_in_inpa(0), world_out_inpa(0), index_out(0), distance(0.)     {};
  relation(uint w_in, uint w_out, uint ind_out, float dist) :
    world_in_inpa(w_in), world_out_inpa(w_out), index_out(ind_out), distance(dist)     {};
};

/**! Insert the inpa_inpa_relations
   @ProgramPoint: if 'PIPE_TYPE==INPA_INPA' the called from teh 'void operator'
   @Description: The stages as:
   #  Test if 'protein_out' has any inparalogs: 
   - If this holds to true, read the 'out_out' inparalogs from the buffer, and continue  
   ##  Iterate through the inparalogs of 'protein_in'
   ###  Test if 'protein_out' has any inparalogs bound to the above inparalog:
   - If this holds to true, continue by retrieving the inpralogs for 'protein_out' from the store
   ####  Iterate over 'the in_out relations' for each given inparalog:
   - If it is not equal to 'protein_out', iterate over 'protein_out's inparalogs
   - Else if it is equal to 'protein_out', add it to the list of relations: jum back to the previous step
   #######	 Iterate over the 'out_out' relations: if ([in_out].ind_out == [out_out].ind_out), add the relation to the list of relations
*/
void pipe_struct::insertInpaInpaRelations(uint my_id, rel_t *buff_in_in, uint taxon_in, uint taxon_out, /*int world_in,*/ uint world_out, uint protein_in, uint protein_out) {
  stack<relation> rel;
  assert(listStructData);
  assert(listTaxa);
  const uint inp_k_length = listStructData->getLength(taxon_in, taxon_in, protein_in);
  for(uint inp_k = 0; inp_k < inp_k_length; inp_k++) { // Iterating thorugh the inparalogs of 'protein_in'
    const uint index_inpa = buff_in_in[inp_k].ind_out;
    if(listTaxa[taxon_in].is_a_local_key(index_inpa)) {
      const uint world_in_inpa = listTaxa[taxon_in].getWorldIndex(index_inpa);
      if(listStructData->has_data(taxon_in, taxon_out, index_inpa)) { // If the inparalog has data to its sisters ortholog:
	//      bool relation_found = false;
	rel_t *buff_in_out = listStructData->getBuffer(taxon_in, taxon_out, index_inpa);
	const uint in_out_k_length = listStructData->getLength(taxon_in, taxon_out, index_inpa);
	for(uint in_out_k = 0; in_out_k < in_out_k_length; in_out_k++) { 	// Iterate over 'the in_out relations' for each given inparalog:
	  const uint index_out = buff_in_out[in_out_k].ind_out;
	  const uint world_out_inpa = listTaxa[taxon_out].getWorldIndex(index_out);
	  const float distance = buff_in_out[in_out_k].distance;
	  if(index_out != protein_out) { // It is not an inparalog-ortholog relation	  	  
	    if(taxon_out > taxon_in) { // at not at the other side of the ortholog pair
	      if(!(listOrtho.hasRelation(world_in_inpa, world_out_inpa))) { // Preventing a dijoint data-set.
		rel.push(relation(world_in_inpa, world_out_inpa, index_out, distance)); // a possible relation, but still nead a confirmation that the 'other side' of the ortholog pair has this relation as wel
	      }
	    }
	  
	  } else { // '(inparalog)=>ortholog' relation
	    if(!(listOrtho.hasRelation(world_in_inpa, world_out))) { // Preventing a dijoint data-set.
	      insert_into_stack_rel(my_id, taxon_in, taxon_out, world_in_inpa, world_out, buff_in_out[in_out_k].distance);
	      //	    relation_found = true;
	    }
	  }
	}
	//	listStructData->free_getBuffer(buff_in_out, taxon_in, taxon_out, protein_in);
      }
    }
#ifndef NDEBUG // Only prints the information if in DEBUG mode
    else {
      fprintf(stderr, "!!\tAn unexpected beaviour, as an inparalogs index is used, bot not defined in the list of its names. This error was generated at line %d in file %s, located in method %s. Please contact the developer at oekseth@gmail.com\n", __LINE__, __FILE__, __FUNCTION__); 
    }
#endif
  }    
  if(!rel.empty()) { // has inpa_inpa_relations
    const uint out_out_k_length = listStructData->getLength(taxon_out, taxon_out, protein_out);
    rel_t *buff_out_out = NULL;
    if(out_out_k_length > 0) { //  out_out[protein_out] has inparalogs:
      // If this holds to true, read the 'out_out' inparalogs from the buffer, and continue  
      buff_out_out = listStructData->getBuffer(taxon_out, taxon_out, protein_out);
      if (buff_out_out != NULL) { // 
	while(!rel.empty()) {
	  const float distance = rel.top().distance;
	  const uint index_out = rel.top().index_out, world_in_inpa = rel.top().world_in_inpa, world_out_inpa = rel.top().world_out_inpa; rel.pop();
	  bool is_inserted = false; // updateted, for use with the verbouse output, it protein used.
	  for(uint out_out_k = 0; out_out_k < out_out_k_length; out_out_k++) { // Iterating thorugh the inparalogs
	    if(buff_out_out[out_out_k].ind_out == index_out) { // ([in_out].ind_out == [out_out].ind_out), add the relation to the list of relations
	      assert(!(listOrtho.hasRelation(world_in_inpa, world_out_inpa))); // { // Preventing a dijoint data-set.
	      insert_into_stack_rel(my_id, taxon_in, taxon_out, world_in_inpa, world_out_inpa, distance);
	      is_inserted=true;
	    }
	  }
	  if(!is_inserted && DEBUG_PRINT_DISCARDED_PAIRS) {
	    printf("(co-orthologs-restricted-discarded-due-to-non-reciprocallity)\tpair (%s, %s)\n", 
		   listTaxa[taxon_in].getProteinNameWorld(world_in_inpa), 
		   listTaxa[taxon_out].getProteinNameWorld(world_out_inpa));
	  }
	}
      }
      //      listStructData->free_getBuffer(buff_out_out, taxon_out, taxon_out, protein_out);
    }
  }
}

/*override*/void* pipe_struct::operator()(void* item) {
  if(false) printf("item(%ld)\n", (long int)item); // Hides compiler message about unused variable
  const clock_t tstart = times(NULL);
  const uint taxon_length = g_taxon_length;
  int my_id = -1;
  uint this_list_pair_pos = 0;
  taxon_pair *bucket = NULL;
  {
    slock_t::scoped_lock updlock(lock_set_ps_id);
    {
      uint i = 0;
      for(i = 0; i< n_threads; i++) {
	if(in_use[i]==false) {in_use[i]= true, my_id=i, i = n_threads;}
      }
    }
    bucket = taxon_pair::get_taxon_pairs(listPair, list_pair_pos);      
    this_list_pair_pos = (uint)list_pair_pos;
  } 
  if(my_id == -1) { // Unlikly errors of importance should be tested- and caught as a safeguard:
    fprintf(stderr, "!!\t\"my_id = -1 \": An error occured at line %d in file %s, found in method %s. Please contact the developer at oekseth@gmail.com if this error is seen.\n",
	    __LINE__, __FILE__, __FUNCTION__);
    assert(my_id != -1);
    return NULL;
  }
  //  if (listPair && protein_end > 0 ) 
  if(bucket) {
    if(PIPE_TYPE == DUMP) {
      mclData[my_id] = mcl_format::init(taxon_length, listTaxa, DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT, DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT,
			  MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT,
			  SORT_ABC_DATA);
      uint debug_length_of_data = 0;
      uint taxon_in = 0, protein_end=0, protein_start = 0; 
      uint bucket_cnt = 0;
      do {
	bucket[bucket_cnt++].getVariables(taxon_in, protein_start, protein_end); // The dataset to work ok    
	if(protein_end != 0) {
	  for(uint protein_in = protein_start; protein_in < protein_end; protein_in++) {
	    produce_row(my_id, protein_in, taxon_in);
	  }
	  debug_length_of_data += listStructData->getProteinLength(taxon_in, true, protein_start, protein_end+1);
	}
      } while(protein_end != 0); 
      taxon_pair::close(bucket);
      {slock_t::scoped_lock updlock(lock_set_ps_id);
	in_use[my_id] = false;
	const clock_t tend = times(NULL);
#ifndef NDEBUG
	mclData[my_id]->log_update_result_cheme(debug_dump_result);
	debug_dump_result.append_meta_data(this_list_pair_pos, protein_end+1-protein_start);
	debug_dump_result.append_inparalog_counter(debug_length_of_data);
#endif
	log->append_measurement(write_resultfile, 1, tstart, tend);
      } // frees id for next run.
      return mclData[my_id]; 
    } else if(PIPE_TYPE == INPA_INPA || PIPE_TYPE == INPA_ORTH) {
      bool inpa_inpa = true;  if(PIPE_TYPE == INPA_ORTH) inpa_inpa = false;
      if(!USE_EVERYREL_AS_ARRNORM_BASIS) {
	assert(!l_arrNorm[my_id]);
	l_arrNorm[my_id] = new list_norm(taxon_length, max_input_value, PRINT_NORMALIXATION_BASIS, DEBUG_NORM);  
	assert(l_arrNorm[my_id]);
      }
      uint taxon_in = 0, protein_end=0, protein_start = 0; 
      uint bucket_cnt = 0;
      do {
	bucket[bucket_cnt++].getVariables(taxon_in, protein_start, protein_end); // The dataset to work ok    
	if(protein_end != 0) {
	  for(uint protein_in = protein_start; protein_in < protein_end; protein_in++) {
	    if(listStructData->has_data(taxon_in, taxon_in, protein_in)) {
	      rel_t *buff_in_in = listStructData->getBuffer(taxon_in, taxon_in, protein_in);
	      if(buff_in_in != NULL) {
		uint size_o = 0;
		const uint world_in = listTaxa[taxon_in].getWorldIndex(protein_in);
		rel_t *ortho_p = listOrtho.getGlobalRow(world_in, size_o);
		if(ortho_p != NULL) { // The protein has orthologs
		  // Builds the set of orthologs
		  for(uint o = 0; o < size_o; o++) {
		    assert(listOrtho.hasRelation(world_in, ortho_p[o].ind_out));
		    const uint taxon_out   =  taxa::getTaxonId(ortho_p[o].ind_out, listTaxa, taxon_length);
		    const uint protein_out = listTaxa[taxon_out].getLocalIndex(ortho_p[o].ind_out);
		    const float sim_score = ortho_p[o].distance;;
		    if(!USE_EVERYREL_AS_ARRNORM_BASIS && l_arrNorm && l_arrNorm[my_id]) l_arrNorm[my_id]->insert(taxon_in, taxon_out, sim_score, protein_in, protein_out, listTaxa); // This is directional, as the other direction later will be added
		    
		    // No need to continue our prcessing with this relation, if the ortholog apir doesnt have any data. (Suprsiingly many times, this if sentences is not visited!)
				
		    if(listStructData->has_data(taxon_out, taxon_in, protein_out)) {
		      if(inpa_inpa) insertInpaInpaRelations(my_id, buff_in_in, taxon_in, taxon_out, ortho_p[o].ind_out, protein_in, protein_out);
		      else 		        insertOrthologs(my_id, buff_in_in, taxon_in, taxon_out, world_in, ortho_p[o].ind_out, protein_in, protein_out);
		    }
	    
		  }	    
		}
	      }
	    }
	  }
	}
      } while(protein_end != 0); 
      taxon_pair::close(bucket);
      {slock_t::scoped_lock updlock(lock_set_ps_id);
	in_use[my_id] = false;
	const clock_t tend = times(NULL);
	log->append_measurement(filter_co_orthologs, 1, tstart, tend);
      } // frees id for next run.
      if(!USE_EVERYREL_AS_ARRNORM_BASIS) {
	bucket_norm *bucket = bucket_norm::init(l_arrNorm[my_id]);
	l_arrNorm[my_id] = NULL;
	return bucket;
      } else {
	bucket_norm *bucket = bucket_norm::init(NULL);
	return bucket;
      }
    }
  }
  if(my_id !=-1){  {slock_t::scoped_lock updlock(lock_set_ps_id); in_use[my_id] = false;}} // frees id for next run.
  return NULL;
}

void pipe_struct::free_mem() {
  if(!USE_EVERYREL_AS_ARRNORM_BASIS) list_norm::close(l_arrNorm, n_threads);
  if(in_use != NULL) {
    delete [] in_use; in_use = NULL;
  }
}

 

pipe_struct::pipe_struct(uint _nthread, uint taxon_start, uint _taxon_length, pipe_t type, bool _USE_EVERYREL_AS_ARRNORM_BASIS, list_norm_t *normArr, log_builder_t *_log,
			 bool _MODE_PAIRWISE_OUTPUT, bool _MODE_INTEGER_OUTPUT, bool _DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT,
			 bool _DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT, id_simil_list &_listOrtho, taxa_t *_listTaxa, stack_rel *&_stackRel, list_file_struct_t *&_listStructData, list_file_parse_t *_listParseData, float _max_input_value, bool _MODE_PAIRWISE_OUTPUT_ABC, // for the abc file: if set, the data out pairwise in stead of as in a row
	      bool _MODE_PAIRWISE_OUTPUT_MCL, // for the mcl file: if set, the data out pairwise in stead of as in a row
	      // The following variables decides what output is to be printed
	      bool _PRINT_IN_ABC_FORMAT,
	      bool _PRINT_IN_MCL_FORMAT,
	      bool _SORT_ABC_DATA, // if true, sorts the abc files before outprint
	      bool _PRINT_NORMALIXATION_BASIS,
			 bool _DEBUG_NORM, float **&arr_avgNorm) :
  tbb::filter(/*serial=*/false),
#ifndef NDEBUG
  lst_elements_evaluated(NULL),
#ifdef USE_MPI
  myrank(0),
#endif
#endif
  number_of_nodes(0),
  //  tbb::filter(/*ordered=*/EXECUTE_IN_SERIAL),
  MODE_PAIRWISE_OUTPUT_ABC(_MODE_PAIRWISE_OUTPUT_ABC),  MODE_PAIRWISE_OUTPUT_MCL(_MODE_PAIRWISE_OUTPUT_MCL), 
  PRINT_IN_ABC_FORMAT(_PRINT_IN_ABC_FORMAT), PRINT_IN_MCL_FORMAT(_PRINT_IN_MCL_FORMAT),
  SORT_ABC_DATA(_SORT_ABC_DATA), PRINT_NORMALIXATION_BASIS(_PRINT_NORMALIXATION_BASIS),
  DEBUG_PRINT_DISCARDED_PAIRS(false), DEBUG_NORM(_DEBUG_NORM),  taxon_length(_taxon_length),
  MODE_PAIRWISE_OUTPUT(_MODE_PAIRWISE_OUTPUT), MODE_INTEGER_OUTPUT(_MODE_INTEGER_OUTPUT),
  DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT(_DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT),
  DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT(_DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT),
  listOrtho(_listOrtho), listTaxa(_listTaxa), stackRel(_stackRel),
  listStructData(_listStructData), listParseData(_listParseData), 
  max_input_value(_max_input_value), 
  log(_log),
  n_threads(_nthread),
  //  g_protein_length(_size_prot),
  g_taxon_length(taxon_length),
  SIZE_BUFFER(0),
  list_pair_pos(0),
  PIPE_TYPE(type),
  USE_EVERYREL_AS_ARRNORM_BASIS(_USE_EVERYREL_AS_ARRNORM_BASIS)
{
  if(PIPE_TYPE == DUMP)  {
    arrAvgNorm = arr_avgNorm;
    log_builder::test_memory_condition_and_if_not_abort(arrAvgNorm!=NULL, __LINE__, __FILE__, __FUNCTION__);
    l_arrNorm = NULL;
    mclData = new mcl_format_t*[n_threads];
    debug_dump_result = pipe_struct_result();
    assert(listStructData);
    debug_dump_result.set_cnt_total_elements_in_list_file_struct(listStructData->getTotalLengthOfData());
#ifndef NDEBUG
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes);
#endif
    //! Expectations: Initiates the storage for metadata
    //! Note: The numbers set are for this node:
    elements_expectation = meta_pipe_struct();
    assert(listTaxa);
    assert(taxon_length);
    assert(taxon_start < (uint)taxon_length);
    uint cnt_names_this = 0;
#ifdef USE_MPI
    bool *taxon_is_relevant = listStructData->get_list_of_nodes_taxa_responsilibties();
    assert(taxon_is_relevant);
    for(uint i = (uint)taxon_start; i < (uint)taxon_length; i++) {
      if(taxon_is_relevant[i]) cnt_names_this += listTaxa[i].total_cnt;
    }
#else
    for(uint i = (uint)taxon_start; i < (uint)taxon_length; i++) {
      cnt_names_this += listTaxa[i].total_cnt;
    }
#endif
    //    const uint cnt_names_this = listTaxa[taxon_length-1].rel_end - listTaxa[taxon_start].rel_start;
    elements_expectation.cnt_names = cnt_names_this;
    elements_expectation.cnt_inpa  = listStructData->getTotalLengthOfData_for_inparalogs_myrank(taxon_start, taxon_length);
#ifdef USE_MPI
    elements_expectation.cnt_ortho = listOrtho.get_total_number_of_pairs(taxon_start, taxon_length, listStructData->get_list_of_nodes_taxa_responsilibties());
#else
    elements_expectation.cnt_ortho = listOrtho.get_total_number_of_pairs(taxon_start, taxon_length, NULL);
#endif
    loint rel_cnt = 0; for(int i =0; i< listTaxa[taxon_length-1].rel_end;i++) rel_cnt+=stackRel[i].unsafe_size();
    elements_expectation.cnt_co_orth = rel_cnt;

    //! Actual evaluated: Initiates the storage for metadata
    assert(!lst_elements_evaluated);
    lst_elements_evaluated = new meta_pipe_struct[n_threads];
    assert(lst_elements_evaluated);
    for(uint i = 0; i < (uint)n_threads; i++) {
      lst_elements_evaluated[i]    = meta_pipe_struct();
    }
#endif

  } else {
    mclData = NULL;
    l_arrNorm = NULL;
    if(((PIPE_TYPE == INPA_INPA) || (PIPE_TYPE == INPA_ORTH))) {
      if(!USE_EVERYREL_AS_ARRNORM_BASIS) l_arrNorm = list_norm::init_list(n_threads); 
      arrAvgNorm = NULL;
    }
  }
  _LENGTH_KEY_NAMES = (int)log10(UINT_MAX)+1; // the maximum numbers of chars it can hold
  _WIDTH_DISTANCE = (int)log10(max_input_value)+1+4; // adds n addition the precision specified: '2' decimals
  
  in_use = new bool[n_threads];
  for(uint i = 0;i<n_threads;i++) in_use[i] = false;
  //  uint taxon_start = 0,
  uint taxon_end = taxon_length;
  if(listParseData != NULL) fprintf(stderr, "in pipe_struct only uses the listStructData, i.e. listParseData should be emptied (set to NULL, but is not..)\n");
  if(true) assert(listParseData==NULL);
  listPair = taxon_pair::init_taxon_pair(taxon_start, taxon_end, taxon_length, n_threads,
					 true, // only the inpralogs
					 NULL, listStructData, listTaxa
					 );
  
} 
#ifdef assert_code
void pipe_struct::assert_private_parts() {
  ;
}
#endif
void pipe_struct::assert_class(const bool print_info) {
  const static char *class_name = "pipe_struct";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  /*
  bool MODE_PAIRWISE_OUTPUT_ABC = true; // for the abc file: if set, the data out pairwise in stead of as in a row
  bool MODE_PAIRWISE_OUTPUT_MCL = true; // for the mcl file: if set, the data out pairwise in stead of as in a row
  // The following variables decides what output is to be printed
  bool PRINT_IN_ABC_FORMAT = true;
  bool PRINT_IN_MCL_FORMAT = true;
  bool SORT_ABC_DATA  = false; // if true, sorts teh abc files before outprint
  bool PRINT_NORMALIXATION_BASIS = false;
  bool DEBUG_NORM = false;
  log_builder_t *log = new log_builder[1];//(NULL, '_', 0, 0, 0);
  bool MODE_PAIRWISE_OUTPUT = true; // If set, prints the data out pairwise in stead of as in a row
  bool MODE_INTEGER_OUTPUT = true; // If set, prints the data out as integers in stead of as in a row
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT = true;
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT = true;
  id_simil_list listOrtho;  
  taxa_t *listTaxa = NULL;
  stack_rel *stackRel = NULL;// Holds the co orthologs for every protein in a stack; intialized in 'ortho_set' 
  // Data arrays to intialize and fill
  list_file_struct_t *listStructData = NULL; // Holds the data about the structs; First Initialises in 'ortho_set', and then filled with data during this'- and the inparalog stepsBBBB
  list_file_parse_t *listParseData = NULL;
  float max_input_value = 0.0;
  //  if(false)   print_constants();
  if(log) {log->free_mem(), delete log, log = NULL;}
*/
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

