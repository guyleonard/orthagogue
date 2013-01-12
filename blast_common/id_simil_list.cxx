#include "id_simil_list.h"
/**! A linked list class
   @Changed: 21.12.2010 by oekseth
*/


/**
   @Description:   Initiates a row of the '**arr'
   @ProgramPoint:  At is starting point
   not malloced before at this point (after the iteration of the orthologs,
   but before the iteration of the inparalogs).
*/
void id_simil_list::initGlobalRow(uint index_in, uint size) {
  arr[index_in] = (rel_t*)malloc(sizeof(rel_t)*(size));
  test_condition_and_if_not_abort(arr[index_in]!=NULL, __LINE__, __FUNCTION__);
  for(uint i = 0; i<size; i++) arr[index_in][i] = rel_t();
  arr[index_in][0] = rel_t(0, (float)(size-1)); // In order to know the size of this list
}

/**! Initalizes the data to zero
   @Comment: Does not de allocate the memory
*/
void id_simil_list::clearLocalData(rel_t *&arr) {
  if(arr != NULL) {
    const uint size = (uint)arr[0].distance;
    for(uint i = 1; i<size; i++) arr[i] = rel_t();
    arr[0].ind_out = 0; // th next pos thereby becomes '1'
  }
}

/**
   @Description:   Initiates a row of the '**arr'
   @ProgramPoint:  At is starting point
   not malloced before at this point (after the iteration of the orthologs,
   but before the iteration of the inparalogs).
*/
void id_simil_list::initLocalRow(uint index_in, uint size, rel_t **&arr) {
  if(arr) {
    arr[index_in] = (rel_t*)malloc(sizeof(rel_t)*(size));
    test_condition_and_if_not_abort(arr[index_in]!=NULL, __LINE__, __FUNCTION__);
    for(uint i = 0; i<size; i++) arr[index_in][i] = rel_t();
    arr[index_in][0] = rel_t(0, (float)(size-1)); // In order to know the size of this list
  } else   log_builder::test_memory_condition_and_if_not_abort(arr!=NULL, __LINE__, __FILE__, __FUNCTION__);
}

//! Returns the size of the local array (given as parameter) at the specified index
uint id_simil_list::getLocalSize(uint index_in, rel_t **&arr) {
  if(arr && arr[index_in] != NULL) 
    return (uint)arr[index_in][0].distance;
  else return 0;
}

//! Returns the size of the global array at the specified index
uint id_simil_list::getGlobalSize(uint index_in) {
  if(arr && arr[index_in] != NULL) 
    return (uint)arr[index_in][0].distance;
  else return 0;
}

/**! Allocates, and intiates the global array to the length specified */
void id_simil_list::initGlobalArr(uint _arr_size, uint _index_size) {
  assert(!arr);
  arr_size = _arr_size;
  arr = NULL;
  initLocalArr(arr, arr_size, _index_size); 
}

/**! Allocates, and intiates the array, given as input parameter, to the length specified */
void id_simil_list::initLocalArr(rel_t **&arr, uint arr_size, uint index_size) {
  assert(arr_size);
  arr = (rel_t**)malloc(sizeof(rel_t*)*arr_size);
  test_condition_and_if_not_abort(arr!=NULL, __LINE__, __FUNCTION__);
  for(uint i = 0; i< arr_size; i++) {
    if(index_size != 0) {
      initLocalRow(i, index_size, arr);
    } else arr[i] = NULL;
  }
}

//! Frees the memory allocated for the global arr
void id_simil_list::freeGlobalArr() {
  if(arr && arr_size) {
    for(uint i = 0; i< arr_size; i++) {
      free(arr[i]); arr[i] = NULL;
    }
    free(arr); arr = NULL;  arr_size = 0;
  }
}

//! Frees the memory allocated for the global arr
void id_simil_list::freeLocalArr(rel_t **&arr, uint arr_size) {
  if(arr != NULL) {
    for(uint i = 0; i< arr_size; i++) {
      free(arr[i]); arr[i] = NULL;
    }
    free(arr); arr = NULL; 
  }
}

//! Inserts an elements into the global array
void id_simil_list::insertGlobalElement(uint protein_in, uint world_protein_out, float sim_score) {
  insertLocalElement(arr, protein_in, world_protein_out, sim_score);
}
/**
   @Description:   Inserts an element into a row of the '**arrOrtho'
   @ProgramPoint:  Is first needed in the 'ortho_set_bin', therefore
   not malloced before at this point (after the iteration of the orthologs,
   but before the iteration of the inparalogs).
*/
void id_simil_list::insertLocalElement(rel_t **&arr, uint protein_in, uint world_protein_out, float sim_score) {
  assert(arr);
  if(arr[protein_in] == NULL) initLocalRow(protein_in, 10, arr);
  const uint next_index = arr[protein_in][0].ind_out+1; // In addition the offset
  if((uint)arr[protein_in][0].distance == next_index) { // at the limit; an enlargment is neccessary
    const uint new_size =   (uint)arr[protein_in][0].distance + BASE_ARR_STEP;
    arr[protein_in] = (rel_t*)realloc(arr[protein_in], (sizeof(rel_t)*(new_size))); // enlarges
    assert(arr[protein_in]);
    for(uint i = next_index; i<new_size; i++) arr[protein_in][i] = rel_t();//    memset(arrOrtho[protein_in]+next_index, 0, new_size);
    arr[protein_in][0].distance = (float)new_size; 
  }
  arr[protein_in][0].ind_out = next_index; // Updates the size (number of ortholog pairs)
  arr[protein_in][next_index] = rel_t(world_protein_out, sim_score); // Sets the data
}

/**
   @Description:  Returns a list (array) of the proteins 'protein_in'
   has ortholog pairs with
   @Param '&size' Sets the size of the array returned
*/
rel_t* id_simil_list::getGlobalRow(uint protein_in, uint &size) {
  return getLocalRow(arr, protein_in, size);
}


/**
   @Description:  Returns a list (array) of the proteins 'protein_in' from the array given as input parameter
   has ortholog pairs with
   @Param '&size' Sets the size of the array returned
*/
rel_t* id_simil_list::getLocalRow(rel_t **&arr, uint protein_in, uint &size) {
  if(arr[protein_in] != NULL) {
    size = arr[protein_in][0].ind_out; // The number of elements
    
    //printf("....%u has %u orthologs...\n", protein_in, size);
    return arr[protein_in]+1; // to avoid the extra field with the size;
    
  } else {
    // printf("....%u has zero orthologs...\n", protein_in);
    size = 0;
    return NULL;
  }
}

//! Returns true if the index given as parameter has data set for ti
bool id_simil_list::hasLocalData(rel_t **&arr, uint index_in) {
  if(arr != NULL) {
    if(arr[index_in] != NULL)
      return true;
  }
  return false;
}

//! Returns true if the index given as parameter has data set for ti
bool id_simil_list::hasGlobalData(uint index_in) {
  return hasLocalData(arr, index_in);
}

//! Returns true if 'index_out' has 'index_in' as a relation
bool id_simil_list::hasRelation(rel_t **&arr, uint index_in, uint index_out) {

  uint size = 0;
  rel_t *row = getLocalRow(arr, index_out, size);
  if(row != NULL) {
    for(uint out = 0; out< size; out++) {
      if(row[out].ind_out == index_in)
	return true;
    }
  }
  return false;
}


/**! Creates a bi-partite matching
   Assumes every key (the array index and the index stored as an elemnt in the array follows the same key-code-rule
*/
void id_simil_list::execInternLocalReciProc(rel_t **&arr, uint arr_size) {
  uint debug_cnt_index_in = 0, debug_cnt_index_out_ok = 0, debug_cnt_index_out_dicard = 0;
  for(uint index_in = 0; index_in< arr_size; index_in++) {
    debug_cnt_index_in++;
    uint size = 0;
    rel_t *row = getLocalRow(arr, index_in, size);
    if(row != NULL) {
      queue<rel_t> list;
      for(uint out = 0; out< size; out++) {
	//! Tests if 'index_in' has a match towards 'row[out].ind_out'
	const uint index_out = row[out].ind_out;
	if(hasRelation(arr, index_in, index_out)) {
	  debug_cnt_index_out_ok++;
	  list.push(row[out]); // ok; is added.
	} else { //! Was not added:
	  debug_cnt_index_out_dicard++;
	  if(DEBUG_PRINT_DISCARDED_PAIRS && listTaxa) {	    
	    const uint taxon_in = taxa::getTaxonId(index_in, listTaxa, taxon_length);
	    const uint taxon_out = taxa::getTaxonId(index_out, listTaxa, taxon_length);
 	    printf("(possible-ortholog-pair-discarded-due-to-non-reciproc)\tpair (%s, %s)\n", 
 		   listTaxa[taxon_in].getProteinNameWorld(index_in), 
 		   listTaxa[taxon_out].getProteinNameWorld(index_out));
	  }
	}
      }
      clearLocalData(arr[index_in]);
      while(!list.empty()) { // has data 
	//! Adds 'index_in' -> 'list.front().ind_out'
	insertLocalElement(arr, index_in, list.front().ind_out, list.front().distance);
	list.pop();
      }
    }
  }
}

/**
   @return the number of pairs
   @remarks Use a '+arr_size' offset to get the size, due to the first reference object.
**/
uint id_simil_list::get_total_number_of_pairs(rel_t **arr, uint arr_size) {
  uint sum = 0;
  if(arr != NULL) {
    for(uint world_in = 0; world_in < arr_size; world_in++) {
      if(arr[world_in] != NULL) {
	sum += get_number_of_pairs(arr[world_in]);
      }
    }
 }
 return sum;
}

/**
   @return the number of pairs
   @remarks Use a '+arr_size' offset to get the size, due to the first reference object.
**/
uint id_simil_list::get_total_number_of_pairs(rel_t **arr, uint arr_size, uint taxon_start, uint taxon_end, taxa *listTaxa, uint taxon_length, bool *list_of_this_node_taxa_responsilibties) {
#ifdef USE_MPI
  assert(list_of_this_node_taxa_responsilibties);
#endif
  uint sum = 0;
  assert(listTaxa);
  assert(taxon_start < taxon_length);
  assert(taxon_end <= taxon_length);
  const uint world_start_index = listTaxa[taxon_start].rel_start;  
  const uint world_end_index   =  listTaxa[taxon_end-1].rel_end;
  assert(world_start_index < arr_size);
  assert(world_end_index   <= arr_size);
  if(arr != NULL) {
    for(int taxon_id = taxon_start; taxon_id < (int)taxon_end; taxon_id++) {
      if(!list_of_this_node_taxa_responsilibties || list_of_this_node_taxa_responsilibties[taxon_id]) {
	for(int world_index_in =listTaxa[taxon_id].rel_start; world_index_in< listTaxa[taxon_id].rel_end;world_index_in++) {  	  
	  if(arr[world_index_in] != NULL) {
	    sum += get_number_of_pairs(arr[world_index_in]);
	  }
	}
      }
    }
 }
 return sum;
}

/**
   @return the sum of distances for pairs
   @remarks 
   # Used verifying that the content of two objects are equal
   # Practical when MPI sending/receving is to be verified.
**/
float id_simil_list::get_total_sum_of_distances_for_pairs(rel_t **arr, uint arr_size) {
  float sum = 0;
  float sum_debug = 0;
  uint cnt_relations = 0, debug_cnt_relations = 0;
 if(arr != NULL) {
    for(uint world_in = 0; world_in < arr_size; world_in++) {
      if(arr[world_in] != NULL) {
	uint size = 0;
	rel_t *row = id_simil_list::getLocalRow(arr, world_in, size);
	for(uint i = 0; i < size; i++) {
	  sum += row[i].distance;
	  cnt_relations++;
	}
#ifndef NDEBUG
	//! Uses an alternative procedure, and verfies that the result is the same:
	if(arr[world_in] != NULL) {
	  uint size_debug = arr[world_in][0].ind_out+1; // The number of elements
	  //	  if(size_debug) debug_cnt_relations += 1; // In order to get the meta-object
	  for(uint i = 1; i < size_debug; i++) {
	    sum_debug += arr[world_in][i].distance;
	    debug_cnt_relations++;
	  }	  
	}
	
	if(!(cnt_relations == debug_cnt_relations)) {
	  fprintf(stderr, "!!\t [%u]\tcnt_relations(%u) != debug_cnt_relations(%u), at line %d in file %s\n", world_in, cnt_relations, debug_cnt_relations, __LINE__, __FILE__);
	  assert(cnt_relations == debug_cnt_relations);
	}
	
	if(!((int)sum_debug == (int)sum)) {
	  fprintf(stderr, "!!\t [%u]\tsum_debug(%f) != sum(%f), at line %d in file %s\n", world_in, sum_debug, sum, __LINE__, __FILE__);
	  assert((int)sum_debug == (int)sum);
	}

#endif
      }
    }
 }
 return sum;
}
//! Prints the array givens as input argument
void id_simil_list::printLocalArr(FILE *f, rel_t **arr, const bool use_names, taxa_t *listTaxa, uint taxon_length, uint arr_size/*, bool to_file*/) {
  assert(f);
  fprintf(f,"#\t\tarr_size=%u (at line %d in file %s)\n", arr_size, __LINE__, __FILE__);
  if(arr != NULL) {
    for(uint world_in = 0; world_in < arr_size; world_in++) {
      if(arr[world_in] != NULL) {
	const uint taxon_in = taxa::getTaxonId(world_in, listTaxa, taxon_length);
	const uint protein_in = listTaxa[taxon_in].getLocalIndex(world_in);
	uint size = 0;
	rel_t *row = getLocalRow(arr, world_in, size);
	if(size > 0) {
	  if(use_names) fprintf(f,"-\t%s has %u orthologs:\t", listTaxa[taxon_in].arrKey[protein_in], arr[world_in][0].ind_out);
	  else fprintf(f,"-\t%u has %u orthologs:\t", world_in, arr[world_in][0].ind_out);
	  for(uint k = 0; k<size; k++) {
	    //if(true)fprintf(f,"arr[%u]= %u\n", k, arr[k]);
	    //	    const uint real_protein_out = row[k].ind_out;
	    const uint taxon_id = taxa::getTaxonId(row[k].ind_out, listTaxa, taxon_length);
	    listTaxa[taxon_id].printProtein(f, row[k].ind_out, use_names, '_');
	    fprintf(f,", ");
	    //	if(use_names) fprintf(f,"%s, ", arrKey[k]);	else fprintf(f,"%u, ", k);
	  }
	  fprintf(f,"\n");
	} else {
	  if(use_names) fprintf(f,"-\t%s has %u orthologs:\t", listTaxa[taxon_in].arrKey[protein_in], arr[world_in][0].ind_out);
	  else fprintf(f,"-\t%u has %u orthologs:\t", world_in, arr[world_in][0].ind_out);
	  fprintf(f,"\n");
	}
      } else fprintf(f,"!!\t'arr[world_in(%u)] == NULL\n", world_in);
    }
  } else fprintf(f,"(arr not set)\n");
}

#ifdef USE_MPI
#include "mpi_id_simil_list.h"
/**
   @brief Performs intra-node building of reciprocal ortholog pairs
   @remarks Core code located in file "mpi_id_simil_list.h", making the seperation of the strictly MPI code verbouse.   
**/
void id_simil_list::mpi_tx_rx_recip_tx_rx(int taxon_length) {
  class mpi_id_simil_list tx_rx = mpi_id_simil_list();
  assert(listTaxa);
  //  const uint total_number_of_proteins = listTaxa[taxon_length-1].total_cnt;
  const uint total_number_of_proteins = listTaxa[taxon_length-1].rel_end;
  //! Perform the the exchange, letting all nodes get a complete list of the orthologs.
  tx_rx.tx_rx_before_recip(this, arr, arr_size, total_number_of_proteins);
  tx_rx.free_memory();
}

  /**
     @brief Sends the co-orthologs to all the nodes.
     @remarks Core code located in file "mpi_id_simil_list.h", making the seperation of the strictly MPI code verbouse.   
  **/
void id_simil_list::mpi_send_co_orthologs_accross_nodes(taxa *listTaxa, int taxon_length, stack_rel *stackRel, list_file_parse<rel> *listStructData) {
  //! Perform the the exchange, letting all nodes get a complete list of the orthologs.
  class mpi_id_simil_list tx_rx = mpi_id_simil_list();
  tx_rx.send_co_orthologs(listTaxa, taxon_length, stackRel, listStructData);
  tx_rx.free_memory();
}

#endif
//! The main test function for this class  
void id_simil_list::assert_class(const bool print_info) {
  const static char *class_name = "id_simil_list";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  // TODO: Add something
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);

}

