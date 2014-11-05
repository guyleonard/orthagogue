#include "taxon_pair.h"
//! Sets the values for this object.
void taxon_pair::set_values(taxon_pair &obj) {
  taxon = obj.get_taxon();
  protein_start = obj.get_protein_start();
  protein_end = obj.get_protein_end();
  listParseData = obj.get_listParseData();
  listStructData = obj.get_listStructData();
  listTaxa = obj.get_listTaxa();
    
}

// Returns the variables
void taxon_pair::getVariables(uint &_taxon, uint &_protein_start, uint &_protein_end) {
  _taxon=taxon, _protein_end = protein_end, _protein_start = protein_start;
}

// Returns true if is has more data
bool taxon_pair::has_data() {return (protein_end != 0);}

//! @return an allocated region: A standard interface to ensure consistency
taxon_pair *taxon_pair::initList(uint size) {
  taxon_pair *pair = (taxon_pair*)malloc(sizeof(taxon_pair)*size);
  log_builder::test_memory_condition_and_if_not_abort(pair!=NULL, __LINE__, __FILE__, __FUNCTION__);
  memset(pair, 0, sizeof(taxon_pair)*size);
  for(uint i = 0; i < size; i++) pair[i] = taxon_pair();  // By default initializes
  return pair;
}

  
//! Prints the data
void taxon_pair::print()  {
  if(protein_end > 0) 
    printf("\ttaxon=%d, protein_start=%u, protein_end=%u\n", taxon, protein_start, protein_end);
}

//! Generates a file for a set of this class:
void taxon_pair::write_class_list(uint index, taxon_pair *pair, uint taxon_start, uint taxon_end,
				  uint taxon_length, uint n_threads, 
				  const bool only_inpa, list_file_parse_t *listParseData, list_file_struct_t *listStructData, taxa_t *listTaxa) {
#ifndef NDEBUG
  assert(taxon_start <= taxon_end);
  assert(taxon_end <= taxon_length);
  // Starts the time measurement in order to estimate the cost of this logging operation.
  struct tms tmsstart; clock_t clock_log;
  if((clock_log = times(&tmsstart)) == -1) // starting values
    fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
  bool create_new_file = true;
  if(listParseData && !only_inpa) create_new_file = true;
  else {create_new_file = false;}
  FILE *f_log = log_builder::get_log_file_pointer(n_threads, "thread_allocation_list", __FILE__, __LINE__, create_new_file);
  if(f_log) {
    if(listParseData) {
      if(!only_inpa) fprintf(f_log, "\n\n(1) Schedules for the Orthologs");
      else           fprintf(f_log, "\n\n(2) Schedules for the Inparalogs");
    } else fprintf(f_log, "\n\n(3) Schedules for the co-orthologs");
    fprintf(f_log, ", using %u threads, and iterating over the taxa ", n_threads);
    for(uint i = taxon_start; i < taxon_end; i++) fprintf(f_log, "%s, ", listTaxa[i].name);
    fprintf(f_log, " resulting in the blocks:\n");
    uint taxon_prev = 0;
    for(uint i = 0; i < index; i++) {
      if(taxon_prev != pair[i].taxon) {taxon_prev = pair[i].taxon; fprintf(f_log, "--------------------------------------------------\n");}
      fprintf(f_log, "\t%u\ttaxon_pair[%s=%d]", i, listTaxa[pair[i].taxon].name, pair[i].taxon);
      if(listParseData && only_inpa) fprintf(f_log, "[%s]", listTaxa[pair[i].taxon].name);
      else fprintf(f_log, "[(all)]");	  
      const uint cnt_proteins_allocated = (pair[i].protein_end - pair[i].protein_start);
      uint cnt_protein_pairs = 0;
      // if(listParseData) {	//	assert(listParseData->has_data(pair[i].taxon));
      // 	listParseData->has_data(pair[i].taxon);
      // }
      if(listParseData) cnt_protein_pairs = listParseData->getProteinLength(pair[i].taxon, only_inpa, pair[i].protein_start, pair[i].protein_end);
      else              cnt_protein_pairs = listStructData->getProteinLength(pair[i].taxon, only_inpa, pair[i].protein_start, pair[i].protein_end);
      //	  printf("cnt_protein_pairs=%u\n", cnt_protein_pairs);
      fprintf(f_log, "\t{starts_at(%d) and ends_in(%d)}\t", pair[i].protein_start, pair[i].protein_end);
      fprintf(f_log, " holding %u proteins, giving %u pairs to be processed.\n", cnt_proteins_allocated, cnt_protein_pairs); 
    }
    // Generates the estimate on the time consumption:
    struct tms tmsend; clock_t end; long clktck = 0;
    if((end = times(&tmsend)) == -1) fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
    if((clktck =  sysconf(_SC_CLK_TCK)) < 0) fprintf(stderr, "!!\tAn error in sys-configuration at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
    const clock_t t_real = end - clock_log;
    const double log_time = t_real/(double)clktck;
    fprintf(f_log, "--\tUsed in total %10.8f seconds for this operation\n", log_time);
    fclose(f_log);
  } 
#endif
}
/**
   @brief Builds the set of blocks to be used during the parsing.
   @param <taxon_start> The id of the first taxa to build the list for.
   @param <taxon_end> The id of the first taxa (and the following onwards) to not include in the list build.
   @param <taxon_length> The total number of taxa in the collection used as input.
   @param <only_inpa> If set, only calculates the size of its inparalogs.
   @param <from_parse> If set, uses data build from the parsing of the blast file.
   @param <listParseData> Holds data from the parsing operation.
   @param <listStructData> Holds data from the filtering operation.
   @param <listTaxa> Holds data about the taxa involved in this operation.
   @remarks 
   # Used by methods originating in class blast_filtering.
   # Depends on scheduling functions found in class list_file_parse.
   # Divid blocks using approximations of the best sizes of data.
   @author: Ole Kristian Ekseth (oekseth).
**/
taxon_pair *taxon_pair::init_taxon_pair(uint taxon_start, uint taxon_end, uint taxon_length
					, uint n_threads
					, const bool only_inpa,
					list_file_parse_t *listParseData, list_file_struct_t *&listStructData, taxa_t *listTaxa) {  
  assert(n_threads > 0);
  if(listTaxa != NULL) {
    uint data_length = 0;
    uint biggest_collection_size = 0;

    //! Get the length of the data to process.    
    if(listParseData) {
      data_length = listParseData->getTotalLengthOfData(taxon_start, taxon_end, only_inpa, biggest_collection_size);
      //printf("\t data_length=%u, listParseData=%p, listStructData=%p, at %s:%d\n", data_length, listParseData, listStructData, __FILE__, __LINE__);
    }
    else if(listStructData) {
      data_length = listStructData->getTotalLengthOfData(taxon_start, taxon_end, only_inpa, biggest_collection_size);
      //printf("\t data_length=%u, range[%u, %u], taxon_length=%u, listParseData=%p, listStructData=%p, at %s:%d\n", data_length, taxon_start, taxon_end, taxon_length, listParseData, listStructData, __FILE__, __LINE__);
    } else {log_builder::test_memory_condition_and_if_not_abort(!listStructData && !listParseData, __LINE__, __FILE__, __FUNCTION__);}
    //printf("\t data_length=%u, listParseData=%p, listStructData=%p, at %s:%d\n", data_length, listParseData, listStructData, __FILE__, __LINE__);    
    bool use_brute_force_count = false;
    if(!data_length) { // then we explcitily update it for the given taxa:      
      biggest_collection_size = 0;
      for(uint taxon_id = min(taxon_start, taxon_length); taxon_id < min(taxon_length, taxon_end); taxon_id++) {
	data_length += listTaxa[taxon_id].total_cnt;
	if((uint)listTaxa[taxon_id].total_cnt > (uint)biggest_collection_size) {
	  biggest_collection_size = listTaxa[taxon_id].total_cnt;
	}
	
      }
      use_brute_force_count = true;
    }

    if(data_length > 0) {
      //printf("\t data_length=%u, at %s:%d\n", data_length, __FILE__, __LINE__);
      //! Estimates preferable size for each chunk of data a thread shall process before "continuing":
      uint size_list = n_threads*2; // TODO: Consider updating this number iaw emprical experiences.
      const uint taxon_set_length = taxon_end-taxon_start;
      if(size_list < taxon_set_length) size_list = taxon_end-taxon_start+1;      
      uint avg_size = data_length / size_list;
      if(!avg_size) avg_size = 1;

      //! One of the main challenges is that the length of the pairs, defined for each inner taxa, differs a lot: Below is the worst case scenario included, thereby avoiding the ned for enlarging the list:
      const uint temp_size = taxon_set_length*(1+(biggest_collection_size / avg_size)); // the upper bound
      size_list = max(size_list, temp_size)+1; // In order to have an empty element stating the end of it.
      taxon_pair_t *pair = taxon_pair::initList(size_list);    
      //printf("\t size_list=%u, range-taxa=[%u, %u], at %s:%d\n", (uint)size_list, (uint)taxon_start, (uint)taxon_end, __FILE__, __LINE__);
      uint index = 0;  // the index to add data
      //      uint cnt_proteins_allocated = 0; // 
      for(uint taxon_id = taxon_start; taxon_id < taxon_end; taxon_id++) {
	//! Only processes if the given "inner" taxon is "myranks" reposponisiblity; when MPI is not used, it's always ;)
	bool taxon_is_our_responsibility = false;
	if(listParseData) taxon_is_our_responsibility = listParseData->taxon_is_myranks_responsibility(taxon_id);
	else if(listStructData) taxon_is_our_responsibility = listStructData->taxon_is_myranks_responsibility(taxon_id);
	if(taxon_is_our_responsibility) {
	  if(listTaxa[taxon_id].total_cnt > 0) {
	    // if(use_brute_force_count) {
	    //   //! Then, given our sparse data-sets, we only use a rough exsitamtion:
	    //   pair[index] = taxon_pair(taxon_id, start_pos, end_pos, listParseData, listStructData, listTaxa);
	    //   continue
	    // }
	    int start_pos = 0; // defining the start of the block
	    const uint max_cnt_proteins = (uint)((float)1.3*avg_size); //1 + listTaxa[taxon_id].total_cnt / n_threads;
	    bool values_changed = true;
	    //	  printf("-\t taxon_id(%u) has total_length(%d) in taxon_pair at line %d\n", taxon_id, listTaxa[taxon_id].total_cnt, __LINE__);

	    //! Divides the elements for a taxon into blocks, depending on the max number of proteins set:
	    while(start_pos < listTaxa[taxon_id].total_cnt && values_changed) { // Iterates through the data for the given inner taxa:
	      int end_pos = 0;
	      if(listParseData)       end_pos = listParseData->getProteinStartOfNextBuffer(start_pos, taxon_id, only_inpa || use_brute_force_count, max_cnt_proteins);
	      else if(listStructData) end_pos = listStructData->getProteinStartOfNextBuffer(start_pos, taxon_id, only_inpa || use_brute_force_count, max_cnt_proteins);
	      const int length_of_data = end_pos - start_pos;
	      if((length_of_data) > 0) {
		//		cnt_proteins_allocated += length_of_data; //(end_pos - start_pos);
		assert(index < size_list);
		if(listTaxa[taxon_id].total_cnt == (end_pos+1)) end_pos = end_pos+1; 
		pair[index] = taxon_pair(taxon_id, start_pos, end_pos, listParseData, listStructData, listTaxa);
		index++, start_pos = end_pos; // swaps
	      } else values_changed = false; // If no data is added.
	    }
	  }
	}
      }
#ifndef NDEBUG
      //! Asserts that all proteins are included in this collection:
      uint sum_inner = 0;
      for(uint taxon_id = taxon_start; taxon_id < taxon_end; taxon_id++) {
	bool taxon_is_our_responsibility = false;
	if(listParseData) taxon_is_our_responsibility = listParseData->taxon_is_myranks_responsibility(taxon_id);
	else if(listStructData) taxon_is_our_responsibility = listStructData->taxon_is_myranks_responsibility(taxon_id);
	if(taxon_is_our_responsibility) {	  
	  sum_inner += listTaxa[taxon_id].total_cnt;
	}
      }
      uint taxon_in = 0, protein_end=0, protein_start = 0; 
      uint bucket_cnt = 0;
      //      printf("at line %d in file %s\n", __LINE__, __FILE__);
      uint diff_total = 0;
      do {
	pair[bucket_cnt++].getVariables(taxon_in, protein_start, protein_end); // The dataset to work ok
	if(protein_end != 0) {
	  const uint diff = protein_end - protein_start;
	  diff_total += diff;	  
	}
      } while(protein_end != 0); 
      //      printf("--\tsum_inner(%u), diff_total(%u) at line %d in file %s\n", sum_inner, diff_total, __LINE__, __FILE__);
      assert(diff_total == sum_inner);      
#endif
#ifdef LOG_WRITE_SCHEDULING_LIST
      taxon_pair::write_class_list(index, pair, taxon_start, taxon_end, taxon_length, n_threads, only_inpa, listParseData, listStructData, listTaxa);
#endif
      return pair;
    } else {      
      return NULL;
    }
  } else {
    assert(listTaxa);
    return NULL;
  }
}


taxon_pair *taxon_pair::get_taxon_pairs(taxon_pair *listPair, uint &list_pair_pos) {
  // Only do the oepration iof it has any data:
  //printf("\t at %s:%d\n", __FILE__, __LINE__);
  if(listPair && listPair[list_pair_pos].has_data()) {
    // First find the number of elements we need:
    uint start = list_pair_pos;
    const uint minimum_number = 10000; // The minimum number of pairs to require.
    uint elements_found = 0; 
    while(listPair[start].has_data() && (elements_found < minimum_number)) {
      elements_found += listPair[start].get_length();
      start++;
    }
    const uint lst_length =  start - list_pair_pos;
    //printf("\t lst_length=%u, at %s:%d\n", lst_length, __FILE__, __LINE__);
    if(lst_length) {
      taxon_pair *lst = (taxon_pair*)malloc(sizeof(taxon_pair)*(lst_length+1));
      log_builder::test_memory_condition_and_if_not_abort(lst!=NULL, __LINE__, __FILE__, __FUNCTION__);
      for(uint i = 0; i < lst_length; i++) {
	lst[i].set_values(listPair[list_pair_pos+i]);	
      }
      lst[lst_length] = taxon_pair(); // Sets the end.
      list_pair_pos += lst_length; // Updates the given variable.
      return lst;
    } else {
      return NULL;
    }
  } else {
    return NULL;
  }
}

taxon_pair::taxon_pair(uint _taxon, uint _protein_start, uint _protein_end,
		       list_file_parse_t *&_listParseData, list_file_struct_t *&_listStructData, taxa_t *_list) :
  taxon(_taxon), protein_start(_protein_start), protein_end(_protein_end),
  listParseData(_listParseData),   listStructData(_listStructData), listTaxa(_list) 
{}
taxon_pair::taxon_pair() : taxon(0), protein_start(0), protein_end(0), listParseData(NULL) {}
 
void taxon_pair::assert_class(const bool print_info) {
  const static char *class_name = "taxon_pair";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
#ifndef NDEBUG
  taxon_pair temp_0 = taxon_pair();
  assert(!temp_0.has_data());
  list_file_parse_t *listParseData = NULL;
  list_file_struct_t *listStructData = NULL;
  taxa_t *listTaxa = NULL;
  
  taxon_pair temp(0, 1, 2, listParseData, listStructData, listTaxa);
  assert(temp.has_data());
  uint t, ps, pe;
  temp.getVariables(t, ps, pe);
  assert(0 == t);
  assert(1 == ps);
  assert(2 == pe);
#endif
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}


