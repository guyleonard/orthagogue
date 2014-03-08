#include "pipe_parse_parse.h"
  
//! Returns the number of newlines (rows) found in the file:
loint pipe_parse_parse::debug_get_cnt_newlines(char *buffer_one, char *logical_end) {    
  assert(buffer_one);
  assert(logical_end);
  //  assert(buffer_one < (logical_end-2));
  assert(buffer_one <= logical_end);
  loint cnt_newlines = 0;
  while (buffer_one < (logical_end-2)) { 
    buffer_one += 2; // If two newlines after each other, drops them.
    char *line_start = buffer_one; // A point of reference.
    buffer_one = strchr(buffer_one, '\n');       // update: This implies that if a blink line is added, no damage is done.
    if(!buffer_one) buffer_one = logical_end;
    else {
      // Requires that the line must consist of enough chars to be consideres a real line:
      if((logical_end-line_start)>50) {
	cnt_newlines++;
      } 

    }
  }
  return cnt_newlines;
}

//! Initialises the data structures with the params given.
void pipe_parse_parse::set_parse_blocks_second_read(parse_read_blocks_t *parseBlocks_) {
  parseBlocks = parseBlocks_;
  if(parseBlocks)
    parseBlocks->init_thread_file_reading(CPU_TOT, FILE_INPUT_NAME);
}


//! Returns the arrNorm strcture
list_norm *pipe_parse_parse::getArrNorm(float max_input_value) {
  if(USE_EVERYREL_AS_ARRNORM_BASIS) {
    if(local_arrNorm) {
      const uint index_root = 0; // The part of the index to 'return, ie, to give a way and set to NULL, thereby not removing it when doing the 'cleaning' of the list_norm object.
      local_arrNorm[index_root]->merge_basis_zero(index_root, CPU_TOT, local_arrNorm, max_input_value);
      class list_norm *element = local_arrNorm[index_root];
      local_arrNorm[index_root] = NULL; // Avoid cleaning of it.
      return element;
    } else {
      fprintf(stderr, "!!\tAn error occured at line %d in function %s for class pipe_parse_parse: local_arrNorm not set. Contact oekseth@gmail.com if this error message is seen!\n", __LINE__, __FUNCTION__);
      return NULL;
    }
  } else return NULL;
}


void swap(int &a, int &b) {const int temp = b; b = a; a = temp;} 
/**
   sort_array -- Sorts the array by moving elements backward
   until the list at index i is sorted.
   @Changed: 22.06.2011 by oekseth
*/
void bubble_sort_array(int *keys, const int *values, const uint size) {
  for(uint i = 0; i < size; i++) {
    for(uint j = size-1; j > i; j--) {      
      if(values[keys[j]] < values[keys[j-1]]) swap(keys[j], keys[j-1]);
    }
  }
}
/**
   Returns the sorted list of the taxa: Alogrithm is slow but task is not critical.
   @Return: array[%d] = <global taxon index>
   @algortihm:
*/
int *getLowest_taxon(uint taxon_length, taxon_list_t listProteins) {
  int *array = (int*)malloc(sizeof(int)*taxon_length);
  int *values = (int*)malloc(sizeof(int)*taxon_length);
  memset(array, 0, taxon_length);   memset(values, 0, taxon_length);
  for(uint i = 0; i < taxon_length; i++) array[i] = i, values[i] = listProteins.getProteinLength(i);
  bubble_sort_array(array, values, taxon_length);
  for(uint i = 1; i < taxon_length; i++) {
    assert(values[array[i-1]] <= values[array[i]]);
  }
  free(values); values = NULL;
  return array;
}

/**! Initiates the hash list, whom is used for retrieving the unique id's for
   the proteins, given as a string of chars, in the input file
   @ProgramPoint: after the first run through the blast file:
   Called from 'initSecondRead'
*/
extern char *FILE_INPUT_NAME;//"all.blast"; // the name (the whole path) of the input file
void pipe_parse_parse::initHashProtein(taxon_list_t *listProteins) {
 if(listProteins->getLength() > 0) {
    taxon_length = listProteins->getLength(); // updates the global constant
    int *sorted_ind_list = getLowest_taxon(taxon_length, listProteins[0]);
    hashProtein = prot_list::init(taxon_length); // Only allocates, but do not initializes.
    for(int taxon_id = 0; taxon_id<taxon_length;taxon_id++) {
      const int hash_index = sorted_ind_list[taxon_id];
      hashProtein[taxon_id] = prot_list(listProteins->getProteinLength((uint)hash_index), // the number of proteins found for the taxon 
					listProteins->getBuffer(hash_index), // The list of protein names for the taxon
					listProteins->getTaxonName(hash_index), // The name of the taxon
					//					listProteins->getName(hash_index), // The name of the taxon
					true); // uses a reciprocal test
    }
    free(sorted_ind_list); sorted_ind_list=NULL;
  } else {
    fprintf(stderr, "!!\tOrthaGogue was not able to find any taxa listed in your input file. The authors guess it's incorrect column-numbering of the taxon- and protein location. In our software the first column is numbered '0'. For some users this is unusual, for others it's not.\n");
    fprintf(stderr, "!!\t- A Tips is to verify the input specifications: During the parsing, the program assumed that the taxon column is '%d' and the protein column is '%d', with a '%c' as a seperator between the taxon and protein label. These settings were given for the file '%s'.\n", INDEX_IN_FILE_FOR_TAXON_NAME, INDEX_IN_FILE_FOR_PROTEIN_NAME, SEPERATOR, FILE_INPUT_NAME);
    fprintf(stderr, "!!\t- As this error is fatal, the program exits: We would be thankful if you contact the authors if this problem is seen, sending an email to %s.\n", "oekseth@gmail.com"); exit(2);
  }

#ifndef NDEBUG
  if(stringBuffer) {
    const loint difference = (total_number_of_chars_processed_in_first_read - stringBuffer->getTotalLengthOfData());
    if(difference) {
      fprintf(stderr, "..--\tOur first_parse found %lld chars, while %lld lines was stored in the stringBuffer (difference = %lld) for file_input='%s'at line %d in method %s located in file %s\n",
	      total_number_of_chars_processed_in_first_read, stringBuffer->getTotalLengthOfData(), difference, FILE_INPUT_NAME, __LINE__, __FUNCTION__, __FILE__); 
      assert(difference == 0);
    }
  }
#endif
}

/**! Initiates data used for the second read
   @ProgramPoint: Called from the outside before initiating the pipe
   @Return: the 'listTaxa' to be used for the global purpose
*/
taxa_t* pipe_parse_parse::initSecondRead(taxon_list_t *listProteins, int &_taxon_length, int updated_cpu_cnt) {
  if(CPU_TOT != updated_cpu_cnt) {
    if(in_use) { free(in_use); } //   delete [] in_use;}
    in_use = (bool*)malloc(sizeof(bool)*updated_cpu_cnt);
    for(int i = 0;i<updated_cpu_cnt;i++) in_use[i] = false;
    CPU_TOT = updated_cpu_cnt; // Updates the variable itself.
  }
  assert(CPU_TOT > 0);
  initHashProtein(listProteins); // a function doing this
  _taxon_length = listProteins->getLength(); 
  FIRST_READ = false;
  first_protein = protein_relation::init(CPU_TOT);
  send_all_data_to_pipe = false;
  //  8912 olekrie   25   0 6234m 6.1g 1648 R 99.5  4.8   1:46.08 orthAgogue                                                                                                                                             
  proteinVector = protein_vector::init(CPU_TOT);
  max_sim_value = (float*)malloc(sizeof(float)*CPU_TOT);
  for(uint i = 0; i < (uint)CPU_TOT; i++) max_sim_value[i] = 0;
  //  9366 olekrie   25   0 10.6g  10g 1648 R 97.8  8.4   1:41.92 orthAgogue                                                                                                                                             
  if(USE_EVERYREL_AS_ARRNORM_BASIS) { 
    local_arrNorm = list_norm::init_list(CPU_TOT, taxon_length, 0, PRINT_NORMALIXATION_BASIS, DEBUG_NORM);  // Initiates the arrNorm  
  }
  // The below cal does not increase the memory signature in any major amount.
  listTaxa = intializeListTaxa();
  parseData = list_file_parse<p_rel>::init(CPU_TOT, hashProtein, taxon_length, listTaxa, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE);
  hashTaxa =  prot_list::init(listTaxa, taxon_length);
  return listTaxa;
}

/**! Initiates the global 'taxa_t listTaxa' structure
   @ProgramPoint: Calleed when, after the first run through the file, the pipe is completed
   @Return: the global lsitTaxa to be used in the processes to follow
*/
taxa_t* pipe_parse_parse::intializeListTaxa() {
  taxa_t *listTaxa = taxa::init(taxon_length);//(taxa_t*)malloc(sizeof(taxa)*taxon_length);
  listTaxa[0] = taxa(hashProtein[0].getName(),
		     hashProtein[0].getLength(), // the total number of proteins from this taxon
		     0,
		     hashProtein[0].getLength(), // the boundary (+1) of the range
		     hashProtein[0].arrKey
		     , SEPERATOR, INDEX_IN_FILE_FOR_PROTEIN_NAME, INDEX_IN_FILE_FOR_TAXON_NAME
		     );
  hashProtein[0].lock_variables_transfered_to_taxa_list();
  for (int i = 1;i<taxon_length;i++) {
    listTaxa[i] = taxa(hashProtein[i].getName(),
		       hashProtein[i].getLength(), // the total number of proteins from this taxon
		       listTaxa[i-1].rel_end,//		     0,
		       listTaxa[i-1].rel_end + hashProtein[i].getLength(), // the boundary (+1) of the range
		       hashProtein[i].arrKey
		       , SEPERATOR, INDEX_IN_FILE_FOR_PROTEIN_NAME, INDEX_IN_FILE_FOR_TAXON_NAME
		       );
    hashProtein[i].lock_variables_transfered_to_taxa_list();
  }
  return listTaxa;
}


/**
   @brief If NDEBUG is not set, asserts the parsing operation done.
   @return the number of elements combined for both input argument, and internally stored data.
   @remarks Generating of data overview de-activated if in NDEBUG mode
**/
loint pipe_parse_parse::assert_parsing(list_file_parse_t *data) {
#ifndef NDEBUG 
  // Writes a log file:
#ifdef USE_MPI
  int myrank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);   // Gets specific data about the reading:
  FILE *f_log = log_builder::get_log_file_pointer(CPU_TOT, "pipe_parse_parse", __FILE__, __LINE__, true, myrank);
#else
  FILE *f_log = log_builder::get_log_file_pointer(CPU_TOT, "pipe_parse_parse", __FILE__, __LINE__, true);
#endif
  if(f_log) {
    fprintf(f_log, "-\tFor file %s:\n", FILE_INPUT_NAME);
    fprintf(f_log, "#\ttotal_number_of_chars_processed_in_first_read = %lld\n", total_number_of_chars_processed_in_first_read);
    fprintf(f_log, "#\ttotal_number_of_chars_processed_in_second_read = %lld\n", total_number_of_chars_processed_in_second_read);
    fprintf(f_log, "#\tpairs_in_file_found = %lld\n", pairs_in_file_found);
    fprintf(f_log, "#\ttotal_number_of_pairs_overlapping = %lld\n", total_number_of_pairs_overlapping);
    fclose(f_log);
  }
  bool internal_has_data = false;
  if(parseData) {
    for(uint i = 0; i < (uint)CPU_TOT && !internal_has_data; i++) {
    //    for(uint i = 0; i < (uint)CPU_TOT; i++) {
      if(parseData[i] && parseData[i]->has_data()) internal_has_data = true;
    }
  }
  if((data && data->has_data()) || internal_has_data) {
    // Verifies that all the data is in the structure:
    total_number_of_chars_processed_in_second_read++; // If the last block, the trailing char did not exists, therefore adjusts it backward again.
    if(data) {
      if(!(pairs_in_file_sent_to_merge == data->getTotalLengthOfData())) {
	fprintf(stderr, "pairs_in_file_sent_to_merge(%d) =! data->getTotalLengthOfData(%d)",
		(int)pairs_in_file_sent_to_merge, (int)data->getTotalLengthOfData());
	assert(pairs_in_file_sent_to_merge == data->getTotalLengthOfData());
      }
    } else assert(pairs_in_file_sent_to_merge == 0);

    // At this point counts hte number of elements found in the structure located in the pipe_parse_parse object:
    loint sum_cnt_internal_object = 0;
    if(parseData) {
      for(uint i = 0; i < (uint)CPU_TOT; i++) {
	if(parseData[i]) sum_cnt_internal_object += parseData[i]->getTotalLengthOfData();
      }
      if(!((sum_cnt_internal_object + pairs_in_file_sent_to_merge) == pairs_in_file_found)) {
	fprintf(stderr, "In 'this' object, we have [pairs_in_file_sent_to_merge(%d) + located_locally(%llu)] != pairs_in_file_found(%d)\n",
		(int)pairs_in_file_sent_to_merge, sum_cnt_internal_object, (int)pairs_in_file_found);
	assert((sum_cnt_internal_object + pairs_in_file_sent_to_merge) == pairs_in_file_found);
      }
    } else {
      if(!(pairs_in_file_sent_to_merge == pairs_in_file_found)) {
	fprintf(stderr, "In 'this' object, we have pairs_in_file_sent_to_merge(%d) != pairs_in_file_found(%d)\n",
		(int)pairs_in_file_sent_to_merge, (int)pairs_in_file_found);
	assert(pairs_in_file_sent_to_merge == pairs_in_file_found);
      }
    }
    // Verifies that all data have been read:
    if((!parseBlocks) && stringBuffer) {
      if(!stringBuffer->has_data()) { // At this point the buffer should be empty
	fprintf(stderr, "!!\tAll the data was not read in the stringBuffer. (This error was found at line %d in method %s location in file %s.)\n",
		__LINE__, __FUNCTION__, __FILE__);
	assert(!stringBuffer->has_data()); // At this point the buffer should be empty
      }
    }

    //
    // Finds the number of lines that should have been read, using a slow, but easy-to-verify procedure:
    loint cnt_newlines = lines_in_file_found_first_read;
    loint cnt_newlines_real = 0;
    FILE *file = fopen(FILE_INPUT_NAME, "r");
    uint cnt_chars = 0; // In order to break when the file length is read.
    if(file) {
      if(reading_file_start_position) fseek(file, reading_file_start_position, SEEK_SET);
      char c = '_';
      do { 
	c = fgetc (file);
	if (c == '\n') {
	  cnt_newlines_real++;
	} 
	cnt_chars++;
      } while (c != EOF && (cnt_chars < reading_file_length));
      fclose (file);
    } else fprintf(stderr, "!!\tUnable to open file %s at line %d in method %s located in file %s. Please contact the developer at oekseth@gmail.com\n",
		   FILE_INPUT_NAME, __LINE__, __FUNCTION__, __FILE__);
    
    if(cnt_newlines_real != lines_in_file_found_last_read) { // Prints an informative error if they are not equal:
      fprintf(stderr, "!!\tThe file contained %lld newlines, while our identified only found %lld lines (difference = %lld) for file_input='%s'at line %d in method %s located in file %s. Please contact the developer at oekseth@gmail.com\n",
	      cnt_newlines_real, lines_in_file_found_last_read, (cnt_newlines_real - lines_in_file_found_last_read), FILE_INPUT_NAME, __LINE__, __FUNCTION__, __FILE__);
      //! Note: Is not a critical error, as MPI-blast often generates errnous files.
      //      assert(cnt_newlines_real == lines_in_file_found_last_read);
    } 

    // Verifyies that the same number of liens are read in both the first- and the last read::
    if(lines_in_file_found_first_read != lines_in_file_found_last_read) {      
      fprintf(stderr, "!!\t Read %llu lines in the first read, while %llu lines in the last read, implying a difference of %d lines: in total there should have been %llu newlines in the file.. (This error was found at line %d in method %s location in file %s.)\n",
	      lines_in_file_found_first_read, lines_in_file_found_last_read, (int)(lines_in_file_found_first_read - lines_in_file_found_last_read), cnt_newlines, __LINE__, __FUNCTION__, __FILE__);
      //! Note: Is not a critical error, as MPI-blast often generates errnous files.
      // assert(lines_in_file_found_first_read == lines_in_file_found_last_read);
    }

    //
    // Control that the whole file was prcessed, comparing chars read to the file size,
    // verifying that all the data in the file have been read:
    loint real_file_size = reading_file_length;
    if(!real_file_size) real_file_size = buffer_string_list::get_file_size(FILE_INPUT_NAME);
    const int difference = (int)(real_file_size - total_number_of_chars_processed_in_second_read);
    const int allowable_difference = 10;
    if((difference<-1*allowable_difference) || (difference>allowable_difference) ) {
      fprintf(stderr, "!\tdifference(%d) \treal_file_size(%lld) != total_number_of_chars_processed_in_second_read(%lld) with a difference(%lld). (This error is found at line %d in file %s. If seen, please contact the developer at oekseth@gmail.com.)\n",
	      difference, real_file_size, total_number_of_chars_processed_in_second_read, (real_file_size-  total_number_of_chars_processed_in_second_read),
	      __LINE__, __FILE__);
      assert(real_file_size == total_number_of_chars_processed_in_second_read);
    }
    const int difference_parsing = (int)(total_number_of_chars_processed_in_first_read - total_number_of_chars_processed_in_second_read);
    if(difference_parsing<-1*allowable_difference || difference_parsing>allowable_difference) {
      fprintf(stderr, "!!\ttotal_number_of_chars_processed_in_first_read(%lld) != total_number_of_chars_processed_in_second_read(%lld) with a difference(%lld). (This error is found at line %d in file %s. If seen, please contact the developer at oekseth@gmail.com.)\n",
	      total_number_of_chars_processed_in_first_read, total_number_of_chars_processed_in_second_read, (total_number_of_chars_processed_in_first_read-  total_number_of_chars_processed_in_second_read),
	      __LINE__, __FILE__);
      assert(total_number_of_chars_processed_in_first_read == total_number_of_chars_processed_in_second_read);
    }

    loint tot_elements = 0;
    if(data) tot_elements += data->getTotalLengthOfData();
    if(parseData) {
      for(uint i = 0; i < (uint)CPU_TOT; i++) {
	tot_elements += parseData[i]->getTotalLengthOfData();
      }
    }
    if(tot_elements != pairs_in_file_found) {
      printf("!!\ttot_elements(%d) != pairs_in_file_found(%d) in pipe_parse_parse.cxx \n", (int)tot_elements, (int)pairs_in_file_found);
      assert(tot_elements == pairs_in_file_found);
    }
    return tot_elements;
  }
#endif
  return 0;
}
/**! Copies the content of 'this' into the argument, as the contains much
   more data than the elements of 'this, i.e. data >> list[..]
   - As the name indicates, it in addtion frees the memory allocated
   specifially for this class.
   @Changed: 29.12.2010 by oekseth.
*/
list_file_parse_t *pipe_parse_parse::free_memory(list_file_parse_t *data) {  
  const loint tot_elements = assert_parsing(data); // Asserts the parsing.
  if(parseData) {
    if(data) {
      for(uint i = 0; i < (uint)CPU_TOT; i++) {
	protein_relation empty = protein_relation();
	uint cnt_elements_in_all_relation_lists = 0;
	data->merge_data(parseData[i], empty, cnt_elements_in_all_relation_lists); // Note that the function called only has one parameter set
	list_file_parse<p_rel>::close(parseData[i], false); // Ensures safeness, but may be commented
      }
    } else { // input not set:
      for(uint i = 0; i < (uint)CPU_TOT; i++) {
	if(!data) { // The first element found is used as a basis
	  data = parseData[i];	  
	} else { // When the first element is found, the rest follows the same standard procedure
	  protein_relation empty = protein_relation();
	  uint cnt_elements_in_all_relation_lists = 0;
	  data->merge_data(parseData[i], empty, cnt_elements_in_all_relation_lists); // Note that the function called only has one parameter set
	  list_file_parse<p_rel>::close(parseData[i], false); // Ensures safeness, but may be commented
	}
      }
    }
  }
#ifndef NDEBUG
  if(data && data->has_data()) assert(tot_elements == data->getTotalLengthOfData()); // The 'root' should now consist of the sum of elements found before the merging.
#endif

  //  if(parseBlocks && data && data->has_data()) {
  bool internal_has_data = false;
  if(parseData) {
    for(uint i = 0; i < (uint)CPU_TOT && !internal_has_data; i++) {
      if(parseData[i] && parseData[i]->has_data()) internal_has_data = true;
    }
  }
  if((data && data->has_data()) || internal_has_data) {
    //    if(parseBlocks) {
    if(true) {
      assert(parseBlocks->assert_file_is_read());
    }
  }
  if(stringBuffer) {stringBuffer->free_mem(); delete stringBuffer; stringBuffer = NULL;}
  if(listProteins) {free(listProteins);listProteins= NULL;}
  // Do not remove the content of the class, but mere the pointer to it:
  list_file_parse<p_rel>::close_only_pointer(parseData); 
  protein_vector::free_list(proteinVector, (uint)CPU_TOT);
  proteinVector = NULL;
  free(in_use); in_use = NULL;
  free(first_protein); first_protein = NULL;
  free(max_sim_value); max_sim_value = NULL;
  list_norm::close(local_arrNorm, CPU_TOT);
  if((data && data->has_data()) || internal_has_data) {
    prot_list::close(hashTaxa); // only one elemnt in the list, as only one ontolgy(taxa) is represented.
    prot_list::close(hashProtein, taxon_length);
  }
  return data;
}


//-----------


//! Gets the overlap data, inserts it into the strcture, and returns an updated position in the string
char *pipe_parse_parse::getOverlapAndUpdate(Parse &p, char *&pos_column_start) {
  pos_column_start++; // jumps beyound the tab
  char *pos_column_end = strchr(pos_column_start, '\t'); // After the 3. column
  pos_column_end++; // jumps over the 'tab'
  if(use_improved_overlap_algo) {
    // jumps over column 4, 5, 6:
    assert(pos_column_end);  pos_column_end = strchr(pos_column_end, '\t'); // Now situated at the start of the 5. column
    pos_column_end++; // jumps over the 'tab'
    assert(pos_column_end);  pos_column_end = strchr(pos_column_end, '\t'); // Now situated at the start of the 6. column
    pos_column_end++; // jumps over the 'tab'
    assert(pos_column_end);  pos_column_end = strchr(pos_column_end, '\t'); // Now situated at the start of the 7. column
    pos_column_end++; // jumps over the 'tab'
    const int val_in_first = atoi(pos_column_end);
    assert(pos_column_end);  pos_column_end = strchr(pos_column_end, '\t'); // Now situated at the start of the 8. column
    pos_column_end++; // jumps over the 'tab'
    const int val_in_second = atoi(pos_column_end);
    p.overlap_in = val_in_second - val_in_first; // Sets the overlap for the inner
    assert(pos_column_end);  pos_column_end = strchr(pos_column_end, '\t'); // Now situated at the start of the 9. column
    pos_column_end++; // jumps over the 'tab'
    const int val_out_first = atoi(pos_column_end);
    pos_column_end = strchr(pos_column_end, '\t'); // Now situated at the start of the 10. column
    assert(pos_column_end);
    pos_column_end++; // jumps over the 'tab'
    const int val_out_second = atoi(pos_column_end);
    p.overlap_out = val_out_second - val_out_first; // Sets the overlap for the outer
    //    printf("(%d - %d) && (%d - %d)  ", val_in_second, val_in_first, val_out_second, val_out_first);
  } else p.overlap_in = atoi(pos_column_end);

  return pos_column_end;
}


//! Gets the overlap data, inserts it into the strcture, and returns an updated position in the string
char *pipe_parse_parse::getDistanceUpdate(Parse &p, char *&pos_column_start, float &max_sim_val, char *logical_end, const bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE) {
  // Jumps past 6 Columns:
  char *start = pos_column_start;
  if(!use_improved_overlap_algo) {
    assert(pos_column_start);    pos_column_start = strchr(pos_column_start, '\t'); pos_column_start++;//column to jump over is 1
    assert(pos_column_start);    pos_column_start = strchr(pos_column_start, '\t'); pos_column_start++;//column to jump over is 2
    assert(pos_column_start);    pos_column_start = strchr(pos_column_start, '\t'); pos_column_start++;//column to jump over is 3
    assert(pos_column_start);    pos_column_start = strchr(pos_column_start, '\t'); pos_column_start++;//column to jump over is 4
    assert(pos_column_start);    pos_column_start = strchr(pos_column_start, '\t'); pos_column_start++;//column to jump over is 5
    assert(pos_column_start);    pos_column_start = strchr(pos_column_start, '\t'); pos_column_start++;//column to jump over is 6
    assert(pos_column_start);
  }
  char *before_problem = pos_column_start;
  pos_column_start = strchr(pos_column_start, '\t'); pos_column_start++;//column to jump over is 7
  if(pos_column_start > logical_end) {
    log_builder::throw_warning(blastp_syntax, __LINE__, __FILE__, __FUNCTION__, "Column in blast file not found while calculating the distance in the blast file");
    blast_extractors_t::print_segment(start, before_problem);
    blast_extractors_t::print_segment_stepwise(start, before_problem);
    return NULL;
  } else if (pos_column_start == logical_end) {
    return NULL;
  } else {
    if(pos_column_start != NULL) {
      if(USE_LAST_BLAST_CLOMUN_AS_DISTANCE) { // uses the last blast column as input
	pos_column_start = strchr(pos_column_start, '\t'); pos_column_start++;//column to jump over is 8
	assert(pos_column_start);
	p.distance = atof(pos_column_start);
	if(p.distance > max_sim_val) {
	  max_sim_val = p.distance;
	}
      } else { // uses the second last blast column as input
	// SIMILARITY or DISTANCE:
	double temp_dist = atof(pos_column_start);
	if (temp_dist == 0) p.distance = 0;
	else {
	  p.distance = -1 * log10(temp_dist);
	  if(p.distance > max_sim_val) {
	    max_sim_val = p.distance;
	  }
	}
      }
    }
    return pos_column_start;
  }
}

//! Returns the updated line end of the current line read
char *getUpdatedLineEndPosition(char *&buffer_start, char *logical_end) {
  buffer_start = strchr(buffer_start, '\n');//line_end;//strchr(buffer_one, '\n');
  if(buffer_start == NULL) {
    return logical_end;
  } else {
    return (buffer_start +1);
  } 
}

/**
   @brief Inserts the protein pairs into data structures.
   @param <my_id> The id of the given thread.
   @param <buffer_one> The start position in memory to retrieve data from.
   @param <logical_end> The end position in memory to retrieve data from.
   @param <lines_in_file_found> Sets the number of lines found.     
   @return The total number of proteins found
   @remarks The steps taken for each row in the input given (ie, the algorithm) is:
   -# Jumps to next line if (a) the line does not contain the two identifying labels or (b) the labels found does not match any found in the previous run through the file (ie,  it's not given as a left-right relation)
   -# Get the 4 indexes for the pair found (taxa and protein label).
   -# Updates the data structure: If the relation is new inserts it, else updates the values given.
   @author Ole Kristian Ekseth (oekseth)
**/
uint pipe_parse_parse::parse_blast_blocks_data(int my_id, char *buffer_one, char *logical_end, loint &lines_in_file_found, loint &cnt_overlapping_pairs) {
  // Assumptions:
  assert(parseData);
  assert(parseData[my_id]);
  assert(buffer_one);
  assert(logical_end);
  assert(buffer_one < (logical_end-1));
  loint remaining_chars = (loint)(logical_end - buffer_one);
  uint cnt_inserted_pairs = 0;
  while (buffer_one < (logical_end)) { 
    Parse p = Parse();
    char *start_row = buffer_one;
    assert(buffer_one != logical_end);
    assert(buffer_one < logical_end);
    assert(buffer_one != NULL);
    if(blast_extractors_t::getParseHeaders(p, buffer_one, logical_end, blast_settings))	{
      // Requires that the line must consist of enough chars to be consideres a real line:
      if((logical_end-start_row)>50) {
	lines_in_file_found++;
      }
      buffer_one = getOverlapAndUpdate(p, buffer_one);
      assert(buffer_one);
      buffer_one = getDistanceUpdate(
				     p,
				     buffer_one,
				     max_sim_value[my_id],
				     logical_end, USE_LAST_BLAST_CLOMUN_AS_DISTANCE);
      if(!buffer_one) {
	blast_extractors::print_segment(start_row, logical_end);
	buffer_one = logical_end;
      } else {
	buffer_one = getUpdatedLineEndPosition(buffer_one, logical_end);
      }
      
      bool overlap_is_inserted = false;
      protein_relation prot = proteinVector[my_id].get_protein_indexes(p, overlap_is_inserted, listTaxa); // Gets the indexes for the row
      if(prot.exsists()) { // Must have labels recognised int the first parsing:
	const bool prot_pair_seen_before = first_protein[my_id].is_set();
	if(!prot_pair_seen_before) { // Not set, and thereby the first protein
	  first_protein[my_id].copy(prot);
	} 
	bool is_a_new_protein = false;
	if(parseData[my_id]->insert_new_rel(prot.taxon_in, prot.taxon_out, p, prot)) {
	  cnt_inserted_pairs++;
	  is_a_new_protein = true;
	} else {cnt_overlapping_pairs++;}	
	if(USE_EVERYREL_AS_ARRNORM_BASIS) {
	  assert(local_arrNorm[my_id]);
	  local_arrNorm[my_id]->insert((uint)prot.taxon_in, (uint)prot.taxon_out, p.distance, prot.protein_in, prot.protein_out, listTaxa, !is_a_new_protein);	  
	}
      } 
    } else { // The header is incorrect: Moves to the line end, in order to start the next read
      buffer_one = strchr(buffer_one, '\n');       // update: This implies that if a blink line is added, no damage is done.
    }
    if(buffer_one && (buffer_one < logical_end)) {
      if(*buffer_one =='\n') {buffer_one++;}
    } else buffer_one = logical_end;
    const loint remaining_chars_current = (loint)(logical_end - buffer_one);
    assert(!(remaining_chars_current > remaining_chars));
    p.free_memory();
  }
  return cnt_inserted_pairs;
} 


//! Find the protein identifikators
// Returns the last position before a new protein is found: used to set the block sizes in the enxt run
mem_loc pipe_parse_parse::parse_blast_blocks_ids(int my_id, char *buffer_one, char *logical_end, loint &lines_in_file_found) {
  char *last_block_line = NULL;
  char *buffer_start = buffer_one;
  //  uint cnt_liens = 0;

  assert(buffer_one);
  assert(logical_end);
  assert(buffer_one <= logical_end);
  loint remaining_chars = (loint)(logical_end - buffer_one);
  while (buffer_one < logical_end) {
    Parse p = Parse();
    char *line_start = buffer_one;
    assert(line_start);
    char *temp = blast_extractors::getIDColumn(true, p, buffer_one, logical_end, blast_settings);
    if(temp != NULL) {      //char *name_in = p.name_in; char *taxo_in = p.taxo_in;
      // Requires that the line must consist of enough chars to be consideres a real line:
      if((logical_end-line_start)>50) {
	lines_in_file_found++;
      }
      if(listProteins[my_id]->insert_protein(p.get_taxon_in(), p.get_name_in())) {
	last_block_line = line_start; 	// Sets the end of the last block before a new protein is read. This in order to avoid overlaps in the aprsing of the whole file, and thereby redcusing the mem consumption
      }
      assert(temp < logical_end);
      buffer_one = strchr(temp, '\n');       // update
    } else {
      //! Builds the string to be written, i.e. to improve outpirnt when running mulitple threads.
      char *start = buffer_one;      
      buffer_one = strchr(buffer_one, '\n'); // update
      char *end = logical_end;      
      if(buffer_one) {
	end = buffer_one;
      }
      const int line_length = end - start;
      assert(line_length > 0);
      char *line_block = new char[line_length+1];
      line_block[line_length] = '\0';
      strncpy(line_block, start, line_length);
      //! Print error message and continous
      fprintf(stderr, "!!\t Identified a difficult row in the blastp file:\n%s\n!!\t Suggest you remove it, or fix the problems with the fields of it. To help identify details, we will now include details which should be forwarded to [oekseth@gmail.com]:\n(1)\tId of computer-thread is %d (set number of threads to '1' in order to improve detection of blast errors).\n(2)\tTo estimate the position in the blast-file, we look upon the number of rows found (before this errnous row was dtected) which is %d.\n(3)\tIf this message was not understood, please forward it to the developer at [oekseth@gamil.com], stating your version-number of this software and the code-location of its message, which is [%s]:%s:%d\n", line_block, my_id, (uint)lines_in_file_found, __FUNCTION__, __FILE__, __LINE__);
    }
    if(buffer_one && buffer_one < logical_end) buffer_one += 1;  // , jumping over the line end.
    else buffer_one = logical_end;
    const loint remaining_chars_current = (loint)(logical_end - buffer_one);
    if(remaining_chars_current > remaining_chars) {
      buffer_one = logical_end; 
    }
    else remaining_chars = remaining_chars_current;
    p.free_memory();
  }

  if(last_block_line != NULL) {
    const long int diff = (last_block_line - buffer_start);
    if(diff > 0) return diff-1;
    else return 0;
  } else {return 0;}  
}



/*overrride*/void* pipe_parse_parse::operator()(void* item) {
  int my_id = -1;
  string_section_t *section = NULL;
  loint found_at_index_pos = 0;
  { 
    slock_t::scoped_lock taxalock(taxa_upd);
    for(int i = 0;i<CPU_TOT;i++) {
      if (in_use[i] == false) {in_use[i]=true; my_id=i;i=CPU_TOT;}
    }    
    assert(my_id!= -1);
    
    if(!FIRST_READ) {
      if(stringBuffer) {
	stringBuffer->get_section(section, found_at_index_pos);
	if(!section)  {buffer_string_list::close(stringBuffer);}
      }
    }
  }
  if(FIRST_READ) section = static_cast<string_section_t*>(item);
  else if(!parseBlocks && !section) {
    section = static_cast<string_section_t*>(item);
  }
  parse_send_t     *send_second = NULL; 
  parse_send_first_t *send_first = NULL;
  const clock_t tstart = times(NULL);
  loint lines_in_file_found_for_this_thread = 0;
  if(FIRST_READ) {
    listProteins[my_id] = taxon_list::init();
    uint chars_processed_in_this_run = 0;
    if(section->buffer_main_start && section->buffer_main_logical_end && (section->buffer_main_logical_end - section->buffer_main_start)) {
      chars_processed_in_this_run = (section->buffer_main_logical_end - section->buffer_main_start);
    }
    long long int block_length = parse_blast_blocks_ids(my_id, section[0].buffer_main_start,  section->buffer_main_logical_end, lines_in_file_found_for_this_thread);
#ifndef NDEBUG // If in debug mode, uses an alternative approach
    lines_in_file_found_for_this_thread = debug_get_cnt_newlines(section[0].buffer_main_start,  section->buffer_main_logical_end);    
#endif
    long long int last_line_pos = block_length;
    if(!section->end_of_file) {
      if(last_line_pos < 1) { // In the middle of a data-block
	last_line_pos = 0;
      } else if(!section->end_of_file) {
	last_line_pos+= section->prev_index; // the position of the file
	last_line_pos++; // Updates the position in the file
      } else {
	last_line_pos = section->prev_index + section->getStringLength();
      }
    } else {      
      block_length = section->getStringLength();
    }
    
    send_first = new parse_send_first(listProteins[my_id], /*the pos reading to=*/last_line_pos, section->block_cnt);
    //section->print_data_block(); printf("\n");// Uncomment if data is needed to visually verify.
    send_first->data_resides_in_mem = stringBuffer->update(section,(loint)block_length);   

#ifndef NDEBUG
    { 
      slock_t::scoped_lock taxalock(taxa_upd);
      const loint size_in_buff = stringBuffer->get_length_of_inserted_data(section);
      if(size_in_buff == (chars_processed_in_this_run+1)) {chars_processed_in_this_run++;} 
      total_number_of_chars_processed_in_first_read += chars_processed_in_this_run;// Updates the overview of pairs processed.      
      if(size_in_buff != chars_processed_in_this_run) {
	printf("!!\tsize_in_buff(%llu) != chars_processed_in_this_run(%u) in pipe_parse_parse.cxx at line %d\n",
	       size_in_buff, chars_processed_in_this_run, __LINE__);
	assert(size_in_buff != chars_processed_in_this_run); 
      }
    }
#endif
  } else {
    if(!section && parseBlocks) { // Then the data has not been read into the strcture.
      section = parseBlocks->get_string_section(my_id); 
    } 

    if(section != NULL) {
      // Initializes the specific data:
      // the list containing the  indexes for the protein strings:
      proteinVector[my_id].free_memory();
      proteinVector[my_id] = protein_vector(hashProtein, hashTaxa, taxon_length);   
      // The data for the first protein in this block (used in the next step of this process):
      first_protein[my_id] = protein_relation(); 
      assert(section->is_set());
      const uint tot_elements_start = parseData[my_id]->getTotalLengthOfData();
      uint chars_processed_in_this_run = 0;
      if(section->buffer_main_start && section->buffer_main_logical_end && (section->buffer_main_logical_end - section->buffer_main_start)) {chars_processed_in_this_run = (section->buffer_main_logical_end - section->buffer_main_start);}
      loint this_total_number_of_pairs_overlapping = 0;
      const uint cnt_inserted_pairs = parse_blast_blocks_data(my_id, section->buffer_main_start,  section->buffer_main_logical_end, lines_in_file_found_for_this_thread, this_total_number_of_pairs_overlapping);

#ifndef NDEBUG // If in debug mode, uses an alternative approach
      const loint temp_cnt_newlines_this_block = debug_get_cnt_newlines(section[0].buffer_main_start,  section->buffer_main_logical_end);
      // At some blocks, the newline is missing at the end, therefore a difference of '+1' is correct:
      if((lines_in_file_found_for_this_thread != temp_cnt_newlines_this_block) && (lines_in_file_found_for_this_thread != (temp_cnt_newlines_this_block+1)) ) {
	const int cnt_missing =  temp_cnt_newlines_this_block - lines_in_file_found_for_this_thread;
	printf("!!\t Not all the lines (%d lines is not interpreted) of the blast-file was read, due to errors of the format; please investigage earlier error-messages. If this error is not understood, please contact the developer at [oekseth@gamil.com]. Error at [%s]:%s:%d\n", cnt_missing, __FUNCTION__, __FILE__, __LINE__);
	//	printf("!!\t Prints the block, as lines_in_file_found_for_this_thread(%llu) != temp_cnt_newlines_this_block(%llu) in pipe_parse_parse.cxx at line %d\n", lines_in_file_found_for_this_thread, temp_cnt_newlines_this_block, __LINE__);
	//      section->print_data_block();
	//	assert(!((lines_in_file_found_for_this_thread != temp_cnt_newlines_this_block) && (lines_in_file_found_for_this_thread != (temp_cnt_newlines_this_block+1)) ));
      }
#endif
      send_second = (parse_send_t*)malloc(sizeof(parse_send_t));
      const uint tot_elements = parseData[my_id]->getTotalLengthOfData();
      list_file_parse_t *parse_struct_send = parseData[my_id]->createPacket();
      *send_second = parse_send(parse_struct_send, 
				first_protein[my_id],
				proteinVector[my_id].arrOverlap, max_sim_value[my_id], taxon_length, found_at_index_pos);
      assert((tot_elements - tot_elements_start) == cnt_inserted_pairs);

      if(parse_struct_send) { // The sum must be the same:
	assert(tot_elements == (parseData[my_id]->getTotalLengthOfData() + parse_struct_send->getTotalLengthOfData()));
      } else {assert(tot_elements == parseData[my_id]->getTotalLengthOfData());}
#ifndef NDEBUG
      uint sum_debug_myid = 0;
      for(int taxon =0; taxon< taxon_length; taxon++) {
	for(int prot = 0; prot <hashProtein[taxon].getLength(); prot++) {
	  sum_debug_myid += proteinVector[my_id].arrOverlap[taxon][prot];
	}
      }
      assert(sum_debug_myid == proteinVector[my_id].debug_sum_overlap);
      { 
      	slock_t::scoped_lock taxalock(taxa_upd);
	//! Updates the sum of overlap values:
	debug_sum_inserted_overlap_values += sum_debug_myid;

	total_number_of_pairs_overlapping += this_total_number_of_pairs_overlapping;
	if(chars_processed_in_this_run) chars_processed_in_this_run--; // Removes the last line of the string block.
	total_number_of_chars_processed_in_second_read += chars_processed_in_this_run;// Updates the overview of pairs processed.	
	// Updates the overview of pairs processed:
	if(parse_struct_send) {pairs_in_file_sent_to_merge += parse_struct_send->getTotalLengthOfData();}
	pairs_in_file_found += cnt_inserted_pairs;
      }      
#endif
    } else if(parseBlocks) { 
      if(!parseBlocks->all_data_read(false)) { // The other secions have more data
	send_second = (parse_send_t*)malloc(sizeof(parse_send_t));
	send_second[0] = parse_send();//NULL, first_protein[my_id], NULL, 0);
      }
    }
  }
  bool section_is_set = true;
  if(section) {
    string_section::close(section); // frees the buffer holding the text parsed above
  } else section_is_set = false;

  { // A thread-lock:
    slock_t::scoped_lock taxalock(taxa_upd);
    // Update the lines read:
    if(FIRST_READ) lines_in_file_found_first_read += lines_in_file_found_for_this_thread;
    else lines_in_file_found_last_read += lines_in_file_found_for_this_thread;

    // Deselects the thread id:
    if(!section_is_set && (send_first!= NULL) && (send_second != NULL)) {
#ifndef NDEBUG
      if(!FIRST_READ && (!parseBlocks) && stringBuffer) {
	assert(!stringBuffer->has_data()); // At this point the buffer should be empty
      }
      assert(0==(lines_in_file_found_first_read - lines_in_file_found_for_this_thread));
#endif
      return NULL;
    } else {
      in_use[my_id] = false; // TODO: Verify that this does not cause any thread-collisions.
      logid id = read_first;
      if(!FIRST_READ) id = read_second;
      const clock_t tend = times(NULL); log->append_measurement(id, 1, tstart, tend);
    }
  }
  if(FIRST_READ) return send_first;
  else {
    return send_second;
  }
}

pipe_parse_parse::pipe_parse_parse(loint _reading_file_start_position, loint _reading_file_length, uint _disk_buffer_size, int _taxon_length, const bool _USE_EVERYREL_AS_ARRNORM_BASIS, log_builder_t *_log, char _SEPERATOR, char *_FILE_INPUT_NAME, int _CPU_TOT, bool _USE_LAST_BLAST_CLOMUN_AS_DISTANCE, bool print_norm, bool debug_norm, char *binary_loc, int index_protein, int index_taxon				     
				   , uint _DEFAULT_NUMBER_OF_COLUMNS_IN_NAME, bool _use_improved_overlap_algo
				   ) :
  tbb::filter(  /*seriel=*/ false),   //tbb::filter(  /*seriel=*/ false)
  reading_file_start_position(_reading_file_start_position), reading_file_length(_reading_file_length),
  DEFAULT_NUMBER_OF_COLUMNS_IN_NAME(_DEFAULT_NUMBER_OF_COLUMNS_IN_NAME), use_improved_overlap_algo(_use_improved_overlap_algo), 
  PRINT_NORMALIXATION_BASIS(print_norm), DEBUG_NORM(debug_norm), FILE_BINARY_LOCATION(binary_loc), INDEX_IN_FILE_FOR_PROTEIN_NAME(index_protein),
  INDEX_IN_FILE_FOR_TAXON_NAME(index_taxon), 
  USE_LAST_BLAST_CLOMUN_AS_DISTANCE(_USE_LAST_BLAST_CLOMUN_AS_DISTANCE),
  USE_BEST_BLAST_PAIR_SCORE(false), 
  FILE_INPUT_NAME(_FILE_INPUT_NAME), CPU_TOT(_CPU_TOT), 
  listTaxa(NULL), SEPERATOR(_SEPERATOR), log(_log), 
  send_all_data_to_pipe(true),
  lines_in_file_found_first_read(0), lines_in_file_found_last_read(0),
  pairs_in_file_found(0), pairs_in_file_sent_to_merge(0), total_number_of_pairs_overlapping(0),
  total_number_of_chars_processed_in_first_read(0),
  total_number_of_chars_processed_in_second_read(0),
  FIRST_READ(true),
  taxon_length(_taxon_length),
  //    hash_is_created(_hash_is_created),
  //  n_threads(_n_threads),
  USE_EVERYREL_AS_ARRNORM_BASIS(_USE_EVERYREL_AS_ARRNORM_BASIS), hashTaxa(NULL), parseBlocks(NULL),
  stringBuffer(NULL), debug_sum_inserted_overlap_values(0), debug_sum_inserted_overlap_values_found(0)
{
  in_use = (bool*)malloc(sizeof(bool)*CPU_TOT);
  assert(in_use);
  for(int i = 0;i<CPU_TOT;i++) in_use[i] = false;
  // The below list hold the uniuqe name of the proteins, with
  // distinct list for each thread, as the list is built. The usage is therefore
  // on the parsing' first run. In this contect changing the number of cpu's on the second run is therefore un-problematic. (06.04.2012 by oekseth).
  listProteins = (taxon_list_t**)malloc(sizeof(taxon_list_t*)*CPU_TOT);
  assert(listProteins);
  first_protein = NULL;
  proteinVector = NULL;
  max_sim_value = NULL;
  local_arrNorm = NULL;
  parseData = NULL;
  if(FILE_INPUT_NAME) {
    const loint max_block_cnt = buffer_string_list::get_estimate_of_memory_blocks_cnt(_disk_buffer_size, FILE_INPUT_NAME);
    stringBuffer = new buffer_string_list(max_block_cnt);
  }  
  blast_settings = tsettings_input(INDEX_IN_FILE_FOR_PROTEIN_NAME, INDEX_IN_FILE_FOR_TAXON_NAME, USE_LAST_BLAST_CLOMUN_AS_DISTANCE, SEPERATOR, FILE_INPUT_NAME,DEFAULT_NUMBER_OF_COLUMNS_IN_NAME);
}   

#ifdef assert_code
void pipe_parse_parse::assert_private_parts() {
  // Tests the functions for parsing each row (line) of the balst file:
  bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE = true;
  const bool old_use_last_blast_clomun_as_distance = USE_LAST_BLAST_CLOMUN_AS_DISTANCE;
  const bool old_use_improved_overlap_algo = use_improved_overlap_algo;
  char *str = "name1_taxon1 name2_taxon2\t0\t20\t1\t2\t4\t4\t5\t6\t0.011\t0.013\n"; 
  char *str_start = str, *str_end = str + strlen(str);
  USE_LAST_BLAST_CLOMUN_AS_DISTANCE = false;
  use_improved_overlap_algo = false;
  Parse p; const char SEPERATOR = '_';
  tsettings_input_t blast_settings = tsettings_input(INDEX_IN_FILE_FOR_PROTEIN_NAME, INDEX_IN_FILE_FOR_TAXON_NAME, USE_LAST_BLAST_CLOMUN_AS_DISTANCE, SEPERATOR, FILE_INPUT_NAME,DEFAULT_NUMBER_OF_COLUMNS_IN_NAME);
  assert(blast_extractors::getParseHeaders(p, str_start, str_end, blast_settings));
  assert(0==strcmp("name1", p.getInnerName()));
  assert(0==strcmp("name2", p.getOuterName()));
  assert(0==strcmp("taxon1", p.getTaxonInName()));
  assert(0==strcmp("taxon2", p.getTaxonOutName()));

  str_start = getOverlapAndUpdate(p, str_start);
  if(!use_improved_overlap_algo) assert(p.get_overlap_in() == 20);
  float sim_val = 0;  
  str_start = getDistanceUpdate(p, str_start, sim_val, str_end, USE_LAST_BLAST_CLOMUN_AS_DISTANCE);
  //  p.print();
  assert((int)p.get_distance() == (int)(-1*log10(0.011)));

  {
    use_improved_overlap_algo = true;
    char *str_start = str, *str_end = str + strlen(str);    
    Parse p; 
    assert(blast_extractors::getParseHeaders(p, str_start, str_end, blast_settings));
    assert(0==strcmp("name1", p.getInnerName()));
    assert(0==strcmp("name2", p.getOuterName()));
    assert(0==strcmp("taxon1", p.getTaxonInName()));
    assert(0==strcmp("taxon2", p.getTaxonOutName()));

    str_start = getOverlapAndUpdate(p, str_start);
    if(!use_improved_overlap_algo) assert(p.get_overlap_in() == 20);
    float sim_val = 0;  
    str_start = getDistanceUpdate(p, str_start, sim_val, str_end, USE_LAST_BLAST_CLOMUN_AS_DISTANCE);
    //    p.print();
    assert((int)p.get_distance() == (int)(-1*log10(0.011)));

    USE_LAST_BLAST_CLOMUN_AS_DISTANCE = true;
  }
  USE_LAST_BLAST_CLOMUN_AS_DISTANCE = old_use_last_blast_clomun_as_distance;
  use_improved_overlap_algo = old_use_improved_overlap_algo;
}
#endif
void pipe_parse_parse::assert_class(const bool print_info) {
  const static char *class_name = "pipe_parse_parse";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code 
  uint taxon_length = 2; // n_threads = 2;
  log_builder_t *log = log_builder::init();
  const char SEPERATOR = '_';
  const int CPU_TOT = 2; const bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE = true;
  bool print_norm = false; bool debug_norm = false; char *binary_loc = NULL; int index_protein = 0; int index_taxon = 1;
  
  uint DEFAULT_NUMBER_OF_COLUMNS_IN_NAME = 3; bool use_improved_overlap_algo = true;
  //  const uint file_size = (uint)log_builder::get_file_size(FILE_INPUT_NAME, __LINE__, __FILE__);
  pipe_parse_parse temp(0, 0, 1024*1024, taxon_length, false, log, SEPERATOR, "all.blast", CPU_TOT, USE_LAST_BLAST_CLOMUN_AS_DISTANCE,
			print_norm, debug_norm, binary_loc, index_protein, index_taxon, DEFAULT_NUMBER_OF_COLUMNS_IN_NAME, use_improved_overlap_algo
			);
  temp.assert_private_parts();
  list_file_parse_t *temp_s = NULL;
  temp.free_memory(temp_s);

  log_builder::close(log);
  //  if(log) {log->free_mem(), delete log, log = NULL;}
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

