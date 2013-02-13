#include "blast_filtering.h" 

void blast_filtering::exec_filtering(int n_thread, log_builder_t *log, const bool is_inpa) {
#ifdef USE_MPI
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
#endif

  uint list_length = 0;
  bool is_ortholog = false;
  pipe_t TYPE_OF_OPERATION = INPA;
  if(!is_inpa) {
    is_ortholog = true;
    TYPE_OF_OPERATION = ORTH;
  }
  if(is_ortholog) {
    listOrtho = id_simil_list(listTaxa[taxon_length-1].rel_end, listTaxa);
#ifndef NDEBUG
    
#endif
  }
  const bool is_co_orth = false;
  uint *array = listParseData->get_thread_allocations_for_taxa(list_length, is_ortholog, is_inpa, is_co_orth);
  uint current_position = 0;
  uint taxon_start = 0, taxon_end = 0;
  while(listParseData->get_interval(list_length, array, current_position, taxon_start, taxon_end, is_ortholog, is_inpa, is_co_orth)) {
    if(taxon_end != 0) {
      pipe_bucket buck(is_inpa, taxon_start, taxon_end, taxon_length, log, listParseData, listStructData, listTaxa, n_thread);
      pipe_ortholog_inparalog parse(taxon_length, is_inpa, USE_EVERYREL_AS_ARRNORM_BASIS, n_thread, log,
				    listOrtho, listParseData, AMINO_LIMIT, max_input_value, MIN_SIMILARITY_LIMIT,
				    use_improved_overlap_algo, DEBUG_PRINT_DISCARDED_PAIRS, PRINT_OVERLAP_VALUES_ABOVE,
				    PRINT_NORMALIXATION_BASIS, DEBUG_NORM, listTaxa, MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, FILE_BINARY_LOCATION
				    );
      pipe_merge merge(taxon_length, USE_EVERYREL_AS_ARRNORM_BASIS, TYPE_OF_OPERATION, log, listStructData, listTaxa); 
      if(!USE_EVERYREL_AS_ARRNORM_BASIS)  merge.set_arrNorm(arrNorm);
      pipeline pipe;
      pipe.add_filter(buck);
      pipe.add_filter(parse);
      pipe.add_filter(merge);
      pipe.run(n_thread);  
      pipe.clear();
      if(!USE_EVERYREL_AS_ARRNORM_BASIS)  arrNorm = merge.get_arrNorm();
      buck.free_data();
      parse.free_additional_blocks();
      listStructData = merge.getFileStruct(n_thread, MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, FILE_BINARY_LOCATION);

    }
    /*
    if(is_ortholog && taxon_end == (uint)taxon_length) {
      listOrtho.set_DEBUG_PRINT_DISCARDED_PAIRS(DEBUG_PRINT_DISCARDED_PAIRS, taxon_length, listTaxa);
      listOrtho.execInternGlobalReciProc(listTaxa[taxon_length-1].rel_end);
    }
    */
  }
}

void blast_filtering::exec_co_orth(int n_thread, log_builder_t *log, stack_rel *&stackRel) {
  //! Some initial escpectations:
  assert(log);
  assert(taxon_length);
  assert(stackRel);
  assert(listStructData);
  assert(listParseData);

  //! Starts the "pipelining"
  const bool is_ortholog = false;
  const bool is_inpa = false;
  const bool is_co_orth = true;
  uint list_length = 0;  
  uint *array = listStructData->get_thread_allocations_for_taxa(list_length, is_ortholog, is_inpa, is_co_orth);
  uint current_position = 0;
  uint taxon_start = 0, taxon_end = 0;
  while(listParseData->get_interval(list_length, array, current_position, taxon_start, taxon_end, is_ortholog, is_inpa, is_co_orth)) {
    if(taxon_end != 0) {
      pipe_norm merge(taxon_length, USE_EVERYREL_AS_ARRNORM_BASIS, INPA_INPA, log);  
      merge.arrNorm = arrNorm;
      pipeline pipe; 
      float **temp = NULL;
      if(RESTRICTED_DEFENITION) {
	// listTaxa[taxon_length-1].rel_end,
	pipe_struct parse(n_thread,  taxon_start, taxon_end, INPA_ORTH, USE_EVERYREL_AS_ARRNORM_BASIS, NULL, log,
			  MODE_PAIRWISE_OUTPUT, MODE_INTEGER_OUTPUT, DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT,
			  DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT, listOrtho, listTaxa, stackRel, listStructData, NULL, max_input_value
			  , MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT,
			  PRINT_IN_MCL_FORMAT, SORT_ABC_DATA, PRINT_NORMALIXATION_BASIS, DEBUG_NORM, temp
			  );
	parse.set_DEBUG_PRINT_DISCARDED_PAIRS(DEBUG_PRINT_DISCARDED_PAIRS);
	pipe.add_filter(parse);
	pipe.add_filter(merge);
	pipe.run(n_thread);
	pipe.clear(); 
	parse.finalize_memory(taxon_length); 
      } else {
	// listTaxa[taxon_length-1].rel_end,
	pipe_struct parse(n_thread,  taxon_start, taxon_end, INPA_INPA, USE_EVERYREL_AS_ARRNORM_BASIS, NULL, log,
			  MODE_PAIRWISE_OUTPUT, MODE_INTEGER_OUTPUT, DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT,
			  DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT, listOrtho, listTaxa, stackRel, listStructData, NULL/*listParseData*/, max_input_value
			  , MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT,
			  PRINT_IN_MCL_FORMAT, SORT_ABC_DATA, PRINT_NORMALIXATION_BASIS, DEBUG_NORM, temp
			  );
	parse.set_DEBUG_PRINT_DISCARDED_PAIRS(DEBUG_PRINT_DISCARDED_PAIRS);
	pipe.add_filter(parse);
	pipe.add_filter(merge);
	pipe.run(n_thread);
	pipe.clear(); 
	parse.finalize_memory(taxon_length); 	
      }
      arrNorm = merge.arrNorm;	
    }
  }
}

void blast_filtering::write_result_file(int n_thread, log_builder_t *log, stack_rel *&stackRel, pipe_write &write, float **&arr_avgNorm) {
  // Only gets the data for the inparalogs:
  uint list_length = 0;
  assert(listStructData);
  uint *array = listStructData->get_thread_allocations_for_taxa(list_length, false, true, false);
  uint current_position = 0;
  uint taxon_start = 0, taxon_end = 0;
  while(listStructData->get_interval(list_length, array, current_position, taxon_start, taxon_end, false, true, false)) {
    if(taxon_end != 0) {
      pipe_struct parse(n_thread,  taxon_start, taxon_end, DUMP, USE_EVERYREL_AS_ARRNORM_BASIS, arrNorm, log,MODE_PAIRWISE_OUTPUT, MODE_INTEGER_OUTPUT, DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT,
			DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT, listOrtho, listTaxa, stackRel, listStructData, NULL/*listParseData*/, max_input_value
			, MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT,
			PRINT_IN_MCL_FORMAT, SORT_ABC_DATA, PRINT_NORMALIXATION_BASIS, DEBUG_NORM, arr_avgNorm
			); 
      parse.set_DEBUG_PRINT_DISCARDED_PAIRS(DEBUG_PRINT_DISCARDED_PAIRS);
      pipeline pipe; 
      pipe.add_filter(parse);
      pipe.add_filter(write);
      pipe.run(n_thread);
      pipe.clear();	  
      parse.finalize_memory(taxon_length); 
    }
  
  }
}
void blast_filtering::print_class_info() {
  fprintf(stdout, "\n---------------------------------\n");
  fprintf(stdout, "The purpose of this wrapper module is handling a blast file:\n" \
	  " \t(a) Handles user parameters for reading the blast file.\n" \
	  " \t(b)   The purpose of this wrapper module is handling a the filtering process given an input of a taxa_t- and list_file_parse containers, producing a filtered output.\n");
  fprintf(stdout, "\n---------------------------------\n");
}


/** 
    @brief Updates the command line interface with values to be set by the user:
    @param <cmd> The input list to be made available from the terminal console.
    @parm <first_pass> In order to get the correct order of the fields, i.e. "OUTPUT" before "OPERATIONAL".
**/
void blast_filtering::init_values(cmd_list *cmd, const bool first_pass) {
  class cmd_argument cl2; 
  if(first_pass) {
    cl2 = cmd_argument("Send the output in the form of MCL native matrix (all.mcl) to STDOUT for piping", "P", "pipe", BOOLEAN, &OUTPUT_PIPE_MCI_ALL, "OUTPUT"); 
    cmd->add_cmd_argument(cl2); // OUTPUT_PIPE_MCI_ALL

    cl2= cmd_argument("Sort *.abc files by the scores", "S", "sort", BOOLEAN, &SORT_ABC_DATA, "OUTPUT");
    cmd->add_cmd_argument(cl2);  // SORT_ABC_DATA
  } else {
    cl2 = cmd_argument("Skip normalization of similarity scores", "w", "without_norm", BOOLEAN, &DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT, "OPERATIONAL");
    cmd->add_cmd_argument(cl2);  // DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT

    cl2 = cmd_argument("The e-value cutoff (positive floating number, see examples below), mutually exclusive with -u", "e", "cutoff", FLOAT, &MIN_SIMILARITY_LIMIT, "FILTERING");
    //    cl2 = cmd_argument("Threshold for protein pair similarity (positive floating number, e.g. set to '7' or '7.0' to exclude protein pairs with e-values above '1e-07')", "e", "threshold", FLOAT, &MIN_SIMILARITY_LIMIT, "FILTERING");
    cmd->add_cmd_argument(cl2); // MIN_SIMILARITY_LIMIT      

    cl2 = cmd_argument("Threshold for protein pair overlap (integer [1-100], see examples below)", "o", "overlap", UINT_NOT_NULL, &AMINO_LIMIT, "FILTERING");
    cmd->add_cmd_argument(cl2); // AMINO_LIMIT

    cl2 = cmd_argument("Restricted definition of co-orthologs, relations between inparalogs of orthologs are excluded. Useful to give more weight to orthologs during MCL computation", "C", "strict_coorthologs", BOOLEAN, &RESTRICTED_DEFENITION, "FILTERING");
    cmd->add_cmd_argument(cl2); // AMINO_LIMIT


#ifndef NDEBUG
    cl2 = cmd_argument("Print protein pairs discarded during filtering", "d", "discarded_pairs", BOOLEAN, &DEBUG_PRINT_DISCARDED_PAIRS, "DEBUG");
    cmd->add_cmd_argument(cl2); // 

    // cl2 = cmd_argument("Print meta-information about the parsing- and filtering. Useful when analysing the result", "pbfd", "print_blast_filter_data", BOOLEAN, &blastInfo.build_meta_blast, "OUTPUT");
    // cmd->add_cmd_argument(cl2); // 
#endif
  }
}

//! Initiates the list for parsing input arguments fromthe terminal:
cmd_list* blast_filtering::init_cmd_list(char *DEFAULT_OPTION_NAME, uint &DEFAULT_OPTION_NAME_COUNT, char *FILE_INPUT_NAME) {
  cmd_list *list = new cmd_list(DEFAULT_OPTION_NAME, DEFAULT_OPTION_NAME_COUNT, FILE_INPUT_NAME); 
  return list;
}


void blast_filtering::build_cmd_list(cmd_list *cmd, int argc, char *argv[]) { 
  bool DEBUG_BUILD_MAN_PAGE = false;
#ifdef BUILD_MAN_PAGE
  DEBUG_BUILD_MAN_PAGE = true;
#else 
  DEBUG_BUILD_MAN_PAGE = false;
#endif
  //cmd->print_arguments(stdout, NULL);
  cmd->set_arguments(argv, argc, DEBUG_BUILD_MAN_PAGE);  
  cmd->clean_cmd_argument(); // removes the data about the arrguments added. (Clens memory).  
  cmd->free_mem();
}

//! @return the transformed value.
float transform_threshold_value_to_e_if_last_column_set(bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE, float MIN_SIMILARITY_LIMIT) {
  if(USE_LAST_BLAST_CLOMUN_AS_DISTANCE) {
    //    const float old_threshold = MIN_SIMILARITY_LIMIT;
    if(MIN_SIMILARITY_LIMIT != 0) {
      fprintf(stderr, "!!\t orthAgogue does not allow a threshold value when the last column is used. The threshold value therefore changes from %f into %f. If questions, please read the manual or alternatively contact the developer at [oekseth@gmail.com].\n", MIN_SIMILARITY_LIMIT, 0.0);
      return 0;
//     }
//     if(true) {      
//       const float new_threshold = pow(10, MIN_SIMILARITY_LIMIT);
// #ifndef NDEBUG
//       const float back_subst = -1*log10(new_threshold);
//       assert(back_subst >= 0);
//       assert((uint)back_subst == (int)old_threshold);
// #endif
//       return new_threshold;
//     } else {
//       if(MIN_SIMILARITY_LIMIT) {
// 	const float new_threshold = -1*log10(MIN_SIMILARITY_LIMIT);
//       } else return 0;
    }
  } 
  return MIN_SIMILARITY_LIMIT; // did not change
}
/**
   @brief Sets the user settings from the blast parsing process.
   @param <obj> The object of class tsettings_input containing the data to use.
   @remarks The filtering process will depend upon some of these settings, therby relevant.
**/
void blast_filtering::set_values(tsettings_input obj) {
  DEBUG_NORM = obj.DEBUG_NORM; 
  PRINT_NORMALIXATION_BASIS = obj.PRINT_NORMALIXATION_BASIS; 
  USE_EVERYREL_AS_ARRNORM_BASIS = obj.USE_EVERYREL_AS_ARRNORM_BASIS;
  FILE_BINARY_LOCATION = obj.FILE_BINARY_LOCATION; 
  //  SEPERATOR = _SEPERATOR; 
  CPU_TOT = obj.CPU_TOT;
  
  //! Updates the similarity limit:
  MIN_SIMILARITY_LIMIT = transform_threshold_value_to_e_if_last_column_set(obj.USE_LAST_BLAST_CLOMUN_AS_DISTANCE, MIN_SIMILARITY_LIMIT);
}
//! Initializes values after those set in the blast parsing (dound in library blast_parsing):
void blast_filtering::set_values(bool _DEBUG_NORM, bool _PRINT_NORMALIXATION_BASIS, bool _USE_EVERYREL_AS_ARRNORM_BASIS, char *_FILE_BINARY_LOCATION/*,  char _SEPERATOR*/, int _CPU_TOT, bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE) {
  DEBUG_NORM = _DEBUG_NORM; 
  PRINT_NORMALIXATION_BASIS = _PRINT_NORMALIXATION_BASIS; 
  USE_EVERYREL_AS_ARRNORM_BASIS = _USE_EVERYREL_AS_ARRNORM_BASIS;
  FILE_BINARY_LOCATION = _FILE_BINARY_LOCATION; 
  //  SEPERATOR = _SEPERATOR; 
  CPU_TOT = _CPU_TOT;
  //! Updates the similarity limit:
  MIN_SIMILARITY_LIMIT = transform_threshold_value_to_e_if_last_column_set(USE_LAST_BLAST_CLOMUN_AS_DISTANCE, MIN_SIMILARITY_LIMIT);
}

/**
   @brief Executes the main operation for the filtering:
   @param <log> The log object to store measurements in.
   @param <bp>  The object containing the basis whom the filtering will work on.
**/
void blast_filtering::start_filtering(log_builder_t *log, bp_container_t &bp) {
#ifdef USE_MPI
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
#endif
  FILE *f_log = log_builder::get_log_file_pointer(CPU_TOT, "blast_filtering", __FILE__, __LINE__, true); // The log file.
  log->start_measurement(complete_blast_filtering);

  if(bp.has_data()) {
    //Sets the data:
    bp.get_variables(arrNorm, listTaxa, taxon_length, listParseData, max_input_value);
    if(blastInfo.build_meta_blast) {
      //! Adds meta-data, i.e. before the filtering starts
      blastInfo.FILE_INPUT_NAME = bp.FILE_INPUT_NAME;
      blastInfo.taxon_length = taxon_length;
      blastInfo.cnt_unique_proteins = listTaxa[taxon_length-1].getRelativeEndIndex();
      blastInfo.cnt_relations_orthologs_Before_filtering = listParseData->getTotalLengthOfData_for_orthologs();
      blastInfo.cnt_relations_inparalogs_Before_filtering = listParseData->getTotalLengthOfData_for_inparalogs();
      blastInfo.cnt_relations_Total_Before_filtering = listParseData->getTotalLengthOfData(); 
    }
    // Some preprosessing of values et:
    if(DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT) { // Updates both normalization-variables:
      DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT = true; // If set to false, prints the similarity score for the proteins during the dump to mcl for the mcl format
    } else {
      DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT = false; // If set to false, prints the similarity score for the proteins during the dump to mcl for the mcl format
    }
    if(AMINO_LIMIT >= 100) {
      log_builder::throw_warning(input_terminal,  __LINE__, __FILE__, __FUNCTION__, "Discoved the cut-off value set to 100 or above. Continues exectution, although the results will of less worth");
    }
    mcl_t TYPE_OF_RESULTFILE_TO_STDOUT = none;
    if(OUTPUT_PIPE_MCI_ALL) { TYPE_OF_RESULTFILE_TO_STDOUT = all;} // Prints the mci out result file to stdout.
    log->start_measurement(filter_orthologs);
    uint n_thread = CPU_TOT;     
    task_scheduler_init init(n_thread);
    exec_filtering(n_thread, log, /*is_inpa=*/false);
#ifdef USE_MPI
    if(!arrNorm) arrNorm = new list_norm(taxon_length, max_input_value, PRINT_NORMALIXATION_BASIS, DEBUG_NORM);
    //! Sends- and receives the list of arrInpaLimit* objects accross nodes.
    taxa::mpi_send_and_receive_arrInpaLimit(listTaxa, taxon_length);
#else
    assert(arrNorm != NULL);
#endif
    log->end_measurement(filter_orthologs);
    log->start_measurement(filter_inparalogs);
    if(listStructData) { //---------------- The Inparalog operation:
#ifdef USE_MPI
      if(listParseData && listStructData) {
	//! Updates the list of taxon responsibilities:
	listStructData->copy_list_of_nodes_taxa_responsibilities(listParseData);
      }
#endif
      if(f_log) fprintf(f_log, "before_inpa(%u)", listStructData->getTotalLengthOfData());

      exec_filtering(n_thread, log, /*is_inpa=*/true);

      listOrtho.set_DEBUG_PRINT_DISCARDED_PAIRS(DEBUG_PRINT_DISCARDED_PAIRS, taxon_length, listTaxa);
#ifdef USE_MPI
    int number_of_nodes = 1; MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes);
      if(number_of_nodes > 1) {
	//! A call sending- and receiving the ortholog pairs:
	listOrtho.mpi_tx_rx_recip_tx_rx(taxon_length);
	//! Sends- and receives the list of rel_t* objects accross nodes.
	assert(listStructData);
	listStructData->mpi_make_data_consistent_accross_nodes(taxon_length);	
      } 
#endif
      //#ifndef NDEBUG
      //! Write the data to a file before performing the reciprocal ortholog-pair-matching:
      //! Note: Motiavted by the need to exemplify why higher overlap-threshold may result in more ortholog-pairs:
      if(false) { 
	FILE *f_log = log_builder::get_log_file_pointer("orthologs_not_recip", CPU_TOT, __FILE__, __LINE__);
	assert(f_log);
	listOrtho.printGlobalArr(f_log, true, listTaxa,taxon_length);
	fclose(f_log);
      }
      //#endif
      //      const uint relations_before_recip = listOrtho.get_total_number_of_pairs();
      if(blastInfo.build_meta_blast) {
	//! Adds meta-data, i.e. after the filtering is performed:
	//blastInfo.cnt_relations_orthologs_After_filtering = listStructData->getTotalLengthOfData_for_orthologs();
	blastInfo.cnt_relations_orthologs_Before_filtering = listOrtho.get_total_number_of_pairs();	
      }
      listOrtho.execInternGlobalReciProc(listTaxa[taxon_length-1].rel_end);
      //#ifndef NDEBUG
      //! Write the data to a file before performing the reciprocal ortholog-pair-matching:
      //! Note: Motiavted by the need to exemplify why higher overlap-threshold may result in more ortholog-pairs:
      if(false) { 
	f_log = log_builder::get_log_file_pointer("orthologs_reciprocal", CPU_TOT, __FILE__, __LINE__);
	assert(f_log);
	listOrtho.printGlobalArr(f_log, true, listTaxa,taxon_length);
	fclose(f_log);
      }
      //#endif
      if(f_log) fprintf(f_log, " after_inpa(%u)", listStructData->getTotalLengthOfData());
      log->end_measurement(filter_inparalogs);

      assert(listParseData);
      assert(listStructData);      
      log->start_measurement(filter_co_orthologs);
      {

	//--------------- The inpa_inpa_operation 
	assert(listStructData);
	log_builder::test_memory_condition_and_if_not_abort(listStructData!=NULL, __LINE__, __FILE__, __FUNCTION__);
	stack_rel *stackRel;// Holds the co orthologs for every protein in a stack; intialized in 'ortho_set' 
	stackRel = new stack_rel[listTaxa[taxon_length-1].rel_end]; 		
	exec_co_orth(n_thread, log, stackRel);

#ifndef NDEBUG
	loint rel_cnt = 0; for(int i =0; i< listTaxa[taxon_length-1].rel_end;i++) rel_cnt+=stackRel[i].unsafe_size();
	if(f_log) fprintf(f_log, " co-orthologs(%lld)", rel_cnt);
#endif

	//
	//------------- The File Writing Operation
	log->end_measurement(filter_co_orthologs);
	{
	  log->start_measurement(write_resultfile);
#ifdef USE_MPI
	  if(number_of_nodes>1) {
	    //! Sends- and receives the norm basis:
	    if(!USE_EVERYREL_AS_ARRNORM_BASIS) {
	      if(!arrNorm) {
		arrNorm = new list_norm(taxon_length, max_input_value, PRINT_NORMALIXATION_BASIS, DEBUG_NORM);
		arrNorm->set_DEBUG_NORM(DEBUG_NORM);
	      }
	      arrNorm->mpi_make_data_consistent_accross_nodes(MPI_COMM_WORLD, myrank, number_of_nodes);
	    }
	    
	    // Sends- and receives the co-orthologs (stackRel)
	    id_simil_list::mpi_send_co_orthologs_accross_nodes(listTaxa, taxon_length, stackRel, listStructData);
	  }
#endif

#ifndef NDEBUG
	  uint cnt_co_orthologs = 0;
	  for(int i =0; i< listTaxa[taxon_length-1].rel_end;i++) cnt_co_orthologs+= (uint)stackRel[i].unsafe_size();
#endif

	  if(listParseData) {
	    listParseData->free_memory(true);
	    delete listParseData; listParseData = NULL;
	  }

	  if(blastInfo.build_meta_blast) {
	    //! Adds meta-data, i.e. after the filtering is performed:
	    //blastInfo.cnt_relations_orthologs_After_filtering = listStructData->getTotalLengthOfData_for_orthologs();
	    blastInfo.cnt_relations_orthologs_After_filtering = listOrtho.get_total_number_of_pairs();
	    blastInfo.cnt_relations_inparalogs_After_filtering = listStructData->getTotalLengthOfData_for_inparalogs();
	    blastInfo.cnt_relations_Total_After_filtering = listStructData->getTotalLengthOfData(); 
	  }

	  // The file_write_operation 
	  float **arrAvgNorm = NULL;
	  if(arrNorm) arrAvgNorm = arrNorm->build_basis();
	  else       log_builder::throw_warning(pointer,  __LINE__, __FILE__, __FUNCTION__, "arrNorm not defined");
	  pipe_write write(log, listTaxa, taxon_length, FILE_BINARY_LOCATION, MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT, TYPE_OF_RESULTFILE_TO_STDOUT, SORT_ABC_DATA);
	  write_result_file(n_thread, log, stackRel, write, arrAvgNorm);
	  write.free_mem(SORT_ABC_DATA, FILE_BINARY_LOCATION, CPU_TOT); 
	  log->end_measurement(write_resultfile);
	  // Deallocates the arrAvgNorm;
	  if(arrAvgNorm) {
	    for(uint i = 0; i< (uint)taxon_length; i++) free(arrAvgNorm[i]);
	    free(arrAvgNorm); arrAvgNorm = NULL;
	  }

	  blastInfo.print_meta_data(stdout);
	}  
	if(stackRel){delete [] stackRel, stackRel = NULL;}
      }
    } else fprintf(stderr, "!!\tThe object holding all the not-filtered-away pairs (from the ortholog operation) is empty, implying that no data remains after the filtering process. It this is not what you expected, contact the developer at oekseth@gmail.com, inkluding the following information: Message was generated at line %d in file %s, located in method %s.\n", __LINE__, __FILE__, __FUNCTION__);    
    listOrtho.freeGlobalArr();  

    list_file_parse<rel>::close(listStructData, true);
  } else {
    fprintf(stderr, "!!\tData not set for call on method %s in class blast_filtering at line %d. Probaly wrong usage. Either read the documentation or contact oekseth@gmail.com. Software now aborts.\n", __FUNCTION__, __LINE__);
    exit(2);
  }
  log->end_measurement(complete_blast_filtering);
  if(f_log) {
    fprintf(f_log, "\n");
    fclose(f_log);
    f_log = NULL;
  }
  if(arrNorm) list_norm::close(arrNorm);
}
//! De-allocates the data:
void blast_filtering::free_memory() {
  if(listTaxa) {
    for(int i =0; i< taxon_length; i++) listTaxa[i].free_mem();
    delete [] listTaxa, listTaxa = NULL;
    }
} 
void blast_filtering::close(blast_filtering *&obj) {
  if(obj) {obj->free_memory(), delete obj, obj = NULL;}
}

blast_filtering::blast_filtering(cmd_list *cmd) :
  // Variables given as input:
  OUTPUT_PIPE_MCI_ALL(false), DEBUG_NORM(false),  PRINT_NORMALIXATION_BASIS(false),
  USE_EVERYREL_AS_ARRNORM_BASIS(false),  FILE_BINARY_LOCATION("")/*, SEPERATOR('_')*/, CPU_TOT(1),
  arrNorm(NULL), listTaxa(NULL), taxon_length(0), listParseData(NULL), max_input_value(0),
  // Variables of this scope:
  DEBUG_PRINT_DISCARDED_PAIRS(false), DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT(true),
  DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT(true), MODE_PAIRWISE_OUTPUT(false), MODE_INTEGER_OUTPUT(false),
  PRINT_IN_ABC_FORMAT(true), PRINT_IN_MCL_FORMAT(true), PRINT_OVERLAP_VALUES_ABOVE(false), 
  RESTRICTED_DEFENITION(false), SORT_ABC_DATA(false), AMINO_LIMIT(0), MIN_SIMILARITY_LIMIT(0), listStructData(NULL)
{
  init_values(cmd);
}

void blast_filtering::assert_class(const bool print_info) {
  const static char *class_name = "blast_filtering";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  //#include "pipe_struct.h"
#include "build_string.h"
#include "pipe_ortholog_inparalog.h"
#include "pipe_struct.h"
#include "pipe_bucket.h"
#include "pipe_merge.h"
#include "pipe_norm.h"
#include "pipe_write.h"
  pipe_ortholog_inparalog::assert_class(print_info);
  pipe_struct::assert_class(print_info);
  pipe_bucket::assert_class(print_info);
  pipe_merge::assert_class(print_info);
  pipe_norm::assert_class(print_info);
  pipe_write::assert_class(print_info);
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

