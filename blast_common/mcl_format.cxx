#include "mcl_format.h"
#include <math.h>

//! Updates the result cheme for use in the log:
void mcl_format::log_update_result_cheme(pipe_struct_result &result) {
  for(uint i =0; i<mcl_t_size; i++) {
    const loint length = string[i]->get_current_size_string();
    result.append_length(i, length);
#ifndef NDEBUG
    //! Tests if it forresponds to our calculated insertion:
    if(file_chunk) {
      if(file_chunk->get_length(i) != (uint)length) {
	fprintf(stderr, "%s\tfirst_element_size(%u) != other_element_size(%u), at %s:%d\n", get_file_name(i), file_chunk->get_length(i), (uint)length, __FILE__, __LINE__);
	assert(file_chunk->get_length(i) == (uint)length);
      }
    } 
#endif
  }
}

//! Sets the values for the header using properties of this object.
void mcl_format::set_debug_header_size(uint world_index, char *name, uint type_string, uint type_number) {
#ifndef NDEBUG
  assert(file_chunk);
  if(file_chunk) {
    //! Assumes it is called before the insertion, i.e. we validate the length before our update of the lengths:
    const uint length_string = string[type_string]->get_current_size_string();
    assert(length_string == file_chunk->get_size(type_string));
    const uint length_number = string[type_number]->get_current_size_string();
    assert(length_number == file_chunk->get_size(type_number));
    //! Performs the insertion
    file_chunk->insert_sizes_header(world_index, name, type_string, type_number);
  }
#endif
}


//! Sets the values for the header using properties of this object.
void mcl_format::set_debug_pair_size(uint world_index_out, char *label_in, char *label_out, float sim_score, uint type_string, uint type_number) {
#ifndef NDEBUG
  //! Assumes it is called before the insertion, i.e. we validate the length before our update of the lengths:
  assert(file_chunk);
  if(file_chunk) {
    //! A verification procedure:
    const uint length_string = string[type_string]->get_current_size_string();
    assert(length_string == file_chunk->get_size(type_string));
    const uint length_number = string[type_number]->get_current_size_string();
    assert(length_number == file_chunk->get_size(type_number));

    //! The actual isnertion:
    file_chunk->insert_sizes_for_a_pair(world_index_out, label_in, label_out, sim_score, type_string, type_number);
  }
#endif    
}

//! Sets the values for the header using properties of this object.
void mcl_format::set_debug_pair_size(uint world_index_out, char *label_out, float sim_score, uint type_string, uint type_number) {
#ifndef NDEBUG
  //! Assumes it is called before the insertion, i.e. we validate the length before our update of the lengths:
  assert(file_chunk);
  if(file_chunk) {
    //! A verification procedure:
    const uint length_string = string[type_string]->get_current_size_string();
    assert(length_string == file_chunk->get_size(type_string));
    const uint length_number = string[type_number]->get_current_size_string();
    assert(length_number == file_chunk->get_size(type_number));

    //! The actual isnertion:
    file_chunk->insert_sizes_for_a_pair(world_index_out, label_out, sim_score, type_string, type_number);
  }
#endif    
}

  
//! Initiatlizes the mcl-bunch-classes:
void mcl_format::init_string() {
  string = new mcl_bunch_t*[mcl_t_size];
  log_builder::test_memory_condition_and_if_not_abort(string, __LINE__, __FILE__, __FUNCTION__);    
  if(string) {
    for(uint i =0; i<mcl_t_size; i++) {
      string[i] = mcl_bunch::init_class();
    }
  }
}

/**! Prints the current status
   @Changed: 04.03.2011 by oekseth
*/
void mcl_format::print_size_of_data() {
  assert(string);
  assert(string[inpa_number]);
  assert(string[pair_orth_number]);
  assert(string[orth_inpa_number]);
  assert(string[all_number]);
  int myrank = 0; 
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
#endif
  loint cnt_orth = 0;    if(string[pair_orth]) cnt_orth = string[pair_orth]->get_current_size_string();
  loint cnt_inpa = 0;    if(string[inpa]) cnt_inpa = string[inpa]->get_current_size_string();
  loint cnt_co_orth = 0; if(string[orth_inpa]) cnt_co_orth = string[orth_inpa]->get_current_size_string();
  loint cnt_all = 0;     if(string[all]) cnt_all = string[all]->get_current_size_string();
    
  printf("[myrank=%d]\t The mcl_format-object has orthologs(%lld), inparalogs(%lld), co-orthologs(%lld), all(%lld), in mcl_format.cxx:%d\n", 	 myrank, cnt_orth, cnt_inpa, cnt_co_orth, cnt_all, __LINE__);
}


//! Frees the memory reserved for the strings
void mcl_format::free_string() {
  if(string != NULL) {
    for(uint i =0; i<mcl_t_size; i++) {
      mcl_bunch::free_class(string[i]);
    }
    delete[] string;  string = NULL;
  }
#ifndef NDEBUG
  if(file_chunk) {file_chunk->free_memory(); delete file_chunk; file_chunk=NULL;}
#endif
}


/**! Inserts the headers for inparalogs.
   @Changed: 07.01.2011 by oekseth
*/
void mcl_format::set_header_inpa(uint world_index, char *name) {  
  assert(string);
  assert(string[inpa]);
  assert(string[inpa_number]);
  //! Sets the header size:
  set_debug_header_size(world_index, name, inpa, inpa_number);
  //! Sets the actual header:
  if(PRINT_IN_ABC_FORMAT && !MODE_PAIRWISE_OUTPUT_ABC) string[inpa]->set_header(name);
  if(PRINT_IN_MCL_FORMAT && !MODE_PAIRWISE_OUTPUT_MCL) string[inpa_number]->set_header(world_index);
}

/**
   @Name: insert_inpa(..) -- Inserts inparalogs: Does not check the uniqueness. (The data source to be used has unique values)
   @Changed: 07.01.2011 by oekseth
*/
void mcl_format::insert_inpa(uint world_index, char *name, float sim_score, float div_factor) {
  //! If in debug mode, sums the sizes of the estimated file size:
  set_debug_pair_size(world_index, name, sim_score, inpa, inpa_number);
  //! Inserts the actual data:
  assert(string);
  assert(string[inpa]);
  assert(string[inpa_number]);
  assert(MODE_PAIRWISE_OUTPUT_MCL==false && MODE_PAIRWISE_OUTPUT_ABC == false);
  if(PRINT_IN_ABC_FORMAT) {
    if(DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT) 
      string[inpa]->insert(name, sim_score, div_factor);
    else string[inpa]->insert(name, sim_score, 1);
  }
  if(PRINT_IN_MCL_FORMAT) {
    if(DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT)
      string[inpa_number]->insert(world_index, sim_score, div_factor);
    else       string[inpa_number]->insert(world_index, sim_score, 1);
  }
}

/**
   @Name: insert_inpa(..) Inserts inparalogs for pairwise output: Does not check the uniqueness. (The data source to be used has unique values)
   @Changed: 07.01.2011 by oekseth
*/
void mcl_format::insert_inpa(uint world_index_in, uint world_index_out, char *name_in, char *name_out, float sim_score, float _div_factor) {
  assert(string);
  assert(string[inpa]);
  assert(string[inpa_number]);
  set_debug_pair_size(world_index_out, name_in, name_out, sim_score, inpa, inpa_number);
  float div_factor_mcl = _div_factor;
  float div_factor_abc = _div_factor;
  if(!DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT) div_factor_mcl = 1;
  if(!DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT) div_factor_abc = 1;
  if(PRINT_IN_ABC_FORMAT) {
    if(MODE_PAIRWISE_OUTPUT_ABC)
      string[inpa]->insert(name_in, name_out, sim_score, div_factor_abc);
    else
      string[inpa]->insert(name_out, sim_score, div_factor_abc);
  }
  if(PRINT_IN_MCL_FORMAT) {
    if(MODE_PAIRWISE_OUTPUT_MCL)
      string[inpa_number]->insert(world_index_in, world_index_out, sim_score, div_factor_mcl);
    else {
      string[inpa_number]->insert(world_index_out, sim_score, div_factor_mcl);
    }
  } //else       fprintf(stderr, "...not in mcl..\n");
}


//! Sets the name_index-thing
void mcl_format::set_name_index(char *name, uint world_index) {
  assert(name);
  const uint inserted_size = string[names_index]->set_name_index(name, world_index);
#ifndef NDEBUG
  //! Illustrates usage in other contexts:
  assert(inserted_size == mcl_bunch::get_size_name_index(name, world_index));
  assert(file_chunk);
  //! If in debug mode, appends the number of chars the name index corresponds to into the list:
  if(file_chunk) file_chunk->insert_size(names_index, inserted_size);
  const uint length = (uint)string[names_index]->get_current_size_string();
  uint calculated_size = 0;
  if(file_chunk) calculated_size = file_chunk->get_size(names_index);
  assert(length == calculated_size);
#endif
}

/**
   @Name: set_header_pair_ortho(..) -- Inserts header for ortho pairs.
   @Changed: 07.01.2011 by oekseth
*/
void mcl_format::set_header_pair_ortho(char *name, uint world_index) {
  //! Sets the header size:
  set_debug_header_size(world_index, name, pair_orth, pair_orth_number);
  if(PRINT_IN_ABC_FORMAT && !MODE_PAIRWISE_OUTPUT_ABC)    string[pair_orth]->set_header(name);
  if(PRINT_IN_MCL_FORMAT && !MODE_PAIRWISE_OUTPUT_MCL) {
    string[pair_orth_number]->set_header(world_index);
  }
}

/**
   @Name: insert_pair_ortho --  Inserts ortho pairs: Does not check the uniqueness. (The data source to be used has unique values)
   @Changed: 14.02.2011 by oekseth
*/
void mcl_format::insert_pair_ortho(uint world_index_in, uint world_index_out, char *name_in, char *name_out, float sim_score, float _div_factor) {
  set_debug_pair_size(world_index_out, name_in, name_out, sim_score, pair_orth, pair_orth_number);
  assert(string);
  assert(string[pair_orth]);
  assert(string[pair_orth_number]);
  float div_factor_mcl = _div_factor;
  float div_factor_abc = _div_factor;
  if(!DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT) div_factor_mcl = 1;
  if(!DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT) div_factor_abc = 1;

  if(PRINT_IN_ABC_FORMAT) {
    if(MODE_PAIRWISE_OUTPUT_ABC) 	string[pair_orth]->insert(name_in, name_out, sim_score, div_factor_abc);
    else 	string[pair_orth]->insert(name_out, sim_score, div_factor_abc);
  }
  if(PRINT_IN_MCL_FORMAT) {
    if(MODE_PAIRWISE_OUTPUT_MCL) 	string[pair_orth_number]->insert(world_index_in, world_index_out, sim_score, div_factor_mcl);
    else                              string[pair_orth_number]->insert(world_index_out, sim_score, div_factor_mcl);
  }
}


/**
   @brief Inserts orthologs into both the 'name' and 'number' strings: Checks uniqueness.
   @return '1' if an eleemnt was inserted, else '0'.
   @remarks The return value is used in debug mode in class pipe_struct verifying that the number of inserted realtions correspons to our expectations.
*/  
uint mcl_format::insert_ortho_inpa(char *name, uint world_index, float sim_score, float _div_factor) {
  assert(string);
  assert(string[orth_inpa]);
  assert(string[orth_inpa_number]);
  pair<set<uint>::iterator,bool> ret =  set_orth_inpa.insert(world_index); // set do not allow unequal keys
  if(ret.second) { // set to true if a new element was inserted
    //! If in debug mode, sums the sizes of the estimated file size:
    set_debug_pair_size(world_index, name, sim_score, orth_inpa, orth_inpa_number);
    //! Inserts the actual data:
    float div_factor_mcl = _div_factor;
    float div_factor_abc = _div_factor;
    if(!DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT) div_factor_mcl = 1;
    if(!DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT) div_factor_abc = 1;
    if(PRINT_IN_ABC_FORMAT) string[orth_inpa]->insert(name, sim_score, div_factor_abc);
    if(PRINT_IN_MCL_FORMAT) string[orth_inpa_number]->insert(world_index, sim_score, div_factor_mcl);
    return 1; // One element was inserted.
  } else return 0; // No elements was inserted
}

/**
   @Name: set_header_ortho_inpa(..) -- Inserts header for orth_inpa
   @Changed: 07.01.2011 by oekseth
*/
void mcl_format::set_header_ortho_inpa(char *name, uint world_index) {
  //! Sets the header size:
  set_debug_header_size(world_index, name, orth_inpa, orth_inpa_number);
  assert(string);
  assert(string[orth_inpa]);
  assert(string[orth_inpa_number]);

  if(PRINT_IN_ABC_FORMAT && !MODE_PAIRWISE_OUTPUT_ABC) {
    string[orth_inpa]->set_header(name);
  }
  if(PRINT_IN_MCL_FORMAT && !MODE_PAIRWISE_OUTPUT_MCL)
    string[orth_inpa_number]->set_header(world_index);
}

/**
   @brief Inserts ortho pairs: Does not check the uniqueness. (The data source to be used has unique values)
   @return '1' if an eleemnt was inserted, else '0'.
   @remarks The return value is used in debug mode in class pipe_struct verifying that the number of inserted realtions correspons to our expectations.
**/  
uint mcl_format::insert_ortho_inpa(uint world_index_in, uint world_index_out, char *name_in, char *name_out, float sim_score, float _div_factor) {
  assert(string);
  assert(string[orth_inpa]);
  assert(string[orth_inpa_number]);

  pair<set<uint>::iterator,bool> ret =  set_orth_inpa.insert(world_index_out); // set do not allow unequal keys
  float div_factor_mcl = _div_factor;
  float div_factor_abc = _div_factor;
  if(!DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT) div_factor_mcl = 1;
  if(!DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT) div_factor_abc = 1;
  
  if(ret.second) { // set to true if a new element was inserted
    set_debug_pair_size(world_index_out, name_in, name_out, sim_score, orth_inpa, orth_inpa_number);
    if(PRINT_IN_ABC_FORMAT)  {
      if(MODE_PAIRWISE_OUTPUT_ABC) 
	string[orth_inpa]->insert(name_in, name_out, sim_score, div_factor_abc);
      else string[orth_inpa]->insert(name_out, sim_score, div_factor_abc);
    }
    if(PRINT_IN_MCL_FORMAT)  {
      if(MODE_PAIRWISE_OUTPUT_MCL) 
	string[orth_inpa_number]->insert(world_index_in, world_index_out, sim_score, div_factor_mcl);
      else string[orth_inpa_number]->insert(world_index_out, sim_score, div_factor_mcl);
    }
    return 1; // One element was inserted.
  } else return 0; // No elements was inserted
}

/**!
   @Name: set_line_end(..) -- Sets the line end, and updates the 'all' and 'all_number' strings
   @Changed: 07.01.2011 by oekseth (init)
*/
void mcl_format::set_line_end() {
  assert(string);
  assert(string);
  assert(string[all]);
  assert(string[all_number]);
//   if(valgrind_try_provoking_uninitialised_values()) {
//     fprintf(stdout, "!!\t Values not initialised, at mcl_format:%d\n", __LINE__);
//   }
  bool co_orthologs_has_data = false;
  if(string[orth_inpa]->get_size_of_line() > 0 || string[orth_inpa_number]->get_size_of_line() > 0)
    co_orthologs_has_data = true;
  if(string[inpa_number]->get_size_of_line() > 0 || string[inpa]->get_size_of_line() > 0) {
    string[all]->copy_line(*string[inpa], *string[pair_orth], *string[orth_inpa]);
    string[all_number]->copy_line(*string[inpa_number], *string[pair_orth_number], *string[orth_inpa_number]);
  } else if (string[pair_orth]->get_size_of_line() > 0 || string[pair_orth_number]->get_size_of_line() > 0) {
    string[all]->copy_line(*string[pair_orth], *string[orth_inpa], *string[inpa]);
    string[all_number]->copy_line(*string[pair_orth_number], *string[orth_inpa_number], *string[inpa_number]);
  } else if(co_orthologs_has_data) {
    string[all]->copy_line(*string[orth_inpa], *string[pair_orth],  *string[inpa]);
    string[all_number]->copy_line(*string[orth_inpa_number], *string[pair_orth_number],  *string[inpa_number]);    
  } else {
    set_orth_inpa.clear(); // clears the list
    return;
  }
#ifndef NDEBUG
  if(file_chunk) {
    file_chunk->reset_size(all, string[all]->get_current_size_string());
    file_chunk->reset_size(all_number, string[all_number]->get_current_size_string());
  }
#endif
  //! Note: As the "file_chunk" object is still in the developmental mode, it's not included in the "live run"
  
//   Debug_update_size_list_with_line_end(all, string[all]->get_current_size_string());
//   debug_update_size_list_with_line_end(all_number, string[all_number]->get_current_size_string());
  //3. Add line ends for the reaminder
  assert(string[inpa]);
  assert(string[inpa_number]);
  assert(string[pair_orth]);
  assert(string[pair_orth_number]);
  assert(string[orth_inpa]);
  assert(string[orth_inpa_number]);
  const bool def = true;
  int positions_inserted = 0;
  positions_inserted = string[inpa]->set_line_end(false, def);
#ifndef NDEBUG
  if(file_chunk) file_chunk->insert_size(inpa, positions_inserted);
#endif

  positions_inserted = string[pair_orth]->set_line_end(false, def);
#ifndef NDEBUG
  if(file_chunk) file_chunk->insert_size(pair_orth, positions_inserted);
#endif

  positions_inserted = string[pair_orth_number]->set_line_end(true, def);
#ifndef NDEBUG
  if(file_chunk) file_chunk->insert_size(pair_orth_number, positions_inserted);
#endif

  positions_inserted = string[orth_inpa_number]->set_line_end(true, def);
#ifndef NDEBUG
  if(file_chunk) file_chunk->insert_size(orth_inpa_number, positions_inserted);
#endif

  positions_inserted = string[all_number]->set_line_end(true, false);
#ifndef NDEBUG
  if(file_chunk) file_chunk->insert_size(all_number, positions_inserted);
#endif
  set_orth_inpa.clear(); // clears the list
}


/**
   @return the file names
   @remarks The naming convention follows:
   - Names with integers (numbers in this context) corresponds to mcl type, thereby .mcl'
   - Those files with strings (those other in this contect), the naming covention '.abc' is used
**/
char *mcl_format::get_file_name(uint id) {
  return result_format::get_file_name_result(id); // Defined in file 'enum_mcl.h'
}
#ifdef USE_MPI

/**!
   @Name: file_open(..) -- Opens the file given the identifier, describing the purpose of the file
*/
MPI_File *mcl_format::file_open(char *identifier, char *FILE_BINARY_LOCATION) {
  char *output_file = new char[100];
  assert(output_file);
  memset(output_file, '\0', 100);
  if(FILE_BINARY_LOCATION != NULL && strlen(FILE_BINARY_LOCATION)>2) {
    if(FILE_BINARY_LOCATION[strlen(FILE_BINARY_LOCATION)-1] != '/')
      sprintf(output_file, "%s/%s", FILE_BINARY_LOCATION, identifier);
    else sprintf(output_file, "%s%s", FILE_BINARY_LOCATION, identifier);
  } else sprintf(output_file, "%s", identifier);
  //! Opens the file:
  MPI_File* out_file = new MPI_File;
  assert(out_file);
  MPI_File_open(MPI_COMM_WORLD, output_file, MPI_MODE_CREATE | MPI_MODE_WRONLY,
		MPI_INFO_NULL, out_file);
  //! Some mpi variables:
  int myrank = 0;  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
  /* Set the file view */
  // TODO: The below ode not removed, in order to serve as template if more effective wrting is choosen:
  //  MPI_Offset disp = myrank * 10; // TODO: This value is probably wrong; should be replaced when understanding is increased.
  //  MPI_Offset disp = myrank * 1000*100; // TODO: This value is probably wrong; should be replaced when understanding is increased.
  //  MPI_File_set_view(*out_file, disp, MPI_CHAR, MPI_CHAR,  		    "native", MPI_INFO_NULL);

  if(!out_file) {
    perror(output_file);
    exit(2);
  }
  delete [] output_file;
  return out_file;
}


#else
/**!
   @Name: file_open(..) -- Opens the file given the identifier, describing the purpose of the file
   @Changed: 08.01.2011 by oekseth
*/
FILE *mcl_format::file_open(char *identifier, char *FILE_BINARY_LOCATION) {
  char *output_file = result_format::get_storage_path_for_given_file_id(identifier, FILE_BINARY_LOCATION);
  FILE* out_file = fopen(output_file, "w");
  if(!out_file) {
    perror(output_file);
    exit(2);
  }
  delete [] output_file;
  return out_file;
}
#endif

/**! @return a list of open files */
#ifdef USE_MPI
MPI_File **mcl_format::init_file_array(mcl_t type, char *FILE_BINARY_LOCATION, bool PRINT_IN_ABC_FORMAT, bool PRINT_IN_MCL_FORMAT) 
#else
  FILE **mcl_format::init_file_array(mcl_t type, char *FILE_BINARY_LOCATION, bool PRINT_IN_ABC_FORMAT, bool PRINT_IN_MCL_FORMAT) 
#endif
{
#ifdef USE_MPI
  MPI_File ** file = new MPI_File*[mcl_t_size];
#else
  FILE ** file = new FILE*[mcl_t_size];
#endif
  for(uint i =0; i<mcl_t_size; i++) file[i] = NULL; // Initilies

  if(true) { // TODO: consider having a user option here.
    const uint index = names_index;
    if(index != (uint)type) {
      file[index] = file_open(get_file_name(index), FILE_BINARY_LOCATION);//fopen(file_name, "wb");
    } else {
#ifdef USE_MPI
      file[index] = file_open(get_file_name(index), FILE_BINARY_LOCATION);//fopen(file_name, "wb");
      log_builder::throw_warning(input_terminal, __LINE__, __FILE__, __FUNCTION__,"Option printing data to terminal de-activated when MPI is used (instead the default is used, ie, data is written to a file)");
#else
      file[index] = stdout;
#endif
    }
  }
  for(uint i =0; i<mcl_list_size; i++) {
    // The formal header for the mcl format (the 'numbers' in this context)
    if(PRINT_IN_MCL_FORMAT) {
      const uint index = list_mcl_number[i];
      if(index != (uint)type) {
	file[index] = file_open(get_file_name(index), FILE_BINARY_LOCATION);//fopen(file_name, "wb");
      } else {
#ifdef USE_MPI
	file[index] = file_open(get_file_name(index), FILE_BINARY_LOCATION);//fopen(file_name, "wb");
	log_builder::throw_warning(input_terminal, __LINE__, __FILE__, __FUNCTION__,"Option printing data to terminal de-activated when MPI is used (instead the default is used, ie, data is written to a file)");
#else
	file[index] = stdout;
#endif
      }
    }
    if(PRINT_IN_ABC_FORMAT) {
      const uint index = list_mcl[i];
      if(index != (uint)type) {
	file[index] = file_open(get_file_name(index), FILE_BINARY_LOCATION);//fopen(file_name, "wb");
      } else {
#ifdef USE_MPI
	file[index] = file_open(get_file_name(index), FILE_BINARY_LOCATION);//fopen(file_name, "wb");
	log_builder::throw_warning(input_terminal, __LINE__, __FILE__, __FUNCTION__,"Option printing data to terminal de-activated when MPI is used (instead the default is used, ie, data is written to a file)");
#else
	file[index] = stdout;
#endif
      }
    }
  }
  return file;
}
/**
   @Name: close_init_file_array -- Closes a list of open files
   @Changed: 07.01.2011 by oekseth
*/
#ifdef USE_MPI
void mcl_format::close_init_file_array(MPI_File **&file, mcl_t type, bool SORT_ABC_DATA, char *FILE_BINARY_LOCATION) 
{
  //! Only one node is set to write trailing chars:
  bool append_chars = false;
  int number_of_nodes = 1; MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes); 
  int myrank = 0;  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
  //! The 'root' node is given responsibility writing the header file
  if(number_of_nodes > 1 && myrank == 0) {
    append_chars = true;
  } else if(number_of_nodes == 1) append_chars = true;
  if(file != NULL) {
    const uint index = names_index;
    //! Closes the protein-map:
    if(index != (uint)type && file[index]) {
      MPI_File_close(file[index]);
      delete file[index]; file[index] = NULL;
    }

    //! Closes the set of data files:
    for(uint i =0; i<mcl_list_size; i++) {
      if(append_chars) {
	if(file[list_mcl_number[i]]!= NULL){
	  char str[10]; memset(str, '\0', 10); sprintf(str, "%c", ')');
	  mcl_bunch::write_string_to_file(file[list_mcl_number[i]], str);
	}
	if(file[list_mcl_number[i]]!= NULL) {
	  char str[10]; memset(str, '\0', 10); sprintf(str, "%c", '\n');
	  mcl_bunch::write_string_to_file(file[list_mcl_number[i]], str);
	  //	fputc('\n', file[list_mcl_number[i]]);
	}
      }
      if(file[list_mcl_number[i]] != NULL) {
	MPI_File_close(file[list_mcl_number[i]]);
	delete file[list_mcl_number[i]]; file[list_mcl_number[i]] = NULL; 
      }
      if(file[list_mcl[i]]!= NULL) {
	MPI_File_close(file[list_mcl[i]]);
	delete file[list_mcl[i]]; file[list_mcl[i]] = NULL;
	//	fclose(file[list_mcl[i]]);
	if(SORT_ABC_DATA && append_chars) {
	  char *identifier = get_file_name(list_mcl[i]);
	  assert(identifier);
	  char *output_file = new char[100];
	  assert(output_file);
	  if(FILE_BINARY_LOCATION != NULL && strlen(FILE_BINARY_LOCATION)>2) {
	    if(FILE_BINARY_LOCATION[strlen(FILE_BINARY_LOCATION)-1] != '/')
	      sprintf(output_file, "%s/%s", FILE_BINARY_LOCATION, identifier);
	    else sprintf(output_file, "%s%s", FILE_BINARY_LOCATION, identifier);
	  } else sprintf(output_file, "%s", identifier);
	  char file_n[400]; for(uint i =0; i<400;i++) file_n[i]='\0';
	  sprintf(file_n, "sort %s -n -k3 -o %s.s", output_file, output_file);
	  assert(-1 != system(file_n));
	  remove(output_file); delete [] output_file;
	}
      }
      
    }
    delete [] file; file = NULL;
  }
}
#else
void mcl_format::close_init_file_array(FILE **&file, mcl_t type, bool SORT_ABC_DATA, char *FILE_BINARY_LOCATION) {
  if(file != NULL) {
    const uint index = names_index;
    //! Closes the protein-map:
    if(index != (uint)type && file[index]) {
      fclose(file[index]); file[index] = NULL;
    }

    for(uint i =0; i<mcl_list_size; i++) {
      if(file[list_mcl_number[i]]!= NULL)      fputc(')', file[list_mcl_number[i]]); 
      if(file[list_mcl_number[i]]!= NULL)      fputc('\n', file[list_mcl_number[i]]); 
      if(type!= list_mcl_number[i] && (file[list_mcl_number[i]]!= NULL) ) fclose(file[list_mcl_number[i]]);
      if(type!= list_mcl[i] && (file[list_mcl[i]]!= NULL)) {
	fclose(file[list_mcl[i]]);
	if(SORT_ABC_DATA) {
	  char *identifier = get_file_name(list_mcl[i]);
	  char *output_file = new char[100];
	  if(FILE_BINARY_LOCATION != NULL && strlen(FILE_BINARY_LOCATION)>2) {
	    if(FILE_BINARY_LOCATION[strlen(FILE_BINARY_LOCATION)-1] != '/')
	      sprintf(output_file, "%s/%s", FILE_BINARY_LOCATION, identifier);
	    else sprintf(output_file, "%s%s", FILE_BINARY_LOCATION, identifier);
	  } else sprintf(output_file, "%s", identifier);
	  char file_n[400]; for(uint i =0; i<400;i++) file_n[i]='\0';
	  //sprintf(file_n, "/usr/bin/sort %s -n -k3 -o %s.s", output_file, output_file);
	  sprintf(file_n, "sort %s -n -k3 -o %s.s", output_file, output_file);
	  const bool ret_val_sys = system(file_n);
	  assert(-1 != ret_val_sys);
	  char assert_sys[400]; sprintf(assert_sys, "test -x %s.s", output_file);
	  const int ret_val = system(assert_sys);
	  if(ret_val == -1) {
	    printf("\t Retval after sorting operation '%s' was %d, at [%s]:%s:%d\n", file_n, ret_val, __FUNCTION__, __FILE__, __LINE__);
	  }
	  remove(output_file); 
	  delete [] output_file;
	}
      }
    }
    delete [] file; file = NULL;
  }
}
#endif

/**
   @brief Creates the headers, and writes them to the files 
   @return a checksum corresponding to the total number of chars inserted
*/
#ifdef USE_MPI
uint mcl_format::write_file_headers(MPI_File **file, taxa_t *listTaxa, int taxon_length, bool PRINT_IN_MCL_FORMAT)
#else
  uint mcl_format::write_file_headers(FILE **file, taxa_t *listTaxa, int taxon_length, bool PRINT_IN_MCL_FORMAT)
#endif
{
  uint cnt_chars_inserted = 0; 
  if(file != NULL) {
    const uint size_proteins =  listTaxa[taxon_length-1].rel_end;
    if(PRINT_IN_MCL_FORMAT) {
      cnt_chars_inserted = result_format::get_length_of_mcl_header(size_proteins);
      for(uint i =0; i<mcl_list_size; i++) {
	//! The formal header for the mcl format (the 'numbers' in this context)
	char string[100]; memset(string, '\0', 100);
	//	const char *arr_str = "(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n(mcldoms\n";
	const char *arr_str = result_format::get_mcl_header_format();
	sprintf(string, arr_str, size_proteins, size_proteins);
	mcl_bunch::write_string_to_file(file[list_mcl_number[i]], string);
	// Do not add a header for the abc files, but reained (below) if future needs change:
	//      fprintf(file[list_mcl[i]],        arr_str, size_proteins, size_proteins);
      }

#ifndef NDEBUG
      {
	const uint expected_size =  result_format::get_length_of_mcl_header(size_proteins);
	const uint size_each_file = cnt_chars_inserted;
	assert(size_each_file == expected_size);
      }
#endif

      // Builds a char list of the data to be written to the numbers:
      const uint length_key_names_number = (uint)log10(UINT_MAX)+1; // the maximum numbers of chars it can hold
      const uint size_str_number = length_key_names_number * size_proteins;
      char *str_number = new char[size_str_number];
      assert(str_number);
      if(str_number != NULL) {
	memset(str_number, '\0', size_str_number);
	uint str_pos = 0; // the position in the string to add data
	for(uint index_in = 0;index_in < size_proteins;index_in++)  {// iterating thworugh every protein
	  sprintf(str_number + str_pos, "%u%c ", index_in, MCL_HEAD_SEPERATOR);
	  uint str_len = 1;
	  if(index_in > 0) str_len = (uint)log10(index_in)+1;
	  str_pos += str_len+1; // for the space
	
	}
	//	sprintf(str_number + str_pos, result_format::get_mcl_header_format_before_main());
	char *string = result_format::get_mcl_header_format_before_main();
	assert(string);
	assert(strlen(string)>0);
	strncpy(str_number + str_pos, string, strlen(string));
	assert(strlen(str_number + str_pos) == strlen(result_format::get_mcl_header_format_before_main()));
#ifndef NDEBUG
	{
	  const uint expected_size =  list_file_chunk::get_expected_size_of_header(listTaxa, taxon_length);
	  const uint size_each_file = cnt_chars_inserted + strlen(str_number);
	  if(size_each_file != expected_size) {
	    printf("test:\texpected_size(%u) =?= size_each_file(%u), with str_pos(%u) at %s:%d\n", expected_size, size_each_file, str_pos, __FILE__, __LINE__); 
	    assert(size_each_file == expected_size);
	  }
	}
	const uint cnt_chars_inserted_before_meta = cnt_chars_inserted;
	const uint expected_size_of_this = mcl_list_size*(cnt_chars_inserted + strlen(str_number));
#endif
	str_pos += 22;//
	cnt_chars_inserted = mcl_list_size*(cnt_chars_inserted + strlen(str_number));
	//	cnt_chars_inserted += strlen(str_number);
	for(uint i =0; i<mcl_list_size; i++) {// Writes the data to the number files
	  mcl_bunch::write_string_to_file(file[list_mcl_number[i]], str_number);
	}

#ifndef NDEBUG
	{
	  //! Verifies our assumption:
	  list_file_chunk temp_obj =  list_file_chunk();
	  const uint expected_size =  temp_obj.append_header_sizes(listTaxa, taxon_length);
	  const uint size_each_file = cnt_chars_inserted;
	  if(size_each_file != expected_size) {
	    printf("test:\texpected_size(%u) =?= size_each_file(%u), with cnt_chars_inserted_before_meta(%u) and expected_size_of_this(%u) at %s:%d\n", expected_size, size_each_file, cnt_chars_inserted_before_meta, expected_size_of_this, __FILE__, __LINE__); 
	    assert(size_each_file == expected_size);
	  }
	  temp_obj.free_memory();
	}
#endif	
	if(str_number) {delete [] str_number; str_number = NULL;}

	// Code below made inaccessible at January 17 out the header for the 'abc' format, as the data is not inteded ot be used as input for the mcl directly 
	if(false) {
	  //
	  // Then the list of the given proteins:
	  if(listTaxa != NULL) {
	    const uint length_key_names = 40; // a maximum upper limit
	    const uint size_string_names = length_key_names * size_proteins;
	    char *string_names = new char[size_string_names];
	    char *string_start = string_names;
	    if(string_names != NULL) { // has enough memory
	      memset(string_names, '\0', size_string_names);
	      for(uint i =0; i < (uint)taxon_length; i++) {
		string_names = listTaxa[i].mcl_get_all_header_proteins(string_names);
	      }
	      if((string_names - string_start) > 0) { // data added, no error occured
		*string_names = ')';  string_names++;// Adds the end of the string
		*string_names = '\n';  string_names++;
		*string_names = '(';  string_names++;
		*string_names = '\n';  string_names++;
		*string_names = '\0';
		for(uint i =0; i<mcl_list_size; i++) {// Writes the data to the number files
		  mcl_bunch::write_string_to_file(file[list_mcl[i]], string_start);
		  //		  fputs(string_start, file[list_mcl[i]]);
		}
	      } else fprintf(stderr, "!!\tNo data added. An error occured.\n");
	      delete [] string_start;
	    }  else fprintf(stderr, "!!\tNot enough memory for writing the headers to a name file. An error occured at line %d in class mcl_format. Contact the developer. \n", __LINE__);
	  }
	  else fprintf(stderr, "!!\tListTaxa not set. Wrong usage of the code, or the preceding elements of it. Error at line %d in class mcl_format. Contact the developer for questions.\n", __LINE__);
	}
      } else fprintf(stderr, "!!\tNot enough memory for writing the headers to a number file. An error occured at line %d in class mcl_format. Contact the devloper\n", __LINE__);
    } //else fprintf(stderr, "!!\tFile lsit not initilixed. An error occured, due to wrong usage of this function\n");
  }
 return cnt_chars_inserted;
}

/**
   @Name: write_file(..) -- Writes the data to the files, given the data set.
   @Is_paralell: Negative, only to be run in a serial mode.
*/
#ifdef USE_MPI
void mcl_format::write_file(MPI_File **file) 
#else
  void mcl_format::write_file(FILE **file) 
#endif
{
  if(string!=NULL && file != NULL) {    
    for(uint i =0; i<mcl_list_size; i++) { // Evaluates for the indexes (numbers):
      // Evaluates for the numbers (protein names printed as integers):
      const uint index_mcl = list_mcl_number[i];
      if(file[index_mcl] != NULL) { 
	if(!MODE_PAIRWISE_OUTPUT_MCL) { // With this option, the data are concatenated
	  if(file[index_mcl] != NULL) string[index_mcl]->write_data_to_file(file[index_mcl]);
	} else { // Build the 'all' file directly form these datas	  
	  if (index_mcl == all_number ) { // Writes everything to the all file
	    string[inpa_number]->write_data_to_file(file[index_mcl]);
	    string[pair_orth_number]->write_data_to_file(file[index_mcl]);	    
	    string[orth_inpa_number]->write_data_to_file(file[index_mcl]);
	  } else string[index_mcl]->write_data_to_file(file[index_mcl]);
	}
      }
      // Evaluates for the labels (protein names printed as strings):
      const uint index = list_mcl[i];
      if(file[index] != NULL) { 
	if(!MODE_PAIRWISE_OUTPUT_ABC) { 
	  string[index]->write_data_to_file(file[index]);
	} else {
	  if(index == all) {// Writes everything to the 'all' 
	    string[inpa]->write_data_to_file(file[index]);
	    string[pair_orth]->write_data_to_file(file[index]);
	    string[orth_inpa]->write_data_to_file(file[index]);
	  } else {
	    const uint ret_val = string[index]->write_data_to_file(file[index]);
#ifndef NDEBUG
	    if(file_chunk) {
	      assert(ret_val == file_chunk->get_length(index));
	    }
#endif
	  }
	}
      
      }
    }
    if(file[names_index] != NULL) { // The file holding info of label-interg-data.
      string[names_index]->write_data_to_file(file[names_index]);
    }      
  }
}

mcl_format *mcl_format::init(int taxon_length, taxa_t *listTaxa, bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT, 
			     bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT, bool MODE_PAIRWISE_OUTPUT_ABC,  
			     bool MODE_PAIRWISE_OUTPUT_MCL, bool PRINT_IN_ABC_FORMAT, bool PRINT_IN_MCL_FORMAT, 
			     bool SORT_ABC_DATA) {
  return new mcl_format(taxon_length, listTaxa, DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT, DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT,
			MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT,
			SORT_ABC_DATA);
}
// Class constructor:
mcl_format::mcl_format(  int _taxon_length,
			 taxa_t *_listTaxa, 
			 bool _DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT, 
			 bool _DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT, 
			 bool _MODE_PAIRWISE_OUTPUT_ABC,  // for the abc file: if set, the data out pairwise in stead of as in a row
			 bool _MODE_PAIRWISE_OUTPUT_MCL,  // for the mcl file: if set, the data out pairwise in stead of as in a row
			 // The following variables decides what output is to be printed
			 bool _PRINT_IN_ABC_FORMAT, 
			 bool _PRINT_IN_MCL_FORMAT, 
			 bool _SORT_ABC_DATA  // if true, sorts teh abc files before outprint)
			 ) :
  taxon_length(_taxon_length), listTaxa(_listTaxa), 
  DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT(_DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT), 
  DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT(_DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT), 
  MODE_PAIRWISE_OUTPUT_ABC(_MODE_PAIRWISE_OUTPUT_ABC),  // for the abc file: if set), the data out pairwise in stead of as in a row
  MODE_PAIRWISE_OUTPUT_MCL(_MODE_PAIRWISE_OUTPUT_MCL),  // for the mcl file: if set), the data out pairwise in stead of as in a row
  // The following variables decides what output is to be printed
  PRINT_IN_ABC_FORMAT(_PRINT_IN_ABC_FORMAT), 
  PRINT_IN_MCL_FORMAT(_PRINT_IN_MCL_FORMAT), 
  SORT_ABC_DATA(_SORT_ABC_DATA)  // if true), sorts teh abc files before outprint)
#ifndef NDEBUG
  , file_chunk(NULL)
#endif
{      
  init_string();
#ifndef NDEBUG
  file_chunk = new list_file_chunk(MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT, SORT_ABC_DATA);
#endif
}
mcl_format::mcl_format(mcl_print_settings_t print_settings) 
#ifndef NDEBUG
  : file_chunk(NULL)
#endif
{      
  init_string();
  taxon_length = print_settings.taxon_length, listTaxa = print_settings.listTaxa;
  DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT = print_settings.DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT;
  DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT = print_settings.DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT;
  MODE_PAIRWISE_OUTPUT_ABC = print_settings.MODE_PAIRWISE_OUTPUT_ABC; // for the abc file: if set, the data out pairwise in stead of as in a row
  MODE_PAIRWISE_OUTPUT_MCL = print_settings.MODE_PAIRWISE_OUTPUT_MCL; // for the mcl file: if set, the data out pairwise in stead of as in a row
  PRINT_IN_ABC_FORMAT  = print_settings.PRINT_IN_ABC_FORMAT;
  PRINT_IN_MCL_FORMAT = print_settings.PRINT_IN_MCL_FORMAT;
  SORT_ABC_DATA = print_settings.SORT_ABC_DATA; // if true, sorts teh abc files before outprint
#ifndef NDEBUG
  file_chunk = new list_file_chunk(MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT, SORT_ABC_DATA);
#endif
}

#ifdef assert_code
//! Prints the elements found in the 'orth_inpa list
void mcl_format::print_ortho_inpa_line_list() {
  for (set<uint>::iterator it=set_orth_inpa.begin(); it!=set_orth_inpa.end(); it++) {
    printf("world_index(%u)\n", *it);
  }
}

//! Compares the strings
bool mcl_format::compare_string(char *inp, const bool print_diff, mcl_t type) {
  return string[type]->compare_string(inp, print_diff);
}

#endif
void mcl_format::assert_class(const bool print_info) {
  const static char *class_name = "mcl_format";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  mcl_print_settings_t print_settings;
  // Permutation 1: SEttings for easy visual validation:
  print_settings.DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT = false;
  print_settings.DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT = false;
  print_settings.PRINT_IN_ABC_FORMAT = true,   print_settings.PRINT_IN_MCL_FORMAT = true;
  print_settings.MODE_PAIRWISE_OUTPUT_MCL = false, print_settings.MODE_PAIRWISE_OUTPUT_ABC = false;
  class mcl_format arg = mcl_format(print_settings);
  arg.set_header_ortho_inpa("ole", 32);
  assert(arg.compare_string("ole\t", true, orth_inpa));

  assert(arg.compare_string("32\t", true, orth_inpa_number));   
  arg.insert_ortho_inpa("hans", 50, 4, 1);
  assert(arg.compare_string("32\t50:4.000", true, orth_inpa_number));   

  arg.insert_ortho_inpa("hans", 50, 4, 1); // verifies that duplictes are not inserted
  assert(arg.compare_string("32\t50:4.000", true, orth_inpa_number));
  assert(arg.compare_string("ole\thans:4.000", true, orth_inpa));
  arg.set_line_end();
  assert(arg.compare_string("ole\thans:4.000\t\n", true, all));
  
  arg.free_mem();
  
  arg = mcl_format(print_settings);
  arg.set_header_ortho_inpa("ole", 32);
  assert(arg.compare_string("ole\t", false, orth_inpa));
  assert(arg.compare_string("32\t", true, orth_inpa_number));   
  arg.insert_ortho_inpa("hans", 50, 4, 1);
  assert(arg.compare_string("32\t50:4.000", true, orth_inpa_number));   
  
  arg.insert_ortho_inpa("hans", 50, 4, 1); // verifies that duplictes are not inserted
  assert(arg.compare_string("32\t50:4.000", true, orth_inpa_number));
  assert(arg.compare_string("ole\thans:4.000", false, orth_inpa));
  arg.set_header_inpa(10, "miron");
  arg.insert_inpa(43, "ole", 5.0, 1);
  arg.insert_inpa(44, "ole2", 15.0, 1);
  
  //  assert(arg.compare_string("miron\tole:5.000\tole2:15.000\t", true, inpa));  
  assert(arg.compare_string("miron\tole:5.000\tole2:15.000\thans:4.000\t$", true, inpa));  
  arg.insert_inpa(45, "ole3", 25.0, 1);
  arg.set_line_end();
  char *string = "miron\tole:5.000\tole2:15.000\tole3:25.000\thans:4.000\t";  
  assert(arg.compare_string(string, true, all));
  arg.free_mem();

  // Tests using pairwise out-print:
  print_settings.MODE_PAIRWISE_OUTPUT_MCL = true; print_settings.MODE_PAIRWISE_OUTPUT_ABC = true;
  arg = mcl_format(print_settings);

  arg.set_header_pair_ortho("ole", 32);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_ortho_inpa(0, 1, "ole", "one", 5.0, 1.0);
  arg.insert_ortho_inpa(0, 2, "ole", "two", 5.0, 1.0);
  arg.insert_pair_ortho( 0, 1, "ole", "one", 5.0, 1.0);

  arg.set_line_end();
  assert(arg.compare_string("ole\tone\t5.000\nole\ttwo\t5.000\n", true, orth_inpa));
  assert(arg.compare_string("0\t1\t5.000$\n0\t2\t5.000$\n", true, orth_inpa_number));
  assert(arg.compare_string("ole\tone\t5.000", true, pair_orth));
  assert(arg.compare_string("0\t1\t5.000$", true, pair_orth_number));
  assert(arg.compare_string("", true, all));
  assert(arg.compare_string("", true, all_number));
  arg.free_mem();



  //
  // Tests using row-wise out-print, with focus on building the 'all'
  print_settings.MODE_PAIRWISE_OUTPUT_MCL = false; print_settings.MODE_PAIRWISE_OUTPUT_ABC = false;
  arg = mcl_format(print_settings);

  arg.set_header_inpa(0, "ole");
  assert(arg.compare_string("", true, inpa));
  arg.insert_inpa(0, 1, "ole", "one", 5.0, 1.0);

  arg.set_header_ortho_inpa("ole", 0);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_ortho_inpa(0, 1, "ole", "one", 5.0, 1.0);

  arg.set_header_pair_ortho("ole", 0);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_pair_ortho(0, 1, "ole", "one", 5.0, 1.0);

  arg.set_line_end();
  assert(arg.compare_string("ole\tone:5.000\t\n", true, inpa));
  assert(arg.compare_string("0\t1:5.000\t$\n", true, inpa_number));
  assert(arg.compare_string("ole\tone:5.000\t\n", true, orth_inpa));
  assert(arg.compare_string("0\t1:5.000\t$\n", true, orth_inpa_number));
  assert(arg.compare_string("ole\tone:5.000\t\n", true, pair_orth));
  assert(arg.compare_string("0\t1:5.000\t$\n", true, pair_orth_number));
  //  arg.print(all);
  assert(arg.compare_string("ole\tone:5.000\tone:5.000\tone:5.000\t\n", true, all));
  assert(arg.compare_string("0\t1:5.000\t1:5.000\t1:5.000\t$\n", true, all_number));

  //
  // Tests using row-wise out-print, with focus on building the 'all', having multiple rows:
  arg.set_header_inpa(99, "hanne");
  assert(arg.compare_string("", true, inpa));
  arg.insert_inpa(99, 1, "hanne", "one", 5.0, 1.0);

  arg.set_header_ortho_inpa("hanne", 99);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_ortho_inpa(99, 1, "hanne", "one", 5.0, 1.0);

  arg.set_header_pair_ortho("hanne", 99);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_pair_ortho(99, 1, "hanne", "one", 5.0, 1.0);

  arg.set_line_end();
  //  arg.print(all);  arg.print(all_number);
  assert(arg.compare_string("ole\tone:5.000\t\nhanne\tone:5.000\t\n", true, inpa));
  assert(arg.compare_string("0\t1:5.000\t$\n99\t1:5.000\t$\n", true, inpa_number));
  assert(arg.compare_string("ole\tone:5.000\t\nhanne\tone:5.000\t\n", true, orth_inpa));
  assert(arg.compare_string("0\t1:5.000\t$\n99\t1:5.000\t$\n", true, orth_inpa_number));
  assert(arg.compare_string("ole\tone:5.000\t\nhanne\tone:5.000\t\n", true, pair_orth));
  assert(arg.compare_string("0\t1:5.000\t$\n99\t1:5.000\t$\n", true, pair_orth_number));
  //  arg.print(all);
  assert(arg.compare_string("ole\tone:5.000\tone:5.000\tone:5.000\t\nhanne\tone:5.000\tone:5.000\tone:5.000\t\n", true, all));
  assert(arg.compare_string("0\t1:5.000\t1:5.000\t1:5.000\t$\n99\t1:5.000\t1:5.000\t1:5.000\t$\n", true, all_number));

  arg.free_mem();


  //
  // Tests using row-wise out-print, with focus on building the 'all'.
  // -- Here, the ortho_inpa are set to empty:
  print_settings.MODE_PAIRWISE_OUTPUT_MCL = false; print_settings.MODE_PAIRWISE_OUTPUT_ABC = false;
  arg = mcl_format(print_settings);

  arg.set_header_inpa(0, "ole");
  assert(arg.compare_string("", true, inpa));
  arg.insert_inpa(0, 1, "ole", "one", 5.0, 1.0);
  
  arg.set_header_pair_ortho("ole", 0);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_pair_ortho(0, 1, "ole", "one", 5.0, 1.0);

  arg.set_line_end();
  assert(arg.compare_string("", true, orth_inpa));
  assert(arg.compare_string("", true, orth_inpa_number));
  assert(arg.compare_string("ole\tone:5.000\t\n", true, pair_orth));
  assert(arg.compare_string("0\t1:5.000\t$\n", true, pair_orth_number));
  assert(arg.compare_string("ole\tone:5.000\tone:5.000\t\n", true, all));
  assert(arg.compare_string("0\t1:5.000\t1:5.000\t$\n", true, all_number));
  arg.free_mem();

  //
  // Tests using row-wise out-print, with focus on building the 'all'.
  // -- Here, the pair_ortho are set to empty:
  print_settings.MODE_PAIRWISE_OUTPUT_MCL = false; print_settings.MODE_PAIRWISE_OUTPUT_ABC = false;
  arg = mcl_format(print_settings);

  arg.set_header_inpa(0, "ole");
  assert(arg.compare_string("", true, inpa));
  arg.insert_inpa(0, 1, "ole", "one", 5.0, 1.0);
  
  arg.set_header_ortho_inpa("ole", 0);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_ortho_inpa(0, 1, "ole", "one", 5.0, 1.0);
  
  arg.set_line_end();
  assert(arg.compare_string("ole\tone:5.000\t\n", true, orth_inpa));
  assert(arg.compare_string("0\t1:5.000\t$\n", true, orth_inpa_number));
  assert(arg.compare_string("", true, pair_orth));
  assert(arg.compare_string("", true, pair_orth_number));
  assert(arg.compare_string("ole\tone:5.000\tone:5.000\t\n", true, all));
  assert(arg.compare_string("0\t1:5.000\t1:5.000\t$\n", true, all_number));
  arg.free_mem();

  // Tests using row-wise out-print, with focus on building the 'all'.
  // -- Here, only ortho_inpa is set.
  print_settings.MODE_PAIRWISE_OUTPUT_MCL = false; print_settings.MODE_PAIRWISE_OUTPUT_ABC = false;
  arg = mcl_format(print_settings);
  assert(arg.compare_string("", true, inpa));
  
  arg.set_header_ortho_inpa("ole", 0);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_ortho_inpa(0, 1, "ole", "one", 5.0, 1.0);

  arg.set_line_end();
  assert(arg.compare_string("ole\tone:5.000\t\n", true, orth_inpa));
  assert(arg.compare_string("0\t1:5.000\t$\n", true, orth_inpa_number));
  assert(arg.compare_string("", true, pair_orth));
  assert(arg.compare_string("", true, pair_orth_number));
  assert(arg.compare_string("ole\tone:5.000\t\n", true, all));
  assert(arg.compare_string("0\t1:5.000\t$\n", true, all_number));
  arg.free_mem();

  //
  // Tests using row-wise out-print, with focus on building the 'all'.
  // -- Here, the inparalogs are set to empty:
  print_settings.MODE_PAIRWISE_OUTPUT_MCL = false; print_settings.MODE_PAIRWISE_OUTPUT_ABC = false;
  arg = mcl_format(print_settings);

  arg.set_header_ortho_inpa("ole", 0);
  arg.set_header_inpa(0, "ole");

  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_ortho_inpa(0, 1, "ole", "one", 5.0, 1.0);

  arg.set_header_pair_ortho("ole", 0);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_pair_ortho(0, 1, "ole", "one", 5.0, 1.0);

  arg.set_line_end();
  assert(arg.compare_string("ole\tone:5.000\t\n", true, orth_inpa));
  assert(arg.compare_string("0\t1:5.000\t$\n", true, orth_inpa_number));
  assert(arg.compare_string("ole\tone:5.000\t\n", true, pair_orth));
  assert(arg.compare_string("0\t1:5.000\t$\n", true, pair_orth_number));
  assert(arg.compare_string("ole\tone:5.000\tone:5.000\t\n", true, all));
  assert(arg.compare_string("0\t1:5.000\t1:5.000\t$\n", true, all_number));
  arg.free_mem();

  //
  // Tests using row-wise out-print, having multiple pairs:
  print_settings.MODE_PAIRWISE_OUTPUT_MCL = false; print_settings.MODE_PAIRWISE_OUTPUT_ABC = false;
  arg = mcl_format(print_settings);
  arg.set_header_ortho_inpa("ole", 0);

  arg.set_header_inpa(0, "ole");
  assert(arg.compare_string("", true, inpa));
  arg.insert_inpa(0, 1, "ole", "one", 5.0, 1.0);
  arg.insert_inpa(0, 2, "ole", "two", 5.0, 1.0);
  arg.insert_inpa(0, 3, "ole", "three", 5.0, 1.0);

  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_ortho_inpa(0, 1, "ole", "one", 5.0, 1.0);
  arg.insert_ortho_inpa(0, 2, "ole", "two", 5.0, 1.0);
  arg.insert_ortho_inpa(0, 3, "ole", "three", 5.0, 1.0);

  arg.set_header_pair_ortho("ole", 0);
  assert(arg.compare_string("", true, orth_inpa));
  arg.insert_pair_ortho(0, 1, "ole", "one", 5.0, 1.0);
  arg.insert_pair_ortho(0, 2, "ole", "two", 5.0, 1.0);
  arg.insert_pair_ortho(0, 3, "ole", "three", 5.0, 1.0);

  arg.set_line_end();
  assert(arg.compare_string("ole\tone:5.000\ttwo:5.000\tthree:5.000\t\n", true, inpa));
  assert(arg.compare_string("0\t1:5.000\t2:5.000\t3:5.000\t$\n", true, inpa_number));
  assert(arg.compare_string("ole\tone:5.000\ttwo:5.000\tthree:5.000\t\n", true, orth_inpa));
  assert(arg.compare_string("0\t1:5.000\t2:5.000\t3:5.000\t$\n", true, orth_inpa_number));
  assert(arg.compare_string("ole\tone:5.000\ttwo:5.000\tthree:5.000\t\n", true, pair_orth));
  assert(arg.compare_string("0\t1:5.000\t2:5.000\t3:5.000\t$\n", true, pair_orth_number));
  assert(arg.compare_string("ole\tone:5.000\ttwo:5.000\tthree:5.000\tone:5.000\ttwo:5.000\tthree:5.000\tone:5.000\ttwo:5.000\tthree:5.000\t\n", true, all));
  assert(arg.compare_string("0\t1:5.000\t2:5.000\t3:5.000\t1:5.000\t2:5.000\t3:5.000\t1:5.000\t2:5.000\t3:5.000\t$\n", true, all_number));
  arg.free_mem();
  // Tests the creation of pair lists:...

#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
  //  if(false)   print_constants();
}
