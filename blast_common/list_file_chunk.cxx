#include "list_file_chunk.h"

//! @return the size (length) given the index;
uint list_file_chunk::get_size(uint index) {
  if(!list) return 0;
  assert(index < mcl_t_size);
  return list[index];
}

//! Compare lists using the first list as basis, and assuming a length defined in file enum_mcl.h
bool list_file_chunk::compare_with_list(loint *two) {
  assert(list);
  assert(two);
  bool equal = true;
  for(uint i = 0; i < mcl_t_size; i++) {
    if(list[i] != (uint)two[i]) {
      fprintf(stderr, "%s\tfirst_element_size(%u) != other_element_size(%u), at %s:%d\n", result_format::get_file_name_result(i), list[i], (uint)two[i], __FILE__, __LINE__);
      equal = false;
      assert(list[i] == (uint)two[i]);
    }
  }
  return equal;
}

//! @return the file size, given the index:
uint list_file_chunk::get_file_size(uint index, char *FILE_BINARY_LOCATION) {
  assert(index < mcl_t_size);
  char *file_name = result_format::get_storage_path_for_given_file_id(result_format::get_file_name_result(index), FILE_BINARY_LOCATION, SORT_ABC_DATA);
  struct stat sb;
  if(stat(file_name, &sb) == -1) { 
#ifdef USE_MPI
    int myrank = 0;  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
    fprintf(stderr, "[myrank=%d]!!\tUnable to open file '%s' at %s:%d\n", myrank, file_name, __FILE__, __LINE__);
#else 
    fprintf(stderr, "!!\tUnable to open file '%s' at %s:%d\n", file_name, __FILE__, __LINE__);
#endif
    exit(EXIT_FAILURE);
  }
  const uint file_size = (uint)sb.st_size;

  if(file_name) {delete [] file_name; file_name = NULL;}
  return file_size;
}

/**
   @brief 
   @param <file_chunk_header> Combines the internal list with this parameter to produce the size for comparison.
   @return true if all sizes corresponds to expected file sizes.
**/
bool list_file_chunk::compare_sizes_with_files(list_file_chunk *file_chunk_header, char *FILE_BINARY_LOCATION) {
  //! How we expect this function to be used:
  assert(file_chunk_header);
  //! Containers for handling- and verification:
  uint *meta = file_chunk_header->get_list();
  bool is_equal = true;
  //! First sums the elements to define the expected size for "all.abc":
  if(list) {
    //! Below our assertion, corresponding to the one in class mcl_format
    if(!list[all_number]) {
      list[all_number] = list[inpa_number] + list[pair_orth_number] + list[orth_inpa_number];
    } else {
      list[all] = list[inpa] + list[pair_orth] + list[orth_inpa];
    }
  }
  //! The verifies all the elements versus the actual file:
  for(uint i = 0; i < mcl_t_size; i++) {
    //! Expectations:
    uint expected_size = 0;
    if(list) expected_size = list[i];
    if(meta) expected_size += meta[i];
    //! The file size:
    const uint file_size = get_file_size(i, FILE_BINARY_LOCATION);
    if(abs((int)expected_size - (int)file_size) > 3) { // Do not bother adding teh last trailing chars:
      fprintf(stderr, "%s\t!!\t expected_size(%u) != file_size(%u) in list_file_chunk:%d\n", result_format::get_file_name_result(i), expected_size, file_size, __LINE__);
      assert(expected_size == file_size);
      is_equal = false;
    }
  }
  return is_equal;
}

//! @return an initiated version of the array of containers holding the number of chars read.
uint *list_file_chunk::init_list() {
  uint *lst = new uint[mcl_t_size];
  assert(lst);
  memset(lst, 0, sizeof(uint)*mcl_t_size);
#ifndef NDEBUG
  for(uint i = 0; i < mcl_t_size; i++) {
    assert(lst[i] == 0);
  }
#endif
  return lst;
}

//! Inserts an element into the list given.
void list_file_chunk::insert_size(uint index, uint value) {
  if(!list) list = init_list();
  assert(list);
  assert(index < mcl_t_size);
  list[index] += value;
}

//! Resets the value of an element into the list given.
void list_file_chunk::reset_size(uint index, uint value) {
  if(!list) list = init_list();
  assert(list);
  assert(index < mcl_t_size);
  //! Resets:
  list[index] = value;
}

/**
   @brief Merges two lists.
   @param <src_obj>  The source of the data. Leaves it un-tuched.
   @remarks Does not copy data if src not set.
**/
void list_file_chunk::merge_lists(list_file_chunk *src_obj) {    
  if(src_obj && src_obj->get_list()) {
    uint *src = src_obj->get_list();
    //! Initiates if needed.
    if(!list) list = init_list();
    assert(list);
    //! Copies:
    for(uint i = 0; i < mcl_t_size; i++) {
      list[i] += src[i];
    }
#ifndef NDEBUG
    //! A pro-forma validation rule:
    for(uint i = 0; i < mcl_t_size; i++) {
      assert(list[i] >= src[i]);
    }
#endif
  }
}

//! @return the expected sice of the mci header for an mci file.
uint list_file_chunk::get_expected_size_of_header(taxa_t *listTaxa, int taxon_length) {
  assert(listTaxa);
  assert(taxon_length > 0);
  const uint size_proteins =  listTaxa[taxon_length-1].rel_end;
  uint cnt_chars_inserted = result_format::get_length_of_mcl_header(size_proteins);
  for(uint index_in = 0;index_in < size_proteins;index_in++) {
    if(index_in > 0) cnt_chars_inserted += (uint)log10(index_in)+1; // '+1' due to the logartihms attributes.
    else cnt_chars_inserted += 1; // a '0' only reserves one location in memory.
    cnt_chars_inserted++; // Making room for the space.
  }

  //! Adds "$\n)\n(mclmatrix\nbegin\n":
  cnt_chars_inserted += strlen(result_format::get_mcl_header_format_before_main());
  return cnt_chars_inserted;
}

/**
   @brief appends the header sizes to the list.
   @return a checksum corresponding to the total number of chars inserted
**/
uint list_file_chunk::append_header_sizes(taxa_t *listTaxa, int taxon_length) {
  const uint cnt_chars_inserted = list_file_chunk::get_expected_size_of_header(listTaxa, taxon_length);
    
  //! Prepares adding data to the list:
  if(!list) list = init_list();
  assert(list);
    
  //! Sets the data only in the "*.mci" files:
  for(uint i = 0; i < mcl_list_size; i++) {
    //! They are expected bein empty
    assert(!list[list_mcl_number[i]]);
    list[list_mcl_number[i]] = cnt_chars_inserted;
  }
  // The checksum: All files has the same attribute:
  const uint total_length = mcl_list_size * cnt_chars_inserted;
  return total_length;
}

//! @brief sets the size consumption for the header, given the settings:
void list_file_chunk::get_size_of_header(char *label, uint world_index, uint &size_string, uint &size_number) {
  size_string = 0, size_number = 0;
  if(PRINT_IN_ABC_FORMAT && !MODE_PAIRWISE_OUTPUT_ABC) size_string = mcl_bunch::get_size_of_header(label);
  if(PRINT_IN_MCL_FORMAT && !MODE_PAIRWISE_OUTPUT_MCL) {
    size_number = mcl_bunch::get_size_of_header(world_index);
  } 
}


//! Sets the values for the header using properties of this object.
void list_file_chunk::insert_sizes_header(uint world_index, char *name, uint type_string, uint type_number) {
  //! Gets the size:  
  uint size_string = 0, size_number = 0;
  get_size_of_header(name, world_index, size_string, size_number);
  //! Inserts the size:
  insert_size(type_string, size_string);
  insert_size(type_number, size_number);
}

/**
   @return the length of pair to be inserted.
   @remarks
   - For simplification assumes that the standard format is used:
   - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
   - Makes is possible using more effective writes when MPI is activated.
**/
void list_file_chunk::get_size_of_inserted_pair(uint world_index_out, char *label_in, char *label_out, float sim_score, uint &size_string, uint &size_number) {
  assert(MODE_PAIRWISE_OUTPUT_MCL == false);
  assert(MODE_PAIRWISE_OUTPUT_ABC == true);
  assert(PRINT_IN_ABC_FORMAT == true);
  assert(PRINT_IN_MCL_FORMAT == true);
  size_string = mcl_bunch::get_size_of_inserted_pair(label_in, label_out, sim_score);
  size_number = mcl_bunch::get_size_of_inserted_pair(world_index_out, sim_score);
}

/**
   @return the length of pair to be inserted.
   @remarks
   - For simplification assumes that the standard format is used:
   - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
   - Makes is possible using more effective writes when MPI is activated.
**/
void list_file_chunk::get_size_of_inserted_pair(uint world_index_out, char *label_out, float sim_score, uint &size_string, uint &size_number) {
  assert(MODE_PAIRWISE_OUTPUT_MCL == false);
  assert(MODE_PAIRWISE_OUTPUT_ABC == true);
  assert(PRINT_IN_ABC_FORMAT == true);
  assert(PRINT_IN_MCL_FORMAT == true);
  size_string = mcl_bunch::get_size_of_inserted_pair(label_out, sim_score);
  size_number = mcl_bunch::get_size_of_inserted_pair(world_index_out, sim_score);
}

//! Sets the values for the header using properties of this object.
void list_file_chunk::insert_sizes_for_a_pair(uint world_index_out, char *label_in, char *label_out, float sim_score, uint type_string, uint type_number) {
  uint size_string = 0, size_number = 0;
  get_size_of_inserted_pair(world_index_out, label_in, label_out, sim_score, size_string, size_number);
  insert_size(type_string, size_string);
  insert_size(type_number, size_number);
}

//! Sets the values for the header using properties of this object.
void list_file_chunk::insert_sizes_for_a_pair(uint world_index_out, char *label_out, float sim_score, uint type_string, uint type_number) {
  uint size_string = 0, size_number = 0;
  get_size_of_inserted_pair(world_index_out, label_out, sim_score, size_string, size_number);
  insert_size(type_string, size_string);
  insert_size(type_number, size_number);
}

//! De-allocates the memory reserved.
void list_file_chunk::free_memory() {
  if(list) {delete [] list; list=NULL;}
}

//! The constructor:
list_file_chunk::list_file_chunk() :
  MODE_PAIRWISE_OUTPUT_ABC(true), MODE_PAIRWISE_OUTPUT_MCL(false),
  PRINT_IN_ABC_FORMAT(true), PRINT_IN_MCL_FORMAT(true),
  list(NULL)
{};

//! The constructor:
list_file_chunk::list_file_chunk(
				 bool _MODE_PAIRWISE_OUTPUT_ABC,  // for the abc file: if set, the data out pairwise in stead of as in a row
				 bool _MODE_PAIRWISE_OUTPUT_MCL,  // for the mcl file: if set, the data out pairwise in stead of as in a row
				 // The following variables decides what output is to be printed
				 bool _PRINT_IN_ABC_FORMAT, 
				 bool _PRINT_IN_MCL_FORMAT,
				 bool _SORT_ABC_DATA
				 )
  :
  MODE_PAIRWISE_OUTPUT_ABC(_MODE_PAIRWISE_OUTPUT_ABC),  // for the abc file: if set), the data out pairwise in stead of as in a row
  MODE_PAIRWISE_OUTPUT_MCL(_MODE_PAIRWISE_OUTPUT_MCL),  // for the mcl file: if set), the data out pairwise in stead of as in a row
  // The following variables decides what output is to be printed
  PRINT_IN_ABC_FORMAT(_PRINT_IN_ABC_FORMAT), 
  PRINT_IN_MCL_FORMAT(_PRINT_IN_MCL_FORMAT),
  SORT_ABC_DATA(_SORT_ABC_DATA),
  list(NULL)
{};
