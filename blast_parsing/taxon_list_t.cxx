#include "taxon_list_t.h" 
#ifdef USE_MPI
#include "taxa_list_of_lists_settings.h"
#endif

//! Copies protein name from 'src' to 'dest', updating the 'dest_size' variable.
void taxon_list::copy_name(char *&dest, loint &size_dest, char *src, const loint size_src) {
  if(src) {
    if(dest) {
      if(size_dest < size_src) {
	delete [] dest; dest = new char[size_src+1]; dest[size_src] = '\0';	
      }
    } else {dest = new char[size_src+1]; dest[size_src] = '\0';}
    strncpy(dest, src, size_src);
  }
}

//! Copies taxon name from 'src' to 'dest'.
void taxon_list::copy_taxon_name(char *&dest, char *src) {    
  if(src) {
    const loint size_src = strlen(src);
    loint size_dest = 0; if(dest) size_dest = strlen(dest);
    copy_name(dest, size_dest, src, size_src);      
  }
}

/**
   @brief Inserts the name in the taxon_data object list *buffTaxon
   @param <taxon> The label of the taxon.
   @param <name> The label of the protein.
   @return True if a new taxon is inserted.
**/
bool taxon_list::insert_protein(char *taxon, char *name) {
  assert(taxon);
  assert(name);
  if(!listTaxon) listTaxon = taxon_data::init(10); // The default number of taxa to expect
  const int length = min((taxon_used+1),taxon_reserved);
  for(int i = 0; i< length; i++) {   // Iterates over the taxon used:
    if(listTaxon[i].is_equal(taxon)) {
      char *name_prot = listTaxon[i].insert_protein(name);
      assert(taxon_used > -1);
      // Not a new taxon; the protein is inserted, therefore returns
      if(name_prot != NULL) {
	const loint name_prot_size = strlen(name_prot);
	// The "last_name" is first updated when "first_name" has changed, due to the properties of the blastp-file used as input.
	if(protein_name_first_read && (*protein_name_first_read != '\0')) { // The protein name is set.
	  listTaxon[taxon_used].get_taxon_name(taxon_name_last_read, taxon_name_last_read_size);
	  // Updates the protein:
	  copy_name(protein_name_last_read, protein_name_last_read_size, name_prot, name_prot_size);
	} else {
	  listTaxon[taxon_used].get_taxon_name(taxon_name_first_read, taxon_name_first_read_size);
	  // Updates the protein:
	  copy_name(protein_name_first_read, protein_name_first_read_size, name_prot, name_prot_size);
	}
      }
      return (name_prot != NULL); // if the protein is set, it's found.
    } 
  }
  insert_taxon(taxon);   // If code reaches this point, the taxon is not found, and thereby a new taxon:  
  assert(taxon_used > -1);
  char *name_prot = listTaxon[taxon_used].insert_protein(name); // the protein is added to the chain.
  // Updates the data whom in the merging to be used, in order deciding if two names are equal, and therefore must be merged into one name
  if(name_prot!= NULL) { // If the name is found:
    const loint name_prot_size = strlen(name_prot);
    if(protein_name_first_read && (*taxon_name_first_read != '\0')) {
      listTaxon[taxon_used].get_taxon_name(taxon_name_last_read, taxon_name_last_read_size);
      // Updates the protein:
      copy_name(protein_name_last_read, protein_name_last_read_size, name_prot, name_prot_size);
    } else {
      listTaxon[taxon_used].get_taxon_name(taxon_name_first_read, taxon_name_first_read_size);
      // Updates the protein:
      copy_name(protein_name_first_read, protein_name_first_read_size, name_prot, name_prot_size);
    }
  }
  return (name_prot != NULL);
}

//! Inserts a new taxon into the list
uint taxon_list::insert_taxon(char *taxon) {
  if(taxon) { // Only inserts if it's set.
    taxon_used++;  enlarge_taxon_list();
    listTaxon[taxon_used].set_taxon_name(taxon);
  }
  return taxon_used;
}

/**
   @brief Enlarges the taxon list if not big enough
   @date 14.01.2011 by oekseth
*/
void taxon_list::enlarge_taxon_list() {
  if(taxon_used >= taxon_reserved) { // As a safety mechanism
    int taxon_reserved_new = 2* taxon_reserved;
    if(!taxon_reserved_new) taxon_reserved_new = 10;
    //    int taxon_reserved_new = 10+taxon_reserved;

    // First: allocates a new list:
    taxon_data_t *listTaxon_new = taxon_data::init(taxon_reserved_new);
    assert(listTaxon_new); // If memory not avaiblable, this will fail, but should not.
    // Second: 'Copies' the data by moving the pointers from the old object- ot the new one:
    if(listTaxon) { // No point in copying if it's empy
      for(int i = 0; i < taxon_reserved; i++) {
	listTaxon[i].take_and_reset_data(listTaxon_new[i]); // After this operation it's safe deallocationg the listTaxon object.
      }
      // Deallocates the lsit, but not the data itself (as we have used them in a different one):
      taxon_data::close(listTaxon);
    }
    listTaxon = listTaxon_new; // Updates the pointers' pointer
    taxon_reserved = taxon_reserved_new;
  }
}

//! If argument set to true, prints the proteins as well
void taxon_list::printf_buffer(const bool print_prots)
{
  if(listTaxon) {
    printf("The buffer contains the following taxon:\n");
    const int length = min((taxon_used+1),taxon_reserved);
    for(int i = 0; i< length; i++) {
      printf("[%d]\t%s\t%u proteins", i, listTaxon[i].get_taxon_name(), listTaxon[i].get_prots_used()+1);
      if(print_prots) printf(";\t"), listTaxon[i].printf_buffer();
      printf("\n--------------------------\n");
    }
  } else printf("(the taxon_data buffer is empty.)\n");
}

/**
   @brief Produces (writes) a log file holding the list of memory allocations.
   @remarks 
   - Useful for optimizing the memory allocation procedures for building of lists for different taxa.
   - De-activated in in NDEBUG mode.
**/
void taxon_list::log_produce_memory_allocations() {
#ifndef NDEBUG
#ifdef LOG_WRITE_MEMORY_ALLOCATIONS_DURING_BLASTP_PARSING
  struct tms tmsstart; clock_t clock_log;
  assert(listTaxon);
  if((clock_log = times(&tmsstart)) == -1) // starting values
    fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
  char sep = '_';
  // The log file, and contents of it:
  if(sizeof(FUNC_LOG_FOLDER_NAME) > 0) sep = '/';
  char str[200]; memset(str, '\0', 200);
  sprintf(str, "%s%c%s.log", FUNC_LOG_FOLDER_NAME, sep, "taxa_allocations");
  FILE *f_log = fopen(str, "w");
  if(f_log) {
    uint prots_reserved = 0, prots_used = 0;
    const int length = min((taxon_used+1),taxon_reserved);
    for(int i = 0; i< length; i++) {
      prots_used+=listTaxon[i].get_prots_used(), prots_reserved+=listTaxon[i].prots_reserved;
    }
    fprintf(f_log, "The list have a total fill-factor of %f, and a total memory consumption of %fMB:\n", (float)prots_used/(float)prots_reserved, (float)sizeof(taxon_data)*(float)prots_reserved/(1024*1024));
    for(int i = 0; i< length; i++) {
      float fill_factor = 1;
      if(listTaxon[i].prots_reserved) {
	fill_factor = (float)(listTaxon[i].get_prots_used()+1)/(float)listTaxon[i].prots_reserved;
      }
      fprintf(f_log, "[%d]\t%s\t with a fill-factor of %f (%u proteins used, %u proteins reserved)\n",
	      i, listTaxon[i].get_taxon_name(),
	      fill_factor,
	      listTaxon[i].get_prots_used()+1,
	      listTaxon[i].get_prots_reserved());
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
#endif
}
//! Returns the number of taxa in teh collection
uint taxon_list::getLength() {
  const int length = min((taxon_used+1),taxon_reserved);
  return (uint)length;
}

//! Returns the number of proteins for the given taxon
uint taxon_list::getProteinLength(uint i) {
  if(listTaxon && ((int)i < taxon_reserved)) {
    return listTaxon[i].prots_used;
  } else return 0;
}

//! Returns the given buffer from 'lisTaxon'
char **taxon_list::getBuffer(uint i) {
  if(listTaxon && ((int)i < taxon_reserved)) {
    return listTaxon[i].get_buffer();
  } else return NULL;
}

//! Returns the given name from listTaxon    
char *taxon_list::getName(uint i) {
  assert(listTaxon);
  if((int)i < taxon_reserved) {
    //    printf("henter taxon(%u) (totalt brukt %u+1 taxa)\n", i, taxon_used);
    return listTaxon[i].get_taxon_name();
  } else {
    fprintf(stderr, "!!\tInternal taxon index %u outside boundary of %d\n", i, taxon_reserved);
    return NULL;
  }
}


/**! Merges two arrays.
   @Precondition: Must be in the same order
   @Procedure:
   1. Iterate through 'this':
   a. If there exists a 'taxon_name' equal to an element in 'this',
   - call 'mergeData'
   - deallocate the memory for the arg' taxon
   /  Else, junp on next element into 1a)
   2. Insert taxa not previously known:
   a. If there exists an element not found in 'this':
   -  enlarge, if necessary, the taxon list
   -  Copy the data pointer to the new 'taxon_list'element
   3. Update the variables
   - 'taxon_name_last_read'
   - 'protein_name_last_read'
*/
void taxon_list::merge(taxon_list *arg) {
  //    1. Iterate through 'this':
  if(arg) {
    if(listTaxon) { // Only a point doing this if it has data
      const int length = min((taxon_used+1),taxon_reserved);
      for(int i = 0; i< length; i++) {
	taxon_data buff = taxon_data(); // An empty object.
	if(arg->hasTaxon(listTaxon[i].get_taxon_name(), buff)) { // iterate thwough all of 'args' possibilities
	  const bool has_overlap = is_overlap(arg[0], listTaxon[i].get_taxon_name());
	  listTaxon[i].mergeData(buff, has_overlap);
	} 
      }
    }
    //  2. Insert taxa not previously known
    for(int i = 0; i< (int)arg->getLength(); i++) {
      taxon_data buff = taxon_data();
      if(!hasTaxon(arg->listTaxon[i].get_taxon_name(), buff)) { // those not found added
	const uint taxon_id = insert_taxon(arg->listTaxon[i].get_taxon_name()); // Gets an index and allocates memory for the list.
	listTaxon[taxon_id].mergeData(arg->listTaxon[i], false);
      } 
    }

    //      3. Update the variables
    // The taxon label:
    copy_name(taxon_name_last_read, taxon_name_last_read_size, arg->get_taxon_name_last_read(), arg->get_taxon_name_last_read_length());
    // The protein label:
    copy_name(protein_name_last_read, protein_name_last_read_size, arg->get_protein_name_last_read(), arg->get_protein_name_last_read_length());
    arg->free_mem();
  }
}

/**
   @brief Merges an element into the list
   @param <arg> The input to merge.
   @remarks Assumes there is no overlap between proteins due to disjoint regions.
*/
void taxon_list::merge(taxon_data &arg) {
  //!   1. Iterate through 'this':
  if(arg.has_data()) {
    bool taxon_inserted = false;
    if(listTaxon) { // Only a point doing this if it has data
      const int length = min((taxon_used+1),taxon_reserved);
      for(int i = 0; i< length; i++) {
	if(arg.is_taxon(listTaxon[i].get_taxon_name())) {
	  listTaxon[i].mergeData(arg, false); // The overlap should be false
	  taxon_inserted = true;
	} 
      }
    }
    if(!taxon_inserted) { // Taxon were not found in 'this' list; inserts.
      //!  2. Insert taxa not previously known
      taxon_data buff = taxon_data();
      assert(!hasTaxon(arg.get_taxon_name(), buff));  // those not found added
      const uint taxon_id = insert_taxon(arg.get_taxon_name()); // Gets an index and allocates memory for the list.
      listTaxon[taxon_id].mergeData(arg, false);
    }
    arg.free_mem();
  }
}


/**! Returns true if the taxon given as argument is the same as the first taxon in both the argument and 'this', and in addition that the proteins are equal, and therefore the last shall be skipped
 */
bool taxon_list::is_overlap(taxon_list arg, char *taxon) {
  if(taxon) {
    if(arg.is_equal_taxon_name_first(taxon)) {
      // Its in the current taxon
      if(arg.is_equal_taxon_name_first(taxon_name_last_read)) {
	if(arg.is_equal_protein_name_first(protein_name_last_read)) 
	  return true; // thye have an overlap;
      }
    }
  }
  return false;
}


//! Returns true if this list has the given taxon set as argument
bool taxon_list::hasTaxon(char *taxon_name, taxon_data &buff) {
  if(listTaxon && taxon_name) { // Requires both some data, and input-argument set.
    const int length = min((taxon_used+1),taxon_reserved);
    for(int i = 0; i< length; i++) {
      if(listTaxon[i].is_taxon(taxon_name)) {
	buff = listTaxon[i];
	return true;
      }
    }
  }
  return false;
}

/**
   @return The taxon index for the given param.
   @remark Used by getProteinLength(char *taxon)
**/
bool taxon_list::getTaxonIndex(char *taxon_name, int &taxon_name_index) {
  if(listTaxon && taxon_name) {
    const int length = min((taxon_used+1),taxon_reserved);
    for(int i = 0; i< length; i++) {
      if(listTaxon[i].is_taxon(taxon_name)) {
	taxon_name_index = i;
	return true;
      }
    }
  } 
  taxon_name_index = INT_MAX;
  return false;
}

/**
   @return The protein length of the taxon given as argument
   @remark Assumes that the taxon exsits: If not, returns the 0-index(0)
**/
bool taxon_list::getProteinLength_n(char *taxon, int &taxon_name_index) {
  if(listTaxon && taxon) {
    int index = 0;
    const bool has_found_index = getTaxonIndex(taxon, index);
    if(has_found_index) {
      taxon_name_index = getProteinLength(index);
      return true;
    } else return false;
  } else return false;
} 

/**
   @brief Speeds up the exectuion for files having a bunch of taxa with close to empty collections.
   @param <threshold> The number of proteins a taxon must have more of in order to be evaluated.
   @param <taxon_length> The variable to update.
**/
void taxon_list::remove_taxa_with_proteins_below_threshold(uint threshold, int &taxon_length) {
  if(threshold) { 
    const int length = min((taxon_used+1),taxon_reserved);
    taxon_data *updated_list = taxon_data::init(length);
    log_builder::test_memory_condition_and_if_not_abort((updated_list!= NULL), __LINE__, __FILE__, __FUNCTION__);
    int updated_list_cnt = -1;
    if(listTaxon) { // Only a point doing this if it has data
      for(int i = 0; i< length; i++) {
	if(listTaxon[i].get_length() > threshold) { // iterate thwough all of 'args' possibilities
	  updated_list[++updated_list_cnt].steal_object(listTaxon[i]);
	} else {
	  listTaxon[i].delete_buffer();
	}
      }
    }
    if((updated_list_cnt>-1) && (updated_list_cnt < taxon_used)) {
      const int length = min((taxon_used+1),taxon_reserved);
      for(int i = updated_list_cnt+1; i< length; i++) {
	updated_list[i] =taxon_data();
      }
    }
    if(updated_list_cnt>-1) {
      //      printf("--\tInserted %d elements of a total length(%d) given a threshold(%u) at line %d in file %s\n", updated_list_cnt, taxon_used+1, threshold, __LINE__, __FILE__); 
      taxon_data::close(listTaxon);
      taxon_used = updated_list_cnt;
      taxon_reserved = length;
      taxon_length = getLength();
      listTaxon =  updated_list;
      assert(listTaxon);
#ifndef NDEBUG
      for(int taxon_id = 0; taxon_id<taxon_length;taxon_id++) {
	assert(updated_list[taxon_id].get_buffer()); 
	if(!(getBuffer(taxon_id))) {
	  //	  printf("[%d]\tbuffer not set\n", taxon_id);
	  assert(getBuffer(taxon_id)!=NULL);
	}
      }
#endif
    } else {
      taxon_data::close(listTaxon);
      taxon_length = 0;
      taxon_used = -1;
      taxon_reserved = 0;
      listTaxon =  NULL;
      const char *error = "The threshold for the number of proteins a taxon must have, set at terminal input, caused every taxon-pair to be discarded, therefore no output is produced";
      log_builder::throw_error_and_abort(error, __LINE__, __FILE__, __FUNCTION__);
    }
  }
}

//!  @return True if the strings given are equal.
bool taxon_list::compare_strings(char *name_1, char *name_2, uint length_string) {
  if(name_1 && name_2) {
    return (0 == strncmp(name_1, name_2, length_string));
  } else { // If both are not set, they are equal.
    if((name_1 == NULL) && (name_2 == NULL)) return true;
    else return false; // only one of them is set ie, they are not equal.
  }
}

//!  @return True if the strings given are equal.
bool taxon_list::compare_strings(char *name_1, char *name_2) {
  if(name_1 && name_2) {
    return (0 == strcmp(name_1, name_2));
  } else { // If both are not set, they are equal.
    if((name_1 == NULL) && (name_2 == NULL)) return true;
    else return false; // only one of them is set ie, they are not equal.
  }
}

/**
   @remarks 
   - Used to verify that the buffers are equal
   - Uses pointers, meaning extra vlidation testing they are set is included.
   @return true if the input buffers are equal on all of their elements.
**/
bool taxon_list::compare_buffers(char **buff_1, char **buff_2, uint length_buffer, uint length_string) {
  if(buff_1 && buff_2) {
    for(uint i = 0; i < length_buffer; i++) {
      if(buff_1[i] && buff_2[i]) {
	if(0 != strncmp(buff_1[i], buff_2[i], length_string)) return false; 
      } else {
	if(buff_1[i] || buff_2[i]) return false; // only one of them is set ie, they are not equal.
      }
    }
    return true;
  } else {
    if(!buff_1 && !buff_2) return true; // If both are not set, they are equal.
    else return false; // only one of them is set ie, they are not equal.
  }
}


#ifdef USE_MPI

//! Initiates an object of type taxa_buffer_list_settings
void build_mpi_obj_of_taxa_buffer_list_settings(taxa_buffer_list_settings_t t_settings, MPI_Comm comm, int myrank, MPI_Datatype &mpi_taxa_buffer_list_settings) {
  //! One value of each type:
  int block_lengths[4] = {1,1,1,1};
  //! The base types:
  MPI_Datatype old_types[4] = {MPI_UNSIGNED,MPI_UNSIGNED,MPI_UNSIGNED,MPI_UNSIGNED};
  //! The location of each element;
  static const int length = 4;
  MPI_Aint indices[length];
  MPI_Get_address(&t_settings.taxon_size, &indices[0]);
  MPI_Get_address(&t_settings.taxonNamesArr_size, &indices[1]);
  MPI_Get_address(&t_settings.total_cnt_proteins, &indices[2]);
  MPI_Get_address(&t_settings.total_cnt_chars_used_representing_proteins, &indices[3]);
  //! Makes the locations relative:
  for(int i = length-1; i > 0; i--) {
    indices[i] = indices[i] - indices[0];
  }
  indices[0] = 0;
  //! Builds the mpi representation of the structure:
  MPI_Type_struct( 4, block_lengths, indices, old_types, &mpi_taxa_buffer_list_settings);
  MPI_Type_commit(&mpi_taxa_buffer_list_settings);
}


//! Builds an mpi object, including it into the two structures given as parameters.
void taxon_list::build_mpi_1d_lists(taxa_buffer_list_settings_t &obj_lists_length, char *&taxonNamesArr, uint *&taxon_prots_used, char *&proteinLabels, uint *&bufferProteinsLength, int myrank, int number_of_nodes) {
    uint taxonNamesArr_size = 0; 
    const uint taxon_size  = getLength();
    taxonNamesArr = NULL; // For storing all the taxon-names in the same 1-d-list, accepting variable-sizes of the char-strings
    uint taxonNamesArr_used = 0;
    taxon_prots_used = NULL; // The number of proteins for each taxa.
    proteinLabels = NULL;// The buffer of protein labels
    uint proteinLabels_char_used = 0; // The number of chars used in total.
    uint bufferProteinsLength_size = 0;   // The size of the list of protein label-sizes.
    bufferProteinsLength = NULL; // The list of protein label-sizes.
    //    if(myrank) {
      //! Taxon-names: Below variables used storing all the taxon-names in the same 1-d-list, accepting variable-sizes of the char-strings
      taxonNamesArr_size = taxon_size*40; 
      taxonNamesArr = new char[taxonNamesArr_size];
      memset(taxonNamesArr, '\0', taxonNamesArr_size);
      //! General data about the taxa:
      taxon_prots_used = new uint[taxon_size];
      for(uint i = 0; i < taxon_size; i++) {
	taxon_prots_used[i]         = listTaxon[i].get_prots_used();
	bufferProteinsLength_size  += listTaxon[i].get_prots_used();
      }
      //! Protein-names: Below variables used storing a 2d-char-list in the same 1-d list, accepting variable-sizes of the char-strings.
      bufferProteinsLength = new uint[bufferProteinsLength_size]; // The list of protein label-sizes.
      memset(bufferProteinsLength, 0, sizeof(uint)*bufferProteinsLength_size);  
      uint proteinLabels_size = 0;  // An upper estimate of the size of the buffers of protein labels

      proteinLabels_size = 40*bufferProteinsLength_size;  // An upper estimate of the size of the buffers of protein labels
      assert(proteinLabels_size);
      proteinLabels = new char[proteinLabels_size];  // The buffer of protein labels
      memset(proteinLabels, '\0', proteinLabels_size);
      proteinLabels_char_used = 0; // The number of chars used in total.
      uint bufferLength_startpos = 0; // The number of chars used for the taxa before the current taxon in working.

      for(uint i = 0; i < taxon_size; i++) {
	//! Taxa-labels: the memcpy-operation:
	const uint start_name = taxonNamesArr_used;
	listTaxon[i].insert_the_name_into_argument(taxonNamesArr, taxonNamesArr_size, taxonNamesArr_used);	
	listTaxon[i].transform_protein_label_buffer_into_1d_list(proteinLabels, proteinLabels_char_used, proteinLabels_size, bufferProteinsLength, bufferProteinsLength_size, bufferLength_startpos);
#ifndef NDEBUG
	//! Taxon-labels: Validate the new-enlarged list:
	assert(0==strcmp(listTaxon[i].get_taxon_name(), taxonNamesArr + start_name));
	//! Protein-labels: Validate the new-built list:
	uint index_startpos = 0, buffer_startpos = 0;
	for(uint out = 0; out < i+1; out++) {
	  listTaxon[out].assert_1d_list(index_startpos, buffer_startpos, proteinLabels, bufferProteinsLength);
	  index_startpos += taxon_prots_used[out];
	}
#endif
	bufferLength_startpos += taxon_prots_used[i];
      }
#ifndef NDEBUG
      //! Validate the new-built list:
      uint index_startpos = 0, buffer_startpos = 0, taxon_name_index = 0;
      for(uint i = 0; i < taxon_size; i++) {
	//! Taxon-labels: Validate the list:
	assert(0==strcmp(listTaxon[i].get_taxon_name(), taxonNamesArr + taxon_name_index));

	//! Protein-labels: Validate the list: Note: 'buffer_startpos' is update by the function called.
	const uint buffer_startpos_pre = buffer_startpos;
	listTaxon[i].assert_1d_list(index_startpos, buffer_startpos, proteinLabels, bufferProteinsLength);

	//! Validate the transformation of this new-built lists into the original type (ie, a set of taxon_data objects).
	taxon_data obj = taxon_data(taxon_prots_used[i]);
	char *listOfProteins = proteinLabels    + buffer_startpos_pre;
	const uint protein_length = taxon_prots_used[i];
	uint *listOfLengths  = bufferProteinsLength + index_startpos;
	obj.set_labels_using_1d_object(protein_length, listOfProteins, listOfLengths);
	obj.setTaxonName(taxonNamesArr + taxon_name_index);
	assert(listTaxon[i].compare_object(obj, true));
	obj.free_mem();
	
	//! Update the variables for the next run:
	taxon_name_index += strlen(taxonNamesArr + taxon_name_index) +1; // Taxon-name-variables
	index_startpos += taxon_prots_used[i]; // Protein-label-variables.
      }      
#endif
//     } else {
//       obj_lists_length   = taxa_buffer_list_settings_t();
//     }
    obj_lists_length.taxon_size = taxon_size;
    obj_lists_length.taxonNamesArr_size = taxonNamesArr_used;
    obj_lists_length.total_cnt_proteins = bufferProteinsLength_size;
    obj_lists_length.total_cnt_chars_used_representing_proteins = proteinLabels_char_used;
    
}
/**
   @brief Sends- and receives the list of taxa- and protein labels accross nodes.
**/
void taxon_list::mpi_make_data_consistent_accross_nodes(MPI_Comm comm, int myrank, int number_of_nodes, const uint LIMIT_MINIMUM_NUMBER_OF_PROTEINS_FOR_EACH_TAXA, int &taxon_length_to_be_updated) {
   uint taxon_size  = getLength();
   if(listTaxon && taxon_size>0) {
    
     MPI_Datatype mpi_taxa_buffer_list_settings; // The type to define the settings about the taxa
     //    MPI_Datatype mpi_list_of_lists; // The type to define the collection of lists
     //! Collect the proteins.
     //    if(true) { // An simple approach, letting myrank=0 get all
     taxa_buffer_list_settings_t t_settings;// The object to hold the settings about the taxa
     build_mpi_obj_of_taxa_buffer_list_settings(t_settings, comm, myrank, mpi_taxa_buffer_list_settings);

     //      build_mpi_obj_of_taxa_buffer_list_settings(t_list_of_lists, comm, myrank, mpi_list_of_lists);
     enum send_types {object_info, object_list_of_lists};
     if(myrank) {
       //! Builds the 1-d-lists
       char *taxonNamesArr; // For storing all the taxon-names in the same 1-d-list, accepting variable-sizes of the char-strings
       uint *taxon_prots_used; // The number of proteins for each taxa.
       char *proteinLabels;// The buffer of protein labels
       uint *bufferProteinsLength; // The list of protein label-sizes.

       //! Transform the list of taxa into 1d, and sets the varaibles:
       build_mpi_1d_lists(t_settings,taxonNamesArr, taxon_prots_used, proteinLabels, bufferProteinsLength, myrank, number_of_nodes);
       //! The object to hold the settings about the taxa:
       taxa_list_of_lists_settings_t list_of_lists = taxa_list_of_lists_settings_t(taxonNamesArr, taxon_prots_used, proteinLabels, bufferProteinsLength);
       
       //! Sends the data of the sizes (meta-data):
       int ret_val=0; ret_val = MPI_Send(&t_settings, 1, mpi_taxa_buffer_list_settings, 0, object_info, comm);
	
       //! Sends the lists holding the real data:
       list_of_lists.send_data(comm, myrank, object_list_of_lists, t_settings);
       list_of_lists.free_memory();
     } else {
       //! Receives the length of the lists:
       for(int sender_id = 1; sender_id < (int)number_of_nodes; sender_id++) {
	 taxa_list_of_lists_settings_t list_of_lists;// The object to hold the settings about the taxa
	 //! Gets the length of the lists:
	 t_settings = taxa_buffer_list_settings_t();
	 MPI_Status status_object_info;
	 int ret_val=0; ret_val = MPI_Recv(&t_settings, 1, mpi_taxa_buffer_list_settings, sender_id, object_info, comm, &status_object_info) ;
	 int cnt_received = 0; //const int cnt_received_ERR = ;
	 MPI_Get_count(&status_object_info, mpi_taxa_buffer_list_settings, &cnt_received);

	 //! Receives the lists, and merges them:
	 list_of_lists = taxa_list_of_lists_settings_t(t_settings); // Allocates memory
	 list_of_lists.receive_and_merge_data(this, comm, myrank, sender_id, object_list_of_lists, t_settings);
	 list_of_lists.free_memory();
       }
     }
      
     //    }
     //! Do some preprocessing
     if(myrank != 0) {
       free_memory(); // TODO: Consider if this as a waste of resources!
       //! Gets the length of the lists:
       t_settings = taxa_buffer_list_settings_t();
       //       MPI_Status status_object_info;
       const int root_id = 0; // The id of the node delegated the root responsibility
       //       int ret_val=0;
       // Receives data from the broadcast:
       MPI_Bcast(&t_settings, 1, mpi_taxa_buffer_list_settings, root_id, comm);
       
       //! Receives the lists, and merges them:
       taxa_list_of_lists_settings_t list_of_lists = taxa_list_of_lists_settings_t(t_settings); // Allocates memory

       list_of_lists.receive_and_merge_data(this, comm, myrank, root_id, object_list_of_lists, t_settings);
       //	  int cnt_received = 0; const int cnt_received_ERR = MPI_Get_count(&status_object_info, mpi_taxa_buffer_list_settings, &cnt_received);
       list_of_lists.free_memory();
     } else { // All responsibility delegated to the roote node of updating the others:
       //! Cleans unneccessay "stuff"
       remove_taxa_with_proteins_below_threshold(LIMIT_MINIMUM_NUMBER_OF_PROTEINS_FOR_EACH_TAXA, taxon_length_to_be_updated);      
       taxon_size  = getLength(); // Updates the internal variable.
      
       //! Builds the 1-d-lists
       char *taxonNamesArr; // For storing all the taxon-names in the same 1-d-list, accepting variable-sizes of the char-strings
       uint *taxon_prots_used; // The number of proteins for each taxa.
       char *proteinLabels;// The buffer of protein labels
       uint *bufferProteinsLength; // The list of protein label-sizes.
      
       //! Transform the list of taxa into 1d, and sets the varaibles:
       t_settings = taxa_buffer_list_settings_t();
       build_mpi_1d_lists(t_settings,taxonNamesArr, taxon_prots_used, proteinLabels, bufferProteinsLength, myrank, number_of_nodes);
       //! The object to hold the settings about the taxa:
       taxa_list_of_lists_settings_t list_of_lists = taxa_list_of_lists_settings_t(taxonNamesArr, taxon_prots_used, proteinLabels, bufferProteinsLength);

       //! Sends the data of the sizes (meta-data):       
       MPI_Bcast(&t_settings, 1, mpi_taxa_buffer_list_settings, myrank, comm);
       
       //! Sends the lists holding the real data to all the nodes:
       list_of_lists.broadcast_send_data(comm, myrank, number_of_nodes, object_list_of_lists, t_settings);
       list_of_lists.free_memory();
     }

    //! Cleans the objects used, reminv memory allocations:
    MPI_Type_free(&mpi_taxa_buffer_list_settings);  
    //    MPI_Type_free(&mpi_list_of_lists);
  }
}
#endif
//! Returns an object of this class.
taxon_list *taxon_list::init() {
  return new taxon_list();
}
//! Frees the memory
void taxon_list::free_mem() {   delete_buffer(); }

/**   @brief Deallocates the memory given the param.**/
void taxon_list::close(taxon_list *&obj) {
  if(obj){obj->free_mem(), delete obj, obj=NULL;}
}
//@argument: 'size' the size of the name-pointer to be (1) allocated, (2) Intialized to '\0'
char* taxon_list::init_name(uint size) {
  char *name = new char[size]; //(char*)malloc(sizeof(char)*size);
  memset(name, '\0', size);
  return name;
}

//! Deletes the pointers/memory references in this structure
void taxon_list::delete_buffer() {
  taxon_data::close(listTaxon, (uint)taxon_reserved);
  if(taxon_name_first_read) {  delete [] taxon_name_first_read; taxon_name_first_read = NULL;}
  if(taxon_name_last_read) {delete [] taxon_name_last_read; taxon_name_last_read = NULL;}
  if(protein_name_first_read) {delete [] protein_name_first_read; protein_name_first_read = NULL;}
  if(protein_name_last_read) {delete [] protein_name_last_read; protein_name_last_read = NULL;}
  taxon_used = -1; taxon_reserved = 0;
}


taxon_list::taxon_list() : listTaxon(NULL), taxon_name_first_read(NULL), taxon_name_first_read_size(0),
			   taxon_name_last_read(NULL), taxon_name_last_read_size(0),
			   protein_name_first_read(NULL), protein_name_first_read_size(0),
			   protein_name_last_read(NULL), protein_name_last_read_size(0),
			   taxon_used(-1), taxon_reserved(0) // The default size of the taxon
{
}

//#ifdef assert_code
//! Verifies the length set of the taxa and the prots
void taxon_list::assert_getLength() {
  taxon_list list = taxon_list();
  char *name[] = {"sara", "helmut", "jon", "lars"};
  uint name_cnt = 0;
  char *taxon[] = {"bergenser", "norsk", "rosenborg", "tiger"};
  uint taxon_cnt = 0;
  //test_cnt++; if(show_positive_results) printf("??\t Running test %u: Verifying the length of the taxons and the proteins\n", test_cnt);
  list.insert_protein(taxon[taxon_cnt++], name[name_cnt++]);
  assert((taxon_cnt == list.getLength()));
  int index_p = 0;
  bool found = list.getProteinLength_n(name[name_cnt-1], index_p);
  assert(name_cnt == (uint)index_p);
  assert(found);
  list.insert_protein(taxon[taxon_cnt++], name[name_cnt]);
  assert((taxon_cnt == list.getLength()));
  found = list.getProteinLength_n(name[name_cnt-1], index_p);
  assert(name_cnt == (uint)index_p);
  list.insert_protein(taxon[taxon_cnt], name[name_cnt++]);
  assert((taxon_cnt+1) == list.getLength()); // pluss one due to the fact that iut starts at '0', and note one
  found = list.getProteinLength_n(taxon[taxon_cnt], index_p);
  assert(1 == index_p);
  assert(found);
  list.delete_buffer();
}

void taxon_list::assert_insert_protein() {
#ifndef NDEBUG
  taxon_list list = taxon_list();
  char *taxon[] = {"bergenser", "norsk", "rosenborg", "tiger"}; 
  char *name[] = {"hanne", "per", "knut", "mor", "einar"}; 
  bool result = list.insert_protein(taxon[0], name[0]);
  assert(result); // verifies that the name was not inserted 
  result = list.insert_protein(taxon[0], name[0]);
  assert(!result); // verifies that the name was not inserted
  assert(1==list.getTaxaLength());
  loint size_taxon = 0;
  if(taxon[0]) size_taxon = strlen(taxon[0]);
  bool equal = compare_strings(list.getTaxonName(0), taxon[0], size_taxon);
  assert(equal);
  loint size_protein = strlen(name[0]);
  equal = compare_buffers(list.getBuffer(0), name, 1, size_protein);
  assert(equal);
  assert(1==list.getProteinLength(0));
  assert(1 == list.insert_taxon(taxon[1])); // the updated taxon index
  int index; const bool found = list.getTaxonIndex(taxon[1], index);
  assert(1 == index); // getting the taxon index
  assert(found);
  list.delete_buffer();
#endif
}

void taxon_list::assert_class(const bool print_info) {
  const static char *class_name = "taxon_list";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
#ifndef NDEBUG
  taxon_list list = taxon_list();
  list.assert_insert_protein();
  list.assert_getLength(); //! Verifies the length set of the taxa and the prots
  list.delete_buffer();
  
  list = taxon_list();
  const uint size_in = 11, size_out = 1;
  for(uint i = 0; i < size_in; i++) {
    for(uint k = 0; k < size_out; k++) {
      char in[20]; for(uint m = 0; m <20;m++) in[m] = '\0';
      char out[20]; for(uint m = 0; m <20;m++) out[m] = '\0';
      sprintf(in, "%u\n", i); 
      sprintf(out, "%u\n", k); 
      list.insert_protein(in, out);
    }
  }
  list.delete_buffer();
#endif
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
