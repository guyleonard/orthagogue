#include "taxa.h"

//! @returns the inparalog limit using the local index
float taxa::getArrInpaLimit(uint local_index) {
  assert(local_index < (uint)total_cnt);
  if(arrInpaLimit) return arrInpaLimit[local_index];
  else return 0;
}

/**! 
   @Name: mcl_get_all_header_proteins(..) 
   -- Produces the mcl header, having a very long string as the input
   @Changed: 07.01.2011 by oekseth (init)
   @Changed: 25.08.2011 by oekseth (validation and small changes)
   
*/
char *taxa::mcl_get_all_header_proteins(char *string) {
  if(name) {
    const uint name_length = strlen(name);
    for(uint i=0; i< (uint)total_cnt; i++) {
      int length =0;
      bool end = false;
      if(arrKey[i]) { // only work on it if data is set
	for(uint k = 0;k < strlen(arrKey[i]) && !end; k++) {
	  if(arrKey[i][k] != '\0') {
	    length = k+1;     
	  } else end = true;
	}
	memcpy(string, arrKey[i], length);
	string += length;
	*string = SEPERATOR;
	string++;
	strncpy(string, name, name_length);
	string += name_length;
	*string = MCL_HEAD_SEPERATOR; // Adding the seperator
	string++;
      }
    }
    return string;
  } else {string = NULL; return string;}
}

//! Inserts the overlap list without allocting memory
void taxa::set_unsafe_overlap_list(uint *&over) {
  assert(!arrOverlap);
  arrOverlap = over;
}

/**
   Copies the overlap list into this
   @remarks Assumes that if the lsit in this object is set, it's set tot the proper length.
**/
void taxa::copy_arr_overlap_into_list(uint *over) {
  if(!arrOverlap) arrOverlap = new uint[total_cnt];
  log_builder::test_memory_condition_and_if_not_abort(arrOverlap!=NULL, __LINE__, __FILE__, __FUNCTION__);
  memcpy(arrOverlap, over, sizeof(uint)*total_cnt);
}


/**
   Copies the inpaLimit list into this
   @remarks Assumes that if the list in this object is set, it's set tot the proper length.
**/
void taxa::copy_arr_inpalimit_into_list(float *inpalist) {
  if(!arrInpaLimit) arrInpaLimit = new float[total_cnt];
  log_builder::test_memory_condition_and_if_not_abort(arrInpaLimit!=NULL, __LINE__, __FILE__, __FUNCTION__);
  memcpy(arrInpaLimit, inpalist, sizeof(float)*total_cnt);
}


/**
   @return The total overlap score
   @remarks: Useful function for validating after an update
**/
uint taxa::get_total_overlap_score() {
  if(!total_cnt || !arrOverlap) return 0;
  uint sum = 0;
  for(uint i = 0; i < (uint)total_cnt; i++) {
    sum += arrOverlap[i];
  }
  return sum;
}

/**
   @return The total overlap score
   @remarks: Useful function for validating after an update
**/
uint taxa::get_total_overlap_score(taxa *&listTaxa, uint index_start, uint index_end_pluss_one) {
  assert(listTaxa);
  assert(index_end_pluss_one);
  uint sum = 0;
  for(uint i = index_start; i < index_end_pluss_one; i++) {
    sum += listTaxa[i].get_total_overlap_score();
  }
  return sum;
}

//! @return the upper estimate for the object size
uint taxa::get_object_size(uint taxon_length) {
  const uint length = rel_end-rel_start;
  uint size = sizeof(taxa);
  size += SIZE_TAXO + (length)*SIZE_PROTEIN; // 'name' + 'arrKey' using approximate size
  size += length*(taxon_length+1); // An approx for the number of ortho- and inparalog elements reserved.    
  return size;
}

//! @return an approximation for the memory consumption.
uint taxa::get_object_memory_size() {
  const uint length = total_cnt; //rel_end-rel_start;
  uint size = sizeof(taxa);
  size += SIZE_TAXO + (length)*SIZE_PROTEIN; // 'name' + 'arrKey' using approximate size
  size += sizeof(float)*length*2; // An approx for the number of overlap- and inparalog elements reserved.    
  return size;
}

/**
   @brief Prints the memory consumption for the object.
   @remarks For anaylsing the memory fingerprint.
**/
void taxa::print_getTotal_memoryConsumption_for_object(taxa *listTaxa, uint taxon_length, FILE *f) {
  assert(f);
  loint size_reserved =0;
  if(taxon_length && listTaxa) {
    for(uint i = 0; i < taxon_length; i++) {
      size_reserved += listTaxa[i].get_object_memory_size();
    }
  }
  const float size_gb = (float)size_reserved/(1024*1024*1024);
  char end_char = ' ';
  // PRints the data, but in order to avoid clutter, newline is avoided if the goal is printing the data to the log file:
  if((f==stdout) || (f == stderr)) end_char = '\n';
  fprintf(f, "- listTaxa has %lld B of memory =~ %.5f GB %c", size_reserved, size_gb, end_char); 
}


/**!
   @Name: setTaxonName -- Copies the taxon name into the argument, and updatevs the arguments internal position.
   @remarks: increments the pointer.
   @Changed: 25.08.2011 by oekseth (documentation).
*/
void taxa::setTaxonName(char *&string, bool increment_pointer) {
  if(name) {
    const uint name_length = strlen(name);
    if(!string) {
      string = new char[name_length+1];
      log_builder::test_memory_condition_and_if_not_abort(string!=NULL, __LINE__, __FILE__, __FUNCTION__);
      assert(string); // Verifies memory is allocated.
      string[name_length] = '\0';
    }
    strncpy(string, name, name_length);
    if(increment_pointer) string += name_length;
  } else name = NULL;
}

/**
   @brief Sets the protein label into the list, and updates the internal position of the string to its end
   @param <local_index> If index is out of bound, the input string is set to NULL.
   @param <string> Set to NULL if index not found.
**/
void taxa::setProteinName(uint local_index, char *&string) {
  uint length = 0;
  bool end = false; // Memory organized continously, someetimes with onlye one '\0' mark before the next protein starts
  if(arrKey && ((int)local_index < total_cnt) && arrKey[local_index]) {
    //  for(uint k = 0;k < strlen(arrKey[i]) && !end; k++) {
    for(uint k = 0;k < strlen(arrKey[local_index]) && !end; k++) {
      if(arrKey[local_index][k] != ' ' && arrKey[local_index][k] != '\0') {
	if(arrKey[local_index][k] != '\0') {
	  length = k+1;
	} else end = true;
      } else end = true;
    }

    strncpy(string, arrKey[local_index], length);
    string += length;
  } else string = NULL; // Not set
}


/**!
   @Name: getCompleteProteinName(..) -- protein label and taxon label.
   Returns the complete protein name, to be used in the mcl output
   @Changed: 25.08.2011 by oekseth (documentation)
*/
char *taxa::getCompleteProteinName(uint local_index) {
  if(local_index < (uint)total_cnt) {
    char* string = new char[40]; 
    log_builder::test_memory_condition_and_if_not_abort(string!=NULL, __LINE__, __FILE__, __FUNCTION__);
    memset(string, '\0', 40);
    char *end_of_string = string;
    if(INDEX_IN_FILE_FOR_TAXON_NAME > INDEX_IN_FILE_FOR_PROTEIN_NAME) {
      setProteinName(local_index, end_of_string);
      *end_of_string++ = SEPERATOR;
      setTaxonName(end_of_string, true);
    } else {
      setTaxonName(end_of_string, true); 
      *end_of_string++ = SEPERATOR;
      setProteinName(local_index, end_of_string);
    }
    return string;
  } //else fprintf(stderr, "!!\tLocal_index(%u) above total_cnt(%u). Returns NULL in %s at line %d\n", local_index, total_cnt, __FILE__, __LINE__);
  return NULL;
}
 
/**Returns the taxon id for the given 'world key' used as parameter */
static uint getTaxonIndex(taxa *listTaxa, uint world_key,uint taxon_length) {
  for(uint taxon_id = 0; taxon_id < taxon_length; taxon_id++) {
    if((listTaxa[taxon_id].getRelativeStartingIndex() <= (int)world_key) && ((int)world_key < listTaxa[taxon_id].getRelativeEndIndex())) {
      return taxon_id;
    }
  }
  fprintf(stderr, "!!\t\tError: Did not find taxon id for world_key=%u. Returns 0\n", world_key);
  return 0;
}

/**     @return the protein label for the given argument id.  */
char *taxa::getArrKey(uint protein_id) {
  assert(arrKey);
  assert((int)protein_id < total_cnt);
  if((int)protein_id < total_cnt) {
    return arrKey[protein_id];
  } else {return NULL;}
}

/**
   @return the protein label for the given argument id.
   @param <world_protein_id> The index given using the global id scheme.
*/
char *taxa::getArrGlobalKey(int world_protein_id) {
  assert(arrKey);
  assert(world_protein_id >= rel_start);
  assert(world_protein_id < rel_end);
  const int protein_id = (int)world_protein_id - rel_start;
  //    fprintf(stderr,"protein_id(%d) = %u - %d\n", protein_id, world_protein_id, rel_start);
  assert(protein_id < total_cnt);
  assert(protein_id >= 0);
  return arrKey[protein_id];
}

/**
   @Description:  Sets the name of the taxon
   @ProgramPoint: In the first (of two) parsing operations
   through the input blast file
*/
void taxa::setName(char *_name) {
  if(_name) {
    if(name) delete [] name; // A simple proecdreu verifying that it is as it should
    const loint size = strlen(_name);
    name = new char[size+1]; name[size] = '\0';
    log_builder::test_memory_condition_and_if_not_abort(name!=NULL, __LINE__, __FILE__, __FUNCTION__);
    strncpy(name, _name, size);
  }
}


/**
   @brief Generate (writes) a log file describing the details of the memory signature for taxa collected during parsing
   @remarks 
   - Useful for analyzing- optimizing algorithmic correctness- and memory allocation procedures.
   - Deacitvated if in NDEBUG mode
**/
void taxa::generate_memory_allocation_overview(char *file_name, taxa *listTaxa, int taxon_length) {
#ifndef NDEBUG
  char str[200]; memset(str, '\0', 200);
  bool path_ok = true;
#ifdef LOG_FOLDER_NAME
  // The default procedure we use:
  struct stat st;
  if(stat(FUNC_LOG_FOLDER_NAME,&st) != 0) { // The path does not already exists
    if(-1 == mkdir(FUNC_LOG_FOLDER_NAME, S_IRWXU)) {path_ok = false;}
  }
#endif
  if(path_ok) {
#ifdef LOG_FOLDER_NAME
    char sep = '_';
    if(sizeof(FUNC_LOG_FOLDER_NAME) > 0) sep = '/';
    sprintf(str, "%s%c%s.log", FUNC_LOG_FOLDER_NAME, sep, "taxa_list");
#else
    sprintf("%s.log", "taxa_list");
#endif
    FILE *f_log = fopen(str, "w");
    if(f_log) {
      struct tms tmsstart; clock_t clock_log;
      if((clock_log = times(&tmsstart)) == -1) // starting values
	fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
      uint size = 0; for(int i = 0; i<taxon_length; i++) size+= listTaxa[i].get_object_size(taxon_length);
      fprintf(f_log, "Prints the %d taxa found for file %s, with the list of taxa object consuming a total size of %f.2MB\n",
	      taxon_length, file_name,  
	      (float)size/(float)(1024*1024));
      for(int i = 0; i<taxon_length; i++) listTaxa[i].print_variables(f_log);
      // Generates the estimate on the time consumption:
      struct tms tmsend; clock_t end; long clktck = 0;
      if((end = times(&tmsend)) == -1) fprintf(stderr, "!!\tAn error in measurement of time consumption at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
      if((clktck =  sysconf(_SC_CLK_TCK)) < 0) fprintf(stderr, "!!\tAn error in sys-configuration at line %d in file %s. Contact oekseth@gmail.com for further details.\n", __LINE__, __FILE__);
      const clock_t t_real = end - clock_log;
      const double log_time = t_real/(double)clktck;
      fprintf(f_log, "--\tUsed in total %10.8f seconds for this operation\n", log_time);
      fclose(f_log);
    } else fprintf(stderr, "!!\tUnable to open file %s at line %d found at location %s. Contact oekseth@gmail.com\n",
		   str, __LINE__, __FILE__);
  } else {
    fprintf(stderr, "!!\tIn file %s, do not regard the directory given as a set of alphanumerical chars, therefore disregards the path given: %s\n", __FILE__, FUNC_LOG_FOLDER_NAME); fflush(stderr);
    assert(false);
  }
#endif
}

/**
   @Description:  Returns a list (array) of the proteins 'protein_in'
   has ortholog pairs with
   @param: arrOrtho: the ortho array
   @Param '&size' Sets the size of the array returned
   @Tested: 20.12.2010
*/
uint* taxa::getOrthoRow(uint **arrOrtho, uint protein_in, uint &size) {
  if(arrOrtho !=NULL) {
    if(protein_in < (uint)total_cnt) {
      if(arrOrtho[protein_in] != NULL) {
	size = arrOrtho[protein_in][0]; // The number of elements
	return arrOrtho[protein_in]+1; // to avoid the extra field with the size;
      } else {
	size = 0;
	return NULL;
      }
    }
  } 
  fprintf(stderr, "!!\t'array provided as input=NULL\n");  return NULL;
}


/**
   Takes a key from teh world index and prints the chars representing it
   @param 'real_key': the index given in the world key codes
   @param 'print_taxon': is set to true, prints the taxon with the 'seperator'
*/
void taxa::printProtein(FILE *f, int real_key, const bool print_taxon, char seperator) {
  assert(f);
  if(print_taxon)
    fprintf(f, "%s%c%s", arrKey[real_key-rel_start], seperator, name);
  else fprintf(f, "%d", real_key);//arrKey[real_key-rel_start]);
}
 
/**! Print the list of orthologs based on the array given as input
 */
void taxa::printArrOrthoList(uint **arrOrtho, const uint start, uint end, const bool use_names, taxa *listTaxa, uint taxon_length) {
  uint i = 0;
  for(uint protein_id = start; protein_id< end; protein_id++) {
    if(arrOrtho[i] != NULL) {
      if(use_names) printf("-\t%s has %u orthologs:\t", arrKey[protein_id], arrOrtho[i][0]);
      else printf("-\t%u has %u orthologs:\t", getWorldIndex(protein_id), arrOrtho[i][0]);
      uint size = 0;
      uint *arr = getOrthoRow(arrOrtho, i, size);
      for(uint k = 0; k<size; k++) {
	//if(true)printf("arr[%u]= %u\n", k, arr[k]);
	const uint real_protein_out = arr[k];
	const uint taxon_id = getTaxonIndex(listTaxa, real_protein_out, taxon_length);
	listTaxa[taxon_id].printProtein(real_protein_out, use_names, '_');
	printf(", ");
	//	if(use_names) printf("%s, ", arrKey[k]);	else printf("%u, ", k);
      }
      printf("\n");
    } /*else {
	if(use_names) {printf("!!\t");printProtein(protein_id+rel_start, true, '_');}
	else printf("%u", getWorldIndex(protein_id));
	printf(" has zero ortholog pairs!\n");
	}*/
    i++;
  }
}


/**
   Print the names of the elements in the arrOrtho
*/
/*
  void taxa::printArrOrtho(const bool use_names, taxa *listTaxa, uint taxon_length) {
  for(uint i = 0; i< (uint)total_cnt; i++) {
  if(arrOrtho[i] != NULL) {
  if(use_names) printf("-\t%s has %u orthologs:\t", arrKey[i], arrOrtho[i][0]);
  else printf("-\t%u has %u orthologs:\t", getWorldIndex(i), arrOrtho[i][0]);
  uint size = 0;
  uint *arr = getArrOrthoRow(i, size);
  for(uint k = 0; k<size; k++) {
  const uint real_protein_out = arr[k];
  const uint taxon_id = getTaxonIndex(listTaxa, real_protein_out, taxon_length);
  listTaxa[taxon_id].printProtein(real_protein_out, use_names, '_');
  printf(", ");
  }
  printf("\n");
  }
  } 
  }
*/
/** @Description:  Initiates the *arrInpaLimit
    @ProgramPoint: Is first needed in the 'ortho_set_bin', therefore
    not malloced before at this point
*/
void taxa::initArrInpaLimit(bool only_if_empty) {
  if(!only_if_empty || (arrInpaLimit == NULL)) {
    arrInpaLimit = (float*)malloc(sizeof(float)*total_cnt);
    for(uint i = 0; i< (uint)total_cnt; i++) arrInpaLimit[i]=0.0;
  }
}

//! @return true if protein is accepted:
bool taxa::aboveInparalogLimit(uint protein_id, const float averaged_similarity, const bool INPARALOG_OPERATION) {
  assert((int)protein_id < total_cnt);
  if(INPARALOG_OPERATION) {
    // Uses the first alternative, as we regard it as most proper:
    if(true) { // Either greather than- or similar
      if(arrInpaLimit[protein_id] <= averaged_similarity) return true; // the sim score is above: returns true
      else {return false;    }
    } else {
      if(arrInpaLimit[protein_id] < averaged_similarity) return true; // the sim score is above: returns true
      else {return false;    }
    }
  } else return true;
}


//! Prints the limits (i.e. the similarity scores) whom the inpralogs must be above
void taxa::printArrInpaLimit(bool use_names) {
  printf("---\nThe limits whom the inparalogs for taxon '%s' must be above is as below:\n", name);
  for(uint i = 0; i< (uint)total_cnt; i++) {
    printf("\t");
    if(use_names) printf("%s ", arrKey[i]);
    else printf("%u ", i+rel_start);
    printf("has limit '%f'\n", arrInpaLimit[i]);
  }
}
/**
   @Description:  Updates in 'sim_score' arbuemnt is higher than the
   previous value set at position 'protein_in'
   @ProgramPoint: In 'ortho_set_bin'		   
*/
void taxa::insertInpaLimit(const uint protein_in, float sim_score) {
  if((int)protein_in < total_cnt) {
    if(arrInpaLimit[protein_in] < sim_score) arrInpaLimit[protein_in] = sim_score;
  } else fprintf(stderr, "taxa\tprotein_in(%u) !<ttoal(%u)\n", protein_in, total_cnt);
}

/** @Description: Frees 'arrInpaLimit'
    @ProgramPoint: Before collecting the data to an output file
*/
void taxa::delete_initArrInpaLimit() {free(arrInpaLimit);}

/** @Description: Frees 'arrInpaLimit'
    @ProgramPoint: Before collecting the data to an output file
*/
/*
  void taxa::delete_initArrOrtho() {
  if(arrOrtho != NULL) {
  for(uint i = 0; i< (uint)total_cnt; i++)
  free(arrOrtho[i]);
  free(arrOrtho); arrOrtho = NULL;
  } 
  }
*/  

void taxa::printOverlap(const bool use_name) {
  if(arrOverlap != NULL) {
    printf("#\tPrinting the protein lengtht for taxon %s\n", name);
    for(uint i = 0; i< (uint)total_cnt; i++) {    
      if(!use_name) printf("\t[%d]\t%d\n", i, arrOverlap[i]);
      else printf("\t%s\t%d\n", arrKey[i], arrOverlap[i]);
    }
  } else printf("\n'arrOverlap' se not set\n");
}

//! Returns a local key to be used in this subset of the data set
//!   @Tested: 20.12.2010
uint taxa::getLocalIndex(uint world_key) {
  if(rel_start > (int)world_key || (rel_end <= (int)world_key)) {
    fprintf(stderr, "!!\tError: world_key=%u, rel_start=%u, giving local_key=%u\n", world_key, rel_start, (world_key - rel_start));
    return UINT_MAX;
  }
  return (world_key - rel_start);
}

//! Returns the unique key tbu in the whole data set
//!   @Tested: 20.12.2010
uint taxa::getWorldIndex(uint local_key) {return (local_key + rel_start);}
 
//! Returns true if the world key given as input belongs to this taxon
bool taxa::isRelation(int world_key) { return (
					       (rel_start <= world_key) &&
					       (rel_end   >  world_key));
}
/*
  bool taxa::hasArrOrthoRow(uint protein_in) {
  if(arrOrtho !=NULL) {
  if(protein_in < (uint)total_cnt) {
  if(arrOrtho[protein_in] != NULL) {
  return true;
  }
  }
  }
  return false;
  }
*/
void taxa::print_variables() {
  print_variables(stdout);
}

void taxa::print_variables(FILE *f) {
  fprintf(f, "%s is in [%d, %d], with a total of %d\n", name, rel_start, rel_end, total_cnt);
}
//! Prints hte key name pair for this taxon
void taxa::printArrKeyList(int taxon_index) {
  printf("\n--\nPrints the key - name pairs for the taxon %s (at taxon_index=%d):\n", name, taxon_index);
  for(uint i =0; i< (uint)total_cnt; i++) {
    printf("[%d]  =  %s\n", i, arrKey[i]);
  }
}

//! Returns the overlap using the local index
uint taxa::getOverlap(uint local_index) {
  if(arrOverlap != NULL)  {
    if(local_index < (uint)total_cnt)  return arrOverlap[local_index];
    else printf("!!\t\tlocal index = %u > total_collection=%u. REturns 0.\n", local_index, total_cnt);
  } else printf("!!\t\t'arrOverlap' not intialized!\n");
  return 0;
}
 
//! Returns the taxon id based on the inputs given
uint taxa::getTaxonId(uint world_key, taxa *listTaxa, uint taxon_length) {
  for(uint i = 0; i< taxon_length; i++) {
    if(listTaxa[i].isRelation(world_key)) return i;
  }
  fprintf(stderr, "!!\tAn error occured: world_key %u did not exixts in the collection of %u taxon. Aborts!\n", world_key, taxon_length); exit(2);
}

//! Returns the name of the protein, given as a world index
char *taxa::getProteinNameWorld(uint world_ind_in) {
  const uint local_ind_in = getLocalIndex(world_ind_in);
  if(local_ind_in != UINT_MAX) return arrKey[local_ind_in];
  else return NULL;
}
/**
   @Description: Frees the memory allocated
   @ProgramPoint: To be done at the end of the program:
*/
void taxa::delete_buffer() {
  //TODO: See if the outcomment below is to be commented in
  if(arrKey) {
    for(loint i = 0; i < (loint)total_cnt; i++) {
      delete [] arrKey[i]; arrKey[i] = NULL;
    }
    free(arrKey); arrKey = NULL;
  }
  if(arrInpaLimit) {
    free(arrInpaLimit); arrInpaLimit = NULL;
  }
  if(arrOverlap) {
    free(arrOverlap); arrOverlap = NULL;
  }
  if(name) {delete [] name; name = NULL;}
}

void taxa::delete_taxa(int &taxon_size, taxa *&list) {
  if(list) {
    for(int i = 0; i < taxon_size; i++) list[i].delete_buffer();
    delete [] list; list = NULL; taxon_size = 0;
  }
}

void taxa::print_class_info(FILE *out) {
  if(out == NULL) out = stdout; 
  if(true) {
    //fprintf(out, "value at line %d in taxa.cxx is %s\n", __LINE__, name);
    fprintf(out, "value at line %d in taxa.cxx is %d\n", __LINE__, total_cnt);
    fprintf(out, "value at line %d in taxa.cxx is %d\n", __LINE__, rel_start);
    fprintf(out, "value at line %d in taxa.cxx is %d\n", __LINE__, rel_end);
  }
  fprintf(out, "-\tObject taxa of name '%s' has %d proteins in interval[%d, %d]\n", name, total_cnt, rel_start, rel_end);
}
void taxa::print_class_info(FILE *out, taxa *listTaxa, int taxon_length) {
  if(out == NULL) out = stdout;
  fprintf(out, "\n-The taxa-list-of-objects listTaxa has %d elements with the main-properties:\n", taxon_length);
  for(int i = 0; i < taxon_length; i++) {
    listTaxa[i].print_class_info(out);
  }
}

taxa::taxa() :
  SEPERATOR('_'), INDEX_IN_FILE_FOR_PROTEIN_NAME(0), INDEX_IN_FILE_FOR_TAXON_NAME(1),
  total_cnt(0),rel_start(0),rel_end(0), name(NULL), arrKey(NULL),
  arrOverlap(NULL),
  //  arrOrtho(NULL),
  arrInpaLimit(NULL)
{}

taxa::taxa(char *_name, int _total_cnt, int _rel_start, int _rel_end, char **_arrKey, char sep, int index_protein, int index_taxon) :
  SEPERATOR(sep), INDEX_IN_FILE_FOR_PROTEIN_NAME(index_protein), INDEX_IN_FILE_FOR_TAXON_NAME(index_taxon),
  total_cnt(_total_cnt), rel_start(_rel_start), rel_end(_rel_end),  name(NULL), arrKey(_arrKey),
  arrOverlap(NULL),  //arrOrtho(NULL),
  arrInpaLimit(NULL)
{
  if(_name != NULL) {
    setName(_name);
  }
  //  else _name = NULL;
}



//! Purpose of this code is to provide a default length of each taxon to be used
taxa * taxa::init_test_array(uint taxon_size, uint protein_size) {
  taxa *list = new taxa[taxon_size];
  log_builder::test_memory_condition_and_if_not_abort(list!=NULL, __LINE__, __FILE__, __FUNCTION__);
  for(uint i = 0; i < taxon_size; i++) {
    list[i] = taxa();
    list[i].total_cnt = protein_size;
  }
  return list;
}

//! @return a char pointer holding the taxa names, given a taxa object.
char **taxa::get_taxa_names(taxa *obj, int taxon_length) {
  if(obj && taxon_length) {
    char **buffer = new char*[taxon_length];
    log_builder::test_memory_condition_and_if_not_abort(buffer!=NULL, __LINE__, __FILE__, __FUNCTION__);
    assert(buffer); // Verifies memory is allocated.
    for(int i =0; i < taxon_length; i++) {
      buffer[i] = NULL;
      obj[i].setTaxonName(buffer[i], false);
    }
    return buffer;
  } else return NULL;
}


//! Sets the proteins above the given threshold.
void taxa::get_cnt_taxons_above_limit(taxa *obj, int taxon_length, int protein_limit,  loint &cnt_proteins_temp, loint &cnt_taxa_temp) {
  cnt_proteins_temp = 0, cnt_taxa_temp = 0;
  if(obj && taxon_length) {
    for(int i =0; i < taxon_length; i++) {
      if((obj[i].total_cnt > protein_limit)) {
	cnt_taxa_temp++; cnt_proteins_temp +=obj[i].total_cnt;
      }	    
    }
  }
}

/**
   @brief Used validating/comparing results, ie, before- and after MPI transfer.
   @return the total score
**/
float taxa::get_total_score_for_arrInpaLimit(taxa *&listTaxa, uint taxon_length) {
  float sum = 0;
  for(uint i = 0; i<taxon_length; i++) {
    for(uint p = 0; p < (uint)listTaxa[i].total_cnt; p++) {
      if(listTaxa[i].arrInpaLimit) {
	sum += listTaxa[i].arrInpaLimit[p];
      }
    }
  }
  return sum;
}

// ---------------------------------------------------------------------------------
// Below sepcific routines used in MPI transfer:
// ---------------------------------------------------------------------------------
#ifdef USE_MPI

//
// ----------- Below calls used for templatisation of the transformation (2d->1d->2d)
//

uint *get_overlap_list(taxa *listTaxa, const uint taxon_index) {
  assert(listTaxa);
  return listTaxa[taxon_index].getArrOverlap();
}

float *get_arrInpa_list(taxa *listTaxa, const uint taxon_index) {
  assert(listTaxa);
  return listTaxa[taxon_index].getArrInpaLimit();
}

float get_list_element_overlap(taxa *listTaxa, const uint taxon, const uint index) {
  assert(listTaxa);
  return listTaxa[taxon].getOverlap(index);
}

float get_list_element_arrInpaLimit(taxa *listTaxa, const uint taxon, const uint index) {
  assert(listTaxa);
  return listTaxa[taxon].getArrInpaLimit(index);
}

void copy_arr_overlap_into_mother_list(taxa *listTaxa, const uint taxon_index, uint *array) {
  assert(listTaxa);
  listTaxa[taxon_index].copy_arr_overlap_into_list(array);
}

void copy_arr_inpalimit_into_mother_list(taxa *listTaxa, const uint taxon_index, float *array) {
  assert(listTaxa);
  listTaxa[taxon_index].copy_arr_inpalimit_into_list(array);
}



/**
   @brief Merges the lists into one:
   @remarks The "OP" argument must correspond to the type given in "array_start".
**/
template<typename OP, typename T> void merge_list(taxa *listTaxa, uint index_start, uint index_end_pluss_one, OP get_list, uint &total_length, T *&array_start) {
  assert(!array_start);
  assert(index_end_pluss_one);
  assert(listTaxa);
  const int index_end = index_end_pluss_one-1;
  assert(listTaxa[index_end].getRelativeEndIndex() > listTaxa[index_start].getRelativeStartingIndex());
  total_length = listTaxa[index_end].getRelativeEndIndex() - listTaxa[index_start].getRelativeStartingIndex();
  array_start = new T[total_length];
  log_builder::test_memory_condition_and_if_not_abort(array_start!=NULL, __LINE__, __FILE__, __FUNCTION__);
  T *array_curr = array_start;
  uint cnt_inserted = 0; // Verifying that the list holds itself inside the limiits given.
  for(uint i = index_start; i < index_end_pluss_one; i++) {
    T *temp = get_list(listTaxa, i);
    const uint length = listTaxa[i].getLength();
    if(length && temp) {
      memcpy(array_curr, temp, sizeof(T)*length);
      array_curr += length; // Updates the pointer
      cnt_inserted += length;
      assert(cnt_inserted <= total_length);
    }
  }
}


/**
   @brief Asserts the merging of the list:
   @remarks: If "before_mpi" is set, uses this element as basis for the comparison; alternative is the internal list found in the listTaxa object used"
**/
template<typename OP, typename OP2, typename T> void assert_merging_of_lists(taxa *listTaxa, uint index_start, uint index_end_pluss_one, OP get_list, OP2 get_list_element, uint total_length, T *array_start, T *&before_mpi) { 
#ifndef NDEBUG  
  assert(array_start);
  T *array_curr = array_start;
  uint cnt_inserted = 0; // Verifying that the list holds itself inside the limiits given.
  T *overlap_before_mpi_curr = before_mpi;
  for(uint taxon_id = index_start; taxon_id < index_end_pluss_one; taxon_id++) {
    T *temp = overlap_before_mpi_curr;
    if(!before_mpi) temp = get_list(listTaxa, taxon_id); //listTaxa[i].getArrOverlap();
    const uint length = listTaxa[taxon_id].getLength();
    if(length && temp) {
      for(uint protein_id = 0; protein_id < length; protein_id++) {
	if(before_mpi) {
	  assert(temp[protein_id] <= get_list_element(listTaxa, taxon_id, protein_id)); 
	} else 	assert(temp[protein_id] == get_list_element(listTaxa, taxon_id, protein_id)); 
      }
      overlap_before_mpi_curr += length; // Updates the pointer
      array_curr += length; // Updates the pointer
      cnt_inserted += length;
      assert(cnt_inserted <= total_length);
    }
  }
  if(!before_mpi) {
    //! Builds an array for final comparison, verifying the end result:
    before_mpi = new T[total_length];
    log_builder::test_memory_condition_and_if_not_abort(before_mpi!=NULL, __LINE__, __FILE__, __FUNCTION__);
    memcpy(before_mpi, array_start, sizeof(T)*total_length);
  } else {
    //! Deallocates the objects:
    delete [] before_mpi; before_mpi = NULL;
  }
#endif
}

//! The expander (of the list into several):
template<typename OP, typename T> void expand_1d_list(taxa *listTaxa, int number_of_nodes, uint index_start, uint index_end_pluss_one, OP copy_local_array_into_mother, uint &total_length, T *&array_start) {
  assert(array_start);
  assert(listTaxa);
  assert(index_end_pluss_one);
  for(uint node_id = 1; node_id < (uint)number_of_nodes; node_id++) {
    T *array_curr = array_start;
    uint cnt_inserted = 0;
    for(uint i = index_start; i < index_end_pluss_one; i++) {
      const uint length = listTaxa[i].getLength();
      copy_local_array_into_mother(listTaxa, i, array_curr);
      //    listTaxa[i].copy_arr_overlap_into_list(array_curr);
      array_curr += length; // Updates the pointer
      cnt_inserted += length;
      assert(cnt_inserted <= total_length);    
    }
  }
}

/**
   @brief Sends- and receives the list of norm_t* objects accross nodes.
**/
void taxa::mpi_send_and_receive_arrOverlap(MPI_Comm comm, int myrank, int number_of_nodes, taxa *&listTaxa, uint index_start, uint index_end_pluss_one) {
  assert(number_of_nodes);
  //! Get the 1d list:
  uint total_length = 0; uint *array_start = NULL; merge_list(listTaxa, index_start, index_end_pluss_one, get_overlap_list, total_length, array_start);
  assert(array_start);
#ifndef NDEBUG
  //! Verfies the copying process, and makes a deep copy of the list before merging:
  uint *overlap_before_mpi= NULL; 
  assert_merging_of_lists(listTaxa, index_start, index_end_pluss_one, get_overlap_list, get_list_element_overlap, total_length, array_start, overlap_before_mpi);
  uint sum_overlaps_before_merging = taxa::get_total_overlap_score(listTaxa, index_start, index_end_pluss_one);
  uint total_sum_before = 0;
  MPI_Allreduce(&sum_overlaps_before_merging, &total_sum_before, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
  const uint sum_before_myrank = sum_overlaps_before_merging;
  //  const uint sum_before = total_sum_before;

  //! Validates that the sum function is correct:
  uint validated_sum = 0;
  int elements_iterated = 0;
  for(uint i = index_start; i < index_end_pluss_one; i++) {
    uint *arr = listTaxa[i].getArrOverlap();
    assert(arr);
    const uint size = (uint)listTaxa[i].total_cnt;
    for(uint out = 0; out < size; out++) {
      validated_sum += arr[out];
      elements_iterated++;
    }
  }
  //  if(myrank == 0) 
  assert(validated_sum == sum_before_myrank);
  assert(elements_iterated == listTaxa[index_end_pluss_one-1].rel_end);
#endif

  //! The mpi sending procedure, summing the values using the same pointer for the task:
  MPI_Allreduce(MPI_IN_PLACE, array_start, total_length, MPI_UNSIGNED, MPI_SUM, comm );

  //! The expander (of the list into several):
  expand_1d_list(listTaxa, number_of_nodes, index_start, index_end_pluss_one, copy_arr_overlap_into_mother_list, total_length, array_start);

#ifndef NDEBUG
  //! Asserts that the sum of the number of inserted elements are greater- or equal the sum before the merging operation:
  assert(sum_overlaps_before_merging <= taxa::get_total_overlap_score(listTaxa, index_start, index_end_pluss_one));
  assert_merging_of_lists(listTaxa, index_start, index_end_pluss_one, get_overlap_list, get_list_element_overlap, total_length, array_start, overlap_before_mpi);  

  //! Verifies the sums are equal for all the nodes.
  uint sum_myrank_after =  get_total_overlap_score(listTaxa, index_start, index_end_pluss_one);
  uint total_sum = 0;
  MPI_Allreduce(&sum_myrank_after, &total_sum, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
  assert(sum_myrank_after <= total_sum);

  assert(sum_myrank_after == total_sum_before); // myrank shall have gotten all the data.

  // TODO: Usnure if problems could arise due to "rounding":
  if(!((sum_myrank_after * number_of_nodes) == total_sum)) {
    char str[100]; memset(str, '\0', 100);
    sprintf(str, "For myrank[%d] arrInpaCommunication, where [sum_myrank_after(%u)*number_of_nodes(%d) =! MPI_SUM(%u)", myrank, sum_myrank_after, number_of_nodes, total_sum);
    log_builder::throw_warning(software_error, __LINE__, __FILE__, __FUNCTION__,str);
    assert((sum_myrank_after * number_of_nodes) == total_sum);
  }
#endif

  //! De-allocates the data when tests have completed:
  delete [] array_start; array_start = NULL; // Frees the memory.
}

/**
   @brief Sends- and receives the list of arrInpaLimit* objects accross nodes.
**/
void taxa::mpi_send_and_receive_arrInpaLimit(taxa *&listTaxa, uint taxon_length) {
  assert(listTaxa);
  assert(taxon_length);
  //#ifndef NDEBUG
  //! Some MPI variables:
  int myrank = 0;  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
  int number_of_nodes = 1; MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes); 
  //#endif


  //! Get the 1d list:
  uint total_length = 0; float *array_start = NULL; merge_list(listTaxa, 0, taxon_length, get_arrInpa_list, total_length, array_start);

#ifndef NDEBUG
  //! Verfies the copying process, and makes a deep copy of the list before merging:
  float *before_mpi= NULL; 
  assert(array_start);
  assert_merging_of_lists(listTaxa, 0, taxon_length, get_arrInpa_list, get_list_element_arrInpaLimit, total_length, array_start, before_mpi);
  const float sum_arrInpa_before_merging = get_total_score_for_arrInpaLimit(listTaxa, taxon_length);
#endif

  //! The mpi sending procedure, summing the values using the same pointer for the task:
  MPI_Allreduce(MPI_IN_PLACE, array_start, total_length, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );

  //! The expander (of the list into several):
  expand_1d_list(listTaxa, number_of_nodes, 0, taxon_length, copy_arr_inpalimit_into_mother_list, total_length, array_start);

#ifndef NDEBUG
  //! Asserts that the sum of the number of inserted elements are greater- or equal the sum before the merging operation:
  assert(sum_arrInpa_before_merging <= get_total_score_for_arrInpaLimit(listTaxa, taxon_length));
  assert_merging_of_lists(listTaxa, 0, taxon_length, get_arrInpa_list, get_list_element_arrInpaLimit, total_length, array_start, before_mpi);

  //! Verifies the sums are equal for all the nodes.
  float sum_myrank =  get_total_score_for_arrInpaLimit(listTaxa, taxon_length);
  float total_sum = 0;
  MPI_Allreduce(&sum_myrank, &total_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
  assert(sum_myrank <= total_sum);
    //! Round to integers in order to avoid rounding errors:
  if(!((int)(sum_myrank * number_of_nodes) == (int)total_sum)) {
    char str[100]; memset(str, '\0', 100);
    sprintf(str, "For myrank[%d] arrInpaCommunication, where [sum_myrank(%f)*number_of_nodes(%d) =! MPI_SUM(%f)", myrank, sum_myrank, number_of_nodes, total_sum);
    log_builder::throw_warning(software_error, __LINE__, __FILE__, __FUNCTION__,str);
    assert((sum_myrank * number_of_nodes) == total_sum);
  }
#endif
  //! De-allocates the data when tests have completed:
  delete [] array_start; array_start = NULL; // Frees the memory.
}

#endif
