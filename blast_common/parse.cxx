#include "parse.h"

//! @return true if it's a char to be used
  bool Parse::is_legal_char(char c) {
  if(c == '\n') return false;
  else if (c == '\0') return false;
  else return true;
}

//! @return the correct length of the protein by adjusting the input given.
  loint Parse::get_correct_length(char *pos_start, loint length_name) {
   bool at_end = false;
  for(loint i = 0; i < length_name && !at_end; i++) {
     if(!is_legal_char(pos_start[i])) {
       length_name = i+1; at_end = true;
     }
   }
   return length_name;
 }

 //! Sets the name:
  void Parse::set_string(char *&string, char *pos_start, loint length_name) {
   length_name = get_correct_length(pos_start, length_name);    // Gets the correct length of the protein
   if(length_name) {
     char *tmp = strchr(pos_start, '\n');
     if(tmp) pos_start = tmp +1; // If a newline is found, we jump past it.
#ifndef NDEBUG
     for(int i = 0; i < (int)strlen(pos_start); i++) {
       assert(pos_start[i] != '\n');
     }
#endif
     if(!string) string = new char[length_name+1]; string[length_name] = '\0';
     strncpy(string, pos_start, length_name);
   } else string = NULL;
 }

  //! @return the inner name
char *Parse::get_name_in() {
  char *pos_start = name_in;
  if(pos_start) {
    char *tmp = strchr(pos_start, '\n');
    if(tmp) pos_start = tmp +1; // If a newline is found, we jump past it.
#ifndef NDEBUG
     for(int i = 0; i < (int)strlen(pos_start); i++) {
       assert(pos_start[i] != '\n');
     }
#endif
    return pos_start;
  } else return NULL;
}
/**@brief Sets the name given the input       **/
void Parse::setName(char *start_pos, char *&name, uint index, uint *array, uint array_size) {
  char *name_start_pos = start_pos;
#ifndef NDEBUG
  // Adds extra debugging information if the software breaks at this point:
  if(!((index == 0) || (array_size >= 2))) {
    char string[100]; memset(string, '\0', 100);
    sprintf(string, "!!\tDid not find enough blocks in the labels (an identifier seperated the seperator given by the user) to process the given blasp-row; we got index(%u) and array_size(%u)", index, array_size);
    log_builder::throw_warning(software_dependencies, __LINE__, __FILE__, __FUNCTION__,string);
    assert((index == 0) || (array_size >= 2)); 
  }
#endif

  if(index!=0) name_start_pos += array[index-1]+1;
  uint name_length = array[index]; // to avoid the trailing chars
  bool legal_format = true;
  if(index!=0) {
    if(name_length == 0) {
      legal_format = false;
    } else {
      name_length -= array[index-1]-1;
      if(index!=array_size-1) name_length -= 2;
    }
  }
  if(index == 0) {
    //! Could have a newline ('\n') residing in it: evaluates- and if neccessary removes:
    uint new_length = name_length;
    char *new_name_start_pos_start = name_start_pos;
    uint cnt_changed = 0;
    for(uint i = 0; i < (uint)name_length; i++) {
      if(name_start_pos[i] == '\n') {new_name_start_pos_start++, cnt_changed++; new_length--;}
    }
    assert(cnt_changed <= 1);
    if(cnt_changed == 1) {
      //! Updates, ie, 'removes' the newline found:
      name_start_pos = new_name_start_pos_start;
      name_length = new_length;
    }
  }
#ifndef NDEBUG
  //! Asserts that no '\n' is included in the string sent:
  uint correct_pos = 0;
  uint found_newline_at_index = UINT_MAX;
  for(uint i = 0; i < (uint)name_length && (i == correct_pos); i++) {
    if(name_start_pos[i] != '\n') correct_pos++;
    else found_newline_at_index = i;
  }
  if(!(correct_pos == (uint)name_length)) {
    fprintf(stderr, "!!\t found_newline_at_index(%u), given a total_length_of_string(%u), for '%s', at %s:%d\n",
	    found_newline_at_index, name_length, name_start_pos, __FILE__, __LINE__);
	    
    assert(correct_pos == (uint)name_length);
  } else assert(found_newline_at_index == UINT_MAX);
#endif

  if(legal_format) {
    // Rellocates if it's not long enough:
    if(name && !(name_length < strlen(name))) {delete [] name; name = new char[name_length+1];}
    else if(!name) {name = new char[name_length+1];}
    assert(name);
    name[name_length] = '\0'; // Sets the end of the string.
    strncpy(name, name_start_pos, name_length);
  } else {name = NULL;}

#ifndef NDEBUG
  if(name == NULL) {
    fprintf(stdout, "!!\tName not set for index(%u) ", index);
    if(legal_format && !(name_length < 20)) {
      fprintf(stdout, "Error is caused by name beeing too long (%u chars), but the limit is < %d. As chars they are:", name_length, (int)SIZE_PROTEIN);
      for(uint i = 0; i < name_length; i++) printf("[%u]=%c,\t", i, name_start_pos[i]); printf("\n");
    }
    fprintf(stdout, "!!\t at line %d in %s: ",  __LINE__, __FILE__);
    for(uint i = 0; i < array_size; i++) printf("[%u]=%u\t", i, array[i]); printf("\n");
  }
#endif
}

//! Returns an indexed string based on a line from the blast file.
uint *Parse::buildArray(char *start_pos, char *&pos_column_end, char seperator, uint &size_array, char *logical_end) {
  assert(start_pos);
  assert(logical_end);
  size_array = 10; uint *array = new uint[size_array];
  char *array_pos = start_pos;
  for(uint i = 0; i< size_array; i++) array[i] = 0;
  uint array_index = 0;
  while((*array_pos  != '\t') && ((*array_pos  !=' ') || (array_pos >= logical_end)) ) { // continues until column end is reached
    if(*array_pos == seperator) {
      const uint distance = array_pos-start_pos; 
      array[array_index++] = distance;
      if(array_index == (size_array-1)) { // enlarge the array: '-1' due to the need of adding the last size to the end
	const uint new_size = 2*size_array;
	uint *new_array = new uint[new_size];
	for(uint i = 0; i< size_array; i++) new_array[i] = array[i]; // copies the data
	for(uint i = size_array; i< new_size; i++) new_array[i] = 0;

	delete [] array;
	array = new_array, size_array = new_size;
      }
    }
    array_pos++;
  }
  array[array_index++] = array_pos-start_pos-2; // sets the final length: '-1' to avoid the delimiter
  size_array = array_index;
  pos_column_end = array_pos;
  return array;
}


/**
   @brief Returns the end of the column if the data contains all teh filds: retursn NULL if itdoes not contains all the data
   @param <first_column> True if it's the leftmost protein in the given blast-file line
   @param <start_pos> The memory position to start the work from.
   @param <index_name> The index in the blast file where the protein name is found
   @param <index_taxon> The index in the blast file where the taxon name is found
   @param <seperator>   Field separator in FASTA (blast) headers dividing the taxon from the protein label.
   @param <logical_end> The last address where data may be find; used to avoid errors.
   @param <array_size> An initial guess on how many columns there are in the blast file with regard to the name-label.
**/
char *Parse::set_names(const bool first_column, char *start_pos, uint index_name, uint index_taxon, char seperator, char *logical_end, uint array_size) {
  char *column_end = NULL;
  // The array below used holding an index of the protein label:
  uint *array = buildArray(start_pos, column_end, seperator, array_size, logical_end);
  // array[index] corresponds to the length from the start of the given string
  char *name = NULL, *taxon = NULL; // For a more standarised approach to validating the source.
  if(first_column) {setName(start_pos, name_in, index_name, array, array_size); name = name_in;}
  else             {setName(start_pos, name_out,index_name, array, array_size); name = name_out;}
  if(name != NULL) {// data set for the protein name
    if(first_column) {setName(start_pos, taxon_in,  index_taxon, array, array_size); taxon = taxon_in;}
    else             {setName(start_pos, taxon_out, index_taxon, array, array_size); taxon = taxon_out;}
    if(taxon != NULL) { // data set for the taxon
      const uint size_last = array[array_size-1];
      delete [] array; array = NULL;
      return (start_pos + size_last) +3;
    } 
#ifndef NDEBUG
    else  {
      printf("!!\ttaxon not found in column %d for protein_label(%s) with seperator(%c), index_name(%u), index_taxon(%u), taxon_in(%s), taxon_out(%s), name_in(%s), name_out(%s) at line %d in %s\n", first_column, name, seperator, index_name, index_taxon, taxon_in, taxon_out, name_in, name_out, __LINE__, __FILE__);
    }
#endif
  }
#ifndef NDEBUG
  else  {
    printf("!!\tprotein_label not found in column %d for protein_label(%s) with seperator(%c), index_name(%u), index_taxon(%u), taxon_in(%s), taxon_out(%s), name_in(%s), name_out(%s) at line %d in %s\n", first_column, name, seperator, index_name, index_taxon, taxon_in, taxon_out, name_in, name_out, __LINE__, __FILE__);
  }
  
#endif
  delete [] array; array = NULL;
  return NULL;
}

bool Parse::is_equal(char *str, c_id id) {
  if(id == n_in) return (0 == strcmp(name_in, str));
  else if(id == n_out) return (0 == strcmp(name_out, str));
  else if(id == t_in) return (0 == strcmp(taxon_in, str));
  else if(id == t_out) return (0 == strcmp(taxon_out, str));
  else return false;
}

//! Prints info about the class:
void Parse::print_line(char SEPERATOR) {
  printf("\tInner protein: (%s)[%s]; ", taxon_in, name_in);
  printf("Outer protein: (%s)[%s]; ", taxon_out, name_out);
  printf("Overlap_in: %d; ", overlap_in);
  printf("Overlap_out: %d; ", overlap_out);
  printf("Distance: %f; ", distance);    
  //  printf("  %s%c%s\t%s%c%s\t%d\t%d\t%f\t", name_in, SEPERATOR, taxon_in, name_out, SEPERATOR, taxon_out, overlap_in, overlap_out, get_distance());
}
//! Prints info about the class
void Parse::print(char SEPERATOR) {
  if(true) {
    printf("\tInner protein: (%s)[%s]; ", taxon_in, name_in);
    printf("Outer protein: (%s)[%s]; ", taxon_out, name_out);
    printf("Overlap_in: %d; ", overlap_in);
    printf("Overlap_out: %d; ", overlap_out);
    printf("Distance: %f; ", distance);    
    printf("\n");
  } else {
    //      printf("%20s%c%-5s%20s%c%-5s%-6d%5.7f\n", name_in, SEPERATOR, taxon_in, name_out, SEPERATOR, taxon_out, overlap, distance);
    printf("%20s%c%-5s %20s%c%-5s %-6d %-6d %5.7f\n", name_in, SEPERATOR, taxon_in, name_out, SEPERATOR, taxon_out, overlap_in, overlap_out, distance);
  }
}

Parse::Parse() : name_in(NULL), name_out(NULL), taxon_in(NULL), taxon_out(NULL) { 
  distance = 0.0;
  overlap_in = 0;
  overlap_out = 0;
}

Parse::Parse(char *ni, char *no, char *ti, char *to, short int overl_in, short int overl_out, float di)
  : name_in(NULL), name_out(NULL), taxon_in(NULL), taxon_out(NULL) {
  set_string(name_in, ni, strlen(ni));
  set_string(name_out, no, strlen(no));
  set_string(taxon_in, ti, strlen(ti));
  set_string(taxon_out, to, strlen(to));
  distance = di;
  overlap_in = overl_in;
  overlap_out = overl_out;
}


void Parse::debug_printBlastFormatOne(FILE *file_input, char SEPERATOR) {
  fprintf(file_input,
	  //  "%s%c%s\t%s%c%s\t33.00\t%u\t5\t8\t\t9\t10\t00\t%e\t12\n",
	  "%s%c%s\t%s%c%s\t33.0022\t%u\t%u%s88\t9\t00\t%e\t12\n",
	  name_in, SEPERATOR, taxon_in,
	  name_out, SEPERATOR, taxon_out,
	  overlap_in, overlap_out,
	  //	  	  "\t",
	  "\t11\t55\t",
	  //	  	  "\t00\t",
	  distance
	  );
}

void Parse::debug_printBlastFormatTwo(FILE *file_input, char SEPERATOR) {
  fprintf(file_input,
	  "%s%ciijunkii%c%s\t%s%ciijunkii%c%s\t33.00\t%u\t%u\t8\t\t9\t10\t00\t%e\t12\n",
	  taxon_in, SEPERATOR, SEPERATOR, name_in,
	  taxon_out, SEPERATOR, SEPERATOR, name_out,
	  overlap_in, overlap_out, distance
	  );
}

void Parse::assert_class(bool print_info) {
  const static char *class_name = "Parse";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  Parse p = Parse();
  // One tests, for the case of symbolity, as this class only contains setters- and getters.
  p.setDistance(5.0);
  assert(p.getDistance() == 5.0);
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

