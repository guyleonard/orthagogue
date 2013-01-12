#include "blast_extractors.h"

/**
   @brief Prints the segment given.
   @param <start_p> The Position in the memory string to start the reading.
   @param <end> The Position in the memory string to end the reading.
   @author Ole Kristian Ekseth (oekseth)
 **/
void blast_extractors::print_segment(char *start_p, char *end) {
  assert(start_p);
  assert(end);
  assert(end < start_p);
  char *start = start_p;
  while(start < end) putchar(*start), start++;
  printf("\n");
}

/**
   @brief Prints each char seperately.
   @param <start_p> The Position in the memory string to start the reading.
   @param <end> The Position in the memory string to end the reading.
   @remarks Sometimes special chars affects the data. This methods primary usage is for visual validation that the input corresponds the the programmers expecation.
   @author Ole Kristian Ekseth (oekseth)
 **/
void blast_extractors::print_segment_stepwise(char *start_p, char *end) {
  assert(start_p);
  assert(end);
  assert(end < start_p);
  char *start = start_p;
  uint i = 0;
  while(start < end) {
    printf("[%u] = '%c'\n", i, *start);
    i++, start++;
  }
}
//! Returns true if it's a char to be used
bool is_legal_char(char c) {
  if(c == '\n') return false;
  else if (c == '\0') return false;
  else return true;
}

//! Iterates through the chars the string consists of and eleminates possible wild chars lurked into the string:
void copy_name(char *destination, loint destination_size, char *source, loint source_size) {
  loint name_index = 0;
  assert(destination_size >= source_size);
  memset(destination, '\0', destination_size);
  for(loint i = 0; i < source_size; i++) {
    if(is_legal_char(source[i])) destination[name_index++] = source[i];
  }
  
#ifndef NDEBUG
  for(loint i = 0; i < destination_size; i++) {
    assert(destination[i] != '\n');
  }
#endif
}

//! Sets the name of the given column
char *setName(const bool first_column,Parse &p, char *pos_start, char *&pos_seperator){
//char*setName(const bool first_column,Parse &p, char *&pos_start, char *&pos_seperator){
  if(pos_start) {
    if(*pos_start == '\n') pos_start++;
    const uint length_name = pos_seperator - pos_start;
#ifndef NDEBUG
    //! Asserts that no '\n' is included in the string sent:
    for(loint i = 0; i < length_name; i++) {
      assert(pos_start[i] != '\n');
    }
#endif
    if(first_column) {
      p.set_name_in(pos_start, length_name);
      assert(0 == strcmp(p.get_name_in(), pos_start));
    } else {
      p.set_name_out(pos_start, length_name);
      assert(0 == strcmp(p.get_name_out(), pos_start));
    }
    pos_seperator++; // jumps over the seperator
  }
  return pos_seperator;
}

//! Sets the name given the input
void setName(char *start_pos, char *&name, uint index, uint *array, uint array_size) {
  char *name_start_pos = start_pos;
#ifndef NDEBUG
  if(!((index == 0) || (array_size >= 2))) {
    printf("index(%u) and array_size(%u)\n", index, array_size);
  }
  assert((index == 0) || (array_size >= 2)); 
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
  if(legal_format) {
    // Rellocates if it's not long enough:
    if(name && !(name_length < strlen(name))) {delete [] name; name = new char[name_length+1];}
    else if(!name) {name = new char[name_length+1];}
    name[name_length] = '\0'; // Sets the end of the string.
    strncpy(name, name_start_pos, name_length);
  } else {name = NULL;}

#ifndef NDEBUG
  if(name == NULL) {
    fprintf(stdout, "!!\tName not set for index(%u) ", index);
    fprintf(stdout, "!!\t at line %d in %s: ",  __LINE__, __FILE__);
    for(uint i = 0; i < array_size; i++) printf("[%u]=%u\t", i, array[i]); printf("\n");
  }
#endif
}


//! Sets the taxon of the given column
char *setTaxon(const bool first_column,Parse &p,char *&pos_start,char*&pos_seperator) {
  if(first_column) {
    const uint length_name = pos_seperator - pos_start; 
    p.set_taxon_in(pos_start, length_name);
  } else {
    const uint length_name = pos_seperator - pos_start; 
    p.set_taxon_out(pos_start, length_name);
  }
  pos_seperator++; // jumps over the seperator
  while (*pos_seperator == ' ') pos_seperator++;
  return pos_seperator;
}

//! Returns the postion at the end of column
char *getTaxonPartEnd(char *start_pos) {
  char *pos_column_end = strchr(start_pos, '\t');    
  if (pos_column_end == NULL) {
    pos_column_end = strchr(start_pos, ' ');
    if(pos_column_end != NULL) while (*pos_column_end == ' ') pos_column_end++;
  }
  return pos_column_end;
}

//! Returns the postfix of the name
char *getNamePartEnd(char *start_pos, char seperator) {
  char *pos_seperator = strchr(start_pos, seperator); // position of the first position of seperator
  return pos_seperator;
}

// Validates that an enalrged array is correctly enlarged
void assert_enlarge_array(uint *old_arr, uint old_size, uint *new_arr, uint new_size, uint default_value) {
  if(false && (old_arr) && (new_arr) && (new_size) && (default_value)) { ;} // Do not print warnings is 'assert' is deactivated
  for(uint i = 0; i< old_size; i++) assert(old_arr[i] == new_arr[i]);
  for(uint i = old_size; i< new_size; i++) assert(new_arr[i] == default_value);
}


/**
   @brief Gets the labels from a line in the blast file.
   @param <first_column> True if it's the leftmost protein in the given blast-file line
   @param <p> The Parse object to set the names into.
   @param <pos_column_start> The memory position to start the work from.
   @param <logical_end> The last address where data may be find; used to avoid errors.
   @param <b> The tsettings_input_t object iot parse correctly.
   @return The position after the data were the given column were located.
**/
char *blast_extractors::getIDColumn(const bool first_column, Parse &p, char *pos_column_start, char *logical_end, tsettings_input_t b) {
  char *temp = p.set_names(first_column, pos_column_start, b.INDEX_IN_FILE_FOR_PROTEIN_NAME, b.INDEX_IN_FILE_FOR_TAXON_NAME, b.SEPERATOR, logical_end, b.DEFAULT_NUMBER_OF_COLUMNS_IN_NAME);
  return temp;
}

//! For backward compability: parses straing forward, without the suage of indexes to structure the string
char *static_getIDColumn(const bool first_column, Parse &p, char *&pos_column_start, char *logical_end, char seperator) {
  char *pos_seperator = getNamePartEnd(pos_column_start, seperator);
  if (pos_seperator != NULL) {// The name
    pos_seperator = setName(first_column,  p, pos_column_start, pos_seperator);
    char *pos_column_end = getTaxonPartEnd(pos_seperator); // returns the postion at the end of column
    if(pos_column_end != NULL) {
      pos_column_end = setTaxon(first_column, p, pos_seperator, pos_column_end); //! sets the inner taxon
      return pos_column_end;
    }
  } 
  return logical_end;
}


/**
   @brief Gets the labels for a protein pair, using a line in the blast file.
   @param <p> The Parse object to set the names into.
   @param <pos_column_start> The memory position to start the work from.
   @param <logical_end> The last address where data may be find: used to avoid errors.
   @param <blast_settings> The tsettings_input_t object iot parse correctly.
   @return True if data were read correctly.
**/
bool blast_extractors::getParseHeaders(Parse &p, char *&pos_column_start, char *logical_end, tsettings_input_t blast_settings) {
  char *temp = getIDColumn(true, p, pos_column_start, logical_end, blast_settings);
  if(temp != NULL) {
    pos_column_start = temp; // updates the position of the string
    while((*pos_column_start == '\t') ||  (*pos_column_start == ' ')) pos_column_start++; // jump over white spaces
    temp= getIDColumn(false, p, pos_column_start, logical_end, blast_settings);
    if(temp != NULL) {
      pos_column_start = temp;
      return true;
    }
#ifndef NDEBUG
    else  printf("Did not find the 2-second id-column (label) at line %d in %s\n", __LINE__, __FILE__);
#endif
  } 
#ifndef NDEBUG
  else  printf("Did not find the 1-first id-column (label) at line %d in %s\n", __LINE__, __FILE__);
#endif
  return false;
}

