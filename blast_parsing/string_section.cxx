#include "string_section.h"

  // -----------------------------------------------------------
  // Standard proecedures for getting the results:

//! Returns the length of the string
mem_loc string_section::getStringLength() {
  return (buffer_main_logical_end - buffer_main_start);
}

/**
   @brief An assert function returning true if input equals
   @return true if input equals the content of this object's 'buffer_main_start' variable.
**/
bool string_section::is_equal(char *inp, long int inp_length) {
#ifndef NDEBUG
  for(uint i = 0; i < inp_length; i++) {
    assert(inp[i] == buffer_main_start[i]);
  }
  return true;
#else 
  if(inp && inp_length) {;}
  return true;
#endif
}
 
//! Returns the number of chars containing data.
long int string_section::get_data_length() {
  if(buffer_main_start != NULL) return (buffer_main_logical_end - buffer_main_start);
  return 0;
}

//! Prints the data block given the argument.
void string_section::print_data_block(long int end_pos) {
  if(is_set()) {
    string_section::print_segment(buffer_main_start, buffer_main_start + end_pos);
  } else fprintf(stderr, "!!\tData block not set in %s at line %d.\n", __FILE__, __LINE__);
  // else fprintf(stderr, "!!\tData block not set.\n");
}

//! Prints the data block
void string_section::print_data_block() {
  if(is_set()) {
    string_section::print_segment(buffer_main_start, buffer_main_logical_end);
  } else fprintf(stderr, "!!\tData block not set in %s at line %d.\n", __FILE__, __LINE__);
}

//! Prints the segment given the inputs
void string_section::print_segment(char *m_start, char *end) {
  char *start =m_start;
  //  printf("[\n");
  //while(start != end) printf("%c", *start++); //putchar(*start++);
  //printf("length=%d in string_section.cxx\n", end-m_start);
  while(start != end) {
    if(*start != '\n') putchar(*start++);
    else {
      if(true) printf("\n");
      else printf("\t<newline>\n");
      start++;
    }
  }
    //  printf("\n]\n");
}

//! @return true if object contains data.
bool string_section::is_set() {
  if(buffer_main_logical_end != NULL && buffer_main_start != NULL) {
    assert(buffer_main_start <= buffer_main_logical_end);
    return true;
  } else return false;
}

/**
   @return the file position in the file where the next new line in the blast file starts.
   @remarks Return 0 if a shift from one type of left-most protein to another is not found.
**/
loint string_section::get_last_line_pos(char *buffer_start, char *last_block_line) {
  loint block_length = 0;
  if(last_block_line != NULL) {
    block_length = (last_block_line - buffer_start);
    if(block_length > 0) block_length--;
  }
  long long int last_line_pos = block_length;
  if(last_line_pos < 1) { // In the middle of a data-block
    last_line_pos = 0;
  } else if(!end_of_file) {
    //    last_line_pos+= prev_index; // the position of the file
    last_line_pos++; // Updates the position in the file
  } else {
    //    last_line_pos = prev_index + getStringLength();
    last_line_pos = getStringLength();
  }
  return last_line_pos;
}

/**
   @brief Tests if block has a shift from one type of left-most protein to another. 
   @return A value greater than 0 if a position is found.
**/  
loint string_section::more_than_one_leftmost_protein(taxon_list_t *listProteins, tsettings_input_t blast_settings) {    
  char* last_block_line = 0;
  char *buffer_one = buffer_main_start, *logical_end = buffer_main_logical_end;
  char *buffer_start = buffer_one;
  uint cnt_liens = 0;
  buffer_one = strchr(buffer_one, '\n'); // As the input might start in the moddle of a row, this ensures that it's moved to the end of the first line.
  // Jumps to the first newline in the file:
  buffer_one = strchr(buffer_one, '\n')+1;       // update, jumping over the line end:
  /**
     A problem when processing data, without knowing the contect of it,
     is that the data before is unknown, ie, a block of siilar pairs may
     be at the lines before the data of this starts, causing dependencies we are not
     interested in; to avoid this, we require that a second block is found, and by
     this ensuring independency: <-- 16.07.2012 by oekseth
  **/
  uint cnt_block_starts = 0; 
  while (buffer_one && buffer_one < logical_end) {
    Parse p = Parse();
    char *line_start = buffer_one;
    char *temp = blast_extractors::getIDColumn(true, p, buffer_one, logical_end, blast_settings);
    if(temp != NULL) {
      if(listProteins->insert_protein(p.get_taxon_in(), p.get_name_in())) {
	if(cnt_block_starts==1) {
	  last_block_line = line_start; 	// Sets the end of the last block before a new protein is read. This in order to avoid overlaps in the aprsing of the whole file, and thereby redcusing the mem consumption	
	  buffer_one = NULL; // TODO: Write some tests ensuring the correctness of this.
	} else {cnt_block_starts++, buffer_one = strchr(temp, '\n')+1;}       // update, jumping over the line end:
      } else buffer_one = strchr(temp, '\n')+1;       // update, jumping over the line end:
    } else {
      buffer_one = strchr(buffer_one, '\n')+1;       // update, jumping over the line end:
    }
    p.free_memory();
    cnt_liens++;
    if(buffer_one == NULL) buffer_one = logical_end;
  }
  //
  // Calculates the block length:
  const loint last_line_pos = get_last_line_pos(buffer_start, last_block_line);
#ifndef NDEBUG
  //! Verifies that our expectations hold:
  if(last_block_line > 0) {
    if(buffer_start[last_line_pos-1] != '\n') {      
      assert(buffer_start[last_line_pos-1] == '\n');
    }
  }
#endif
  return last_line_pos; 
}  

//! Deallocates the memory
void string_section::finalize() {
  if(buffer_main_start) {
    delete [] buffer_main_start, buffer_main_start = NULL;
  }
}


//! @return An object of this given the input params.
string_section *string_section::init(char *buffer_main_start,
				     char *buffer_main_logical_end,
				     long int chars_sendt_to_parsing,
				     int block_cnt, bool end_of_file
				     ) {
  return new string_section(buffer_main_start, buffer_main_logical_end,
			    chars_sendt_to_parsing, block_cnt, end_of_file);  
}


//! Deallocates the memory given the param.
void string_section::close(string_section *&obj) {
  if(obj) {obj->finalize(), delete obj, obj = NULL;}
}


//! The constructor.
string_section::string_section() :
  buffer_main_start(NULL),
  buffer_main_logical_end(NULL), 
  prev_index(0),
  block_cnt(0),
  end_of_file(false)
{}

//! The constructor.
string_section::string_section(char *_buffer_main_start,
			       char *_buffer_main_logical_end,
			       long int _prev_index, int _block_cnt, bool _end_of_file) : 
  buffer_main_start(_buffer_main_start),
  buffer_main_logical_end(_buffer_main_logical_end),
  prev_index(_prev_index),
  block_cnt(_block_cnt),
  end_of_file(_end_of_file)
{}

//! The main test function for this class  
void string_section::assert_class(const bool print_info) {
  const static char *class_name = "string_section";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  // Initial tests:
  char *buff = new char[100];//(char*)malloc(sizeof(char)*100);
  memset(buff, '-', 100);
  string_section temp(buff, buff+100, 0, 0, true);
  assert(temp.is_set());
  temp.finalize();
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

