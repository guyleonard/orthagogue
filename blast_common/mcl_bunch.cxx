#include "mcl_bunch.h"

/**!
   @Name: set_line_end(..) -- Called when no more data is to be added to the current line. Adds line end if data set
   @Changed: 21.01.2011 by oekseth, setting the doaalr sign as an option
*/
uint mcl_bunch::set_line_end(const bool use_dollar_sign, const bool use_test) { 
  uint cnt_elements_inserted = 0;
  if(use_test) {
    if(position_of_data_start > position_of_line_start){
      if(use_dollar_sign) {
	set_char('$');
	cnt_elements_inserted++;
      }
      set_char('\n');
      cnt_elements_inserted++;
      position_of_line_start = current_size_string;
      position_of_data_start = current_size_string;
    }// else printf("(line end not set)\t");
  } else if(position_of_data_start == position_of_line_start){
    if(use_dollar_sign) {
      set_char('$');
      cnt_elements_inserted++;
    }
    set_char('\n');
    cnt_elements_inserted++;
    position_of_line_start = current_size_string;
    position_of_data_start = current_size_string;
  } ;//else set_char('?');
  return cnt_elements_inserted;
}


//! Inserts the name: Returns true if name not NULL
bool mcl_bunch::set_line_start(char *name) {
  if(name != NULL) {
    enlarge(); // a standard call
    const uint size_name = strlen(name);
    strncpy(string + current_size_string, name, size_name);
    current_size_string += size_name;
    if(size_name) assert((int)size_name == (int)(get_size_of_header(name)-1));
    return true;    
  }
  return false;
}

//! Inserts the name as an integer.
void mcl_bunch::set_line_start(uint world_index_in) {
  enlarge(); // a standard call
  uint strl_len = 1;
  if(world_index_in > 0) strl_len = (int)log10(world_index_in)+1;
  assert((strl_len+1) == get_size_of_header(world_index_in));
  
  sprintf(string + current_size_string, "%u", world_index_in);
  current_size_string += strl_len;
}

//! Enlarges the char array, and updates the size variabel
void mcl_bunch::enlarge() {
  uint new_size = current_size_string + 1024;
  if(new_size > size_string) { // the array has to be enlarged
    new_size*=2; // extens it
    char *temp = allocate_string_buffer(new_size); //allocates
    assert(temp);
    if(!size_string) {
      memset(temp,MCL_SEPERATOR, new_size); 
    } else {
      memcpy(temp, string, current_size_string);
      free_string_buffer(string, size_string,			 
			 //reset_internal_poitner=
			 false); // deallocated
    }
    string = temp; size_string = new_size; // updates the size
    memset(string + current_size_string, MCL_SEPERATOR, size_string - current_size_string);
  } 

}

/**
   @Name: enlarged(..) -- Enlarges the string
   @Usage: Called when a holw bunch of chars is copied from another class to this
   @Changed: 07.01.2011 by oekseth
*/
void mcl_bunch::enlarge(uint size_to_be_inserted) { 
  uint new_size = current_size_string + size_to_be_inserted;
  if(new_size > size_string) { // the array has to be enlarged
    new_size*=2; // extens it
    char *temp = allocate_string_buffer(new_size); //allocates
    assert(temp); // Aborts it not set.
    //    char *temp = (char*)tbb::tbb_allocator<char>().allocate(sizeof(char)*new_size);
    memcpy(temp, string, current_size_string);
    free_string_buffer(string, size_string, /*reset_internal_poitner=*/false); // deallocated
    //    tbb::tbb_allocator<char>().deallocate(string, size_string);
    string = temp; size_string = new_size; // updates the size
    memset(string + current_size_string, MCL_SEPERATOR, size_string - current_size_string);
  }
}

//! @return the total number of chars used rerpesenting the header.
uint mcl_bunch::get_size_of_header(char *label) {
  if(label) {
    return strlen(label) + 1; // The length of the name + a seperator char.
  }
  return 0;
}

//! @return the total number of chars used rerpesenting the header.
uint mcl_bunch::get_size_of_header(uint index) {
  if(index > 0) return ((uint)log10(index) + 2); // The length of the name + a seperator char.
  else return 2;
}
//! @return the size of the sim score:
uint get_size_of_sim_score(float sim_score) {
  uint sim_score_length = 1;
  if(sim_score > 1.0) {sim_score_length = (int)log10(sim_score)+1;}
  sim_score_length += 1; // For the '.' seperator between the integers- and decimals.
  sim_score_length += 3; // For the precision we use as default basis
  return sim_score_length;
}

//! @return the length of pair to be inserted.
uint mcl_bunch::get_size_of_inserted_pair(uint world_index_in, uint world_index_out, float sim_score) {
  uint length_of_string = 1;
  if(world_index_in > 0) length_of_string   = (uint)log10(world_index_in)+1;
  if(world_index_out > 0) length_of_string += (uint)log10(world_index_out)+1;
    
  const uint total_length = length_of_string + 3 + get_size_of_sim_score(sim_score); // '+3' due to 3 seperations between the chars.
  return total_length;
}


//! @return the length of pair to be inserted.
uint mcl_bunch::get_size_of_inserted_pair(uint world_index_out, float sim_score) {
  uint length_of_string = 1;
  if(world_index_out > 0) length_of_string = (uint)log10(world_index_out)+1;    
  const uint total_length = length_of_string + get_size_of_sim_score(sim_score)+2; // '+2' for the seperators
  return total_length;
}

/**
   @return the length of pair to be inserted.
   @remarks
   - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
   - Makes is possible using more effective writes when MPI is activated.
**/
uint mcl_bunch::get_size_of_inserted_pair(char *name_in, char *name_out, float sim_score) {
  const uint size_names = strlen(name_in) +1 +strlen(name_out) +1; // The names + 2 seperators
  const uint size_similarity = get_size_of_sim_score(sim_score);
  const uint total = size_names + size_similarity +1; // '+1' for the newline
  return total;
}

  /**
     @return the length of an mcl-row-element
     @remarks
     - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
     - Makes is possible using more effective writes when MPI is activated.
  **/
uint mcl_bunch::get_size_of_inserted_pair(char *name_out, float sim_score) {
  const uint size_names = strlen(name_out) +1; // The name + 1 seperator
  const uint size_similarity = get_size_of_sim_score(sim_score);
  const uint total = size_names + size_similarity +2; // '+1' for the newline
  return total;
}

//! Sets the header
void mcl_bunch::set_header(uint world_index) {
#ifndef NDEBUG
  const uint pos_of_start = current_size_string;
  #endif
  position_of_line_start = current_size_string;
  set_line_start(world_index);
  string[current_size_string] = SEPERATOR_AFTER_ROW_START;
  current_size_string++;
  position_of_data_start = current_size_string;
#ifndef NDEBUG
  //! Validates our expectations regarding the inserted string:
  const uint pos_of_end = current_size_string;
  const uint chars_added = pos_of_end - pos_of_start;
  const uint expected_size = get_size_of_header(world_index);
  assert(chars_added == expected_size);
  #endif
}

/**! Sets the header*/
void mcl_bunch::set_header(char *name) {
#ifndef NDEBUG
  const uint pos_of_start = current_size_string;
#endif
  position_of_line_start = current_size_string;
  set_line_start(name);
  string[current_size_string] = SEPERATOR_AFTER_ROW_START;
  current_size_string++;
  position_of_data_start = current_size_string;
#ifndef NDEBUG
  //! Validates our expectations regarding the inserted string:
  const uint pos_of_end = current_size_string;
  const uint chars_added = pos_of_end - pos_of_start;
  const uint expected_size = get_size_of_header(name);
  assert(chars_added == expected_size);
#endif
}

/**! Prints the current status of the string; useful for dbugging purposes,
   getting a correspondence between expactaitns and known facts about whats inserted.   */
void mcl_bunch::print_data() {
  printf("\tcurrent_size_string(%u), position_of_line_start(%u), position_of_data_start(%u)\n",
	 current_size_string
	 ,position_of_line_start
	 ,position_of_data_start);
}

//! Inserts a pairwise relation using integer data
void mcl_bunch::insert(uint world_index_in, uint world_index_out, float sim_score, float div_factor) {
  //  if(div_factor>0 && sim_score > 0)
  assert(div_factor>0);
  assert(sim_score>0);
{
#ifndef NDEBUG
  const uint pos_of_start = current_size_string;
#endif
    enlarge();
    uint str_len = 1;
    if(world_index_in > 0) str_len = (uint)log10(world_index_in)+1;
    if(world_index_out > 0) str_len += (uint)log10(world_index_out)+1;
    else str_len +=1;
    //! Calcualting the original length, to be used setting the precision, using a minimum of 3 decimals:
    uint original_length = 1;
    if(sim_score > 1.0) {original_length = (int)log10(sim_score)+1;}
    //! Calculates the length of the normalised score:
    const float similarity = sim_score/div_factor;
    uint sim_len = 1;
    if(similarity > 1.0) {sim_len = (int)log10(similarity)+1;}
    //! Increases the precision, making the length of the file independent 
    uint precision = 3;
    if(original_length != sim_len) {
      assert(sim_len < original_length); // Our assumption.
      const uint difference = original_length - sim_len;
      precision += difference;
      sim_len += difference;
    }
    sprintf(string + current_size_string, "%u%c%u%c%-.*f", world_index_in, MCL_SEPERATOR, world_index_out, MCL_SEPERATOR, precision, similarity); // Add the distance
    
    //    sprintf(string + current_size_string, "%u%c%u%c%-.3f", world_index_in, MCL_SEPERATOR, world_index_out, MCL_SEPERATOR, similarity);
#ifndef NDEBUG
  //! Validates our expectations regarding the inserted string:
  const uint pos_of_end = current_size_string;
  const uint chars_added = pos_of_end - pos_of_start;
  const uint expected_size = get_size_of_inserted_pair(world_index_in, world_index_out, sim_score);
  assert(chars_added == expected_size);
#endif
    //! Add the distance
    assert((str_len + sim_len + 6) == get_size_of_inserted_pair(world_index_in, world_index_out, sim_score));
    current_size_string += str_len + sim_len + 6; // '+6' = '+3' (due to 3 seperations between the chars) and '+3' (due to precision)
    set_char('$'); set_char('\n'); position_of_line_start = current_size_string;
  }
}

//! Inserts a pairwise relation using names
void mcl_bunch::insert(char *name_in, char *name_out, float sim_score, float div_factor) {
  assert(sim_score > 0);
  assert(div_factor > 0);
  //  if(div_factor>0      //&& sim_score > 0      ) 
  {
    //#ifndef NDEBUG
      const uint pos_of_start = current_size_string;
      //#endif
    const bool ret_val_in = set_line_start(name_in);
    assert(ret_val_in);
    {
      current_size_string++; //  creates a space
      const bool ret_val_out = set_line_start(name_out);
      assert(ret_val_out);
      {
	current_size_string++; //  creates a space
	insert_sim_score(sim_score, div_factor);

	set_char('\n'); position_of_line_start = current_size_string;
	//#ifndef NDEBUG
	//! Validates our expectations regarding the inserted string:
	const uint pos_of_end = current_size_string;
	const uint chars_added = pos_of_end - pos_of_start;
	const uint expected_size = get_size_of_inserted_pair(name_in, name_out, sim_score);
	assert(chars_added == expected_size);
	//#endif
      }
    }
  }
#ifndef NDEBUG
  if(has_sub_string(name_out) == false) {
    printf("(failed)\t\t inserts %s->%s not inserted, at mcl_bunch:%d\n", name_in, name_out, __LINE__);
    assert(false);
  }
#endif
}


//! Usage:  ...   const float div_factor_out = arrAvgNorm[taxon_in][taxon_out];
void mcl_bunch::insert(uint world_index_in, float sim_score, float div_factor) {
  assert(div_factor > 0);
  assert(sim_score > 0);
  if(div_factor>0 && sim_score > 0) {
#ifndef NDEBUG
      const uint pos_of_start = current_size_string;
#endif
    enlarge();
    uint str_len = 1;
    if(world_index_in > 0) str_len = (uint)log10(world_index_in)+1;
    //! Calcualting the original length, to be used setting the precision, using a minimum of 3 decimals:
    uint original_length = 1;
    if(sim_score > 1.0) {original_length = (int)log10(sim_score)+1;}
    //! Calculates the length of the normalised score:
    const float similarity = sim_score/div_factor;
    uint sim_len = 1;
    if(similarity > 1.0) sim_len = (int)log10(similarity)+1;
    //! Increases the precision, making the length of the file independent 
    uint precision = 3;
    if(original_length != sim_len) {
      assert(sim_len < original_length); // Our assumption.
      const uint difference = original_length - sim_len;
      precision += difference;
      sim_len += difference;
    }
    sprintf(string + current_size_string, "%u:%-.*f%c", world_index_in, precision, similarity, SEPERATOR_IN_ROWS); // Add the distance
    current_size_string += str_len + sim_len + 6;    
#ifndef NDEBUG
    //! Validates our expectations regarding the inserted string:
    const uint pos_of_end = current_size_string;
    const uint chars_added = pos_of_end - pos_of_start;
    const uint expected_size = get_size_of_inserted_pair(world_index_in, sim_score);
    assert(chars_added == expected_size);
#endif
  }
}


// Uses input of type arrKey[index_out], str_len);
void mcl_bunch::insert(char *name, float sim_score, float div_factor) {

  if(div_factor>0 && sim_score > 0) {
    if(set_line_start(name)) {
      //      const uint pos_of_start = current_size_string;
      string[current_size_string++] = ':';
      insert_sim_score(sim_score, div_factor);
      string[current_size_string++] = SEPERATOR_IN_ROWS;
    }
  }
}

/**! Returns the size of the line. Needed to store this string of data   */
uint mcl_bunch::get_size_of_line() {    
  const long int diff = (current_size_string - position_of_line_start);    
  et("_get_size_");
  assert(diff > -1);
  return (uint)diff;
}

//! Asserts that the current string resides inside the line
void mcl_bunch::et(char *s) {    
  assert(!(current_size_string < position_of_line_start));
  if(false && current_size_string < position_of_line_start) {
      
    printf("%s:\tcurrent_size_string(%u) - position_of_line_start(%u)\n", s, current_size_string, position_of_line_start);
    assert(false);
  }
}

//! Initally written for visual asseritng the input with regard to changing the format
void mcl_bunch::validate() {
  printf("Iterates over space <%u, %u]\n", 0, current_size_string);
  for(uint i =0; i<current_size_string; i++) {
    if(string[i] == 0) printf("[%u]  = '' (size_string = %u)\n", i, size_string);
    else if(string[i] == ' ') printf("string[%u]=' '\n", i);
    else printf("%c\n", string[i]);
  }
}

//! Inserts a char at the end of the data block in memory
void mcl_bunch::set_char(char c) {
  if(string != NULL) {
    if(!(current_size_string < (size_string-1))) enlarge();
    assert((current_size_string < (size_string-1)));
    string[current_size_string++] = c;
  } else fprintf(stderr, "An error occured: tried writing data to an empty string!\n");
}

void mcl_bunch::insert_sim_score(float sim, float div_factor) {
  assert(sim > 0);
  assert(div_factor > 0);
#ifndef NDEBUG
  const uint pos_of_start = current_size_string;
#endif
  //  if(sim > 0      && (div_factor > 0)) 
{
    //! Calcualting the original length, to be used setting the precision, using a minimum of 3 decimals:
    uint original_length = 1;
    if(sim > 1.0) {original_length = (int)log10(sim)+1;}
    //! Calculates the length of the normalised score:
    const float similarity = sim/div_factor;
    uint sim_len = 1;
    if(similarity > 1.0) sim_len = (uint)log10(similarity)+1;

    //! Increases the precision, making the length of the file independent 
    uint precision = 3;
    if(original_length != sim_len) {
      assert(sim_len < original_length); // Our assumption.
      const uint difference = original_length - sim_len;
      precision += difference;
    }
    //    char *old_string = string;
    sprintf(string + current_size_string, "%-.*f", precision, similarity); // Add the distance
    assert(string[current_size_string+sim_len + precision] != '\0');
    current_size_string += sim_len + 1 + precision; 
#ifndef NDEBUG
    //! Validates our expectations regarding the inserted string:
    const uint pos_of_end = current_size_string;
    const uint chars_added = pos_of_end - pos_of_start;
    const uint expected_size = get_size_of_sim_score(sim);
    assert(chars_added == expected_size);
#endif 
  }
}

/**! Sets the name_index-thing
   @Changed: 10.02.2011 by oekseth: Added a new default, printing the protein names in 'index=> label", accroding to adjusted spec from Dr . Mironov
*/
uint mcl_bunch::set_name_index(char *name, uint world_index_in) {
  if(name != NULL) {
    enlarge();
    uint str_len = 1;
    if(world_index_in > 0) str_len = (uint)log10(world_index_in)+1;
    str_len += strlen(name);
    if(true) {
      sprintf(string + current_size_string, "%u%c%s", world_index_in, SEPERATOR_IN_ROWS, name); // Add the distance
    } else {
      sprintf(string + current_size_string, "%s%c%u", name, SEPERATOR_IN_ROWS, world_index_in); // Add the distance
    }
    current_size_string += str_len +  1;
    position_of_data_start = current_size_string;
    position_of_line_start = current_size_string;
    const uint cnt_line_end = set_line_end(false, false);
    const uint inserted_cnt = str_len+cnt_line_end +1; // '+1' due to the seperator.
    assert(inserted_cnt == get_size_name_index(name, world_index_in));
    return inserted_cnt;
  }
  return 0;
}

//! @return the number of chars used representing the given string
uint mcl_bunch::get_size_name_index(char *name, uint world_index_in) {
  if(name) {
    uint string_length = 1;
    if(world_index_in > 0) string_length = (uint)log10(world_index_in)+1;
    string_length += strlen(name);
    string_length += 2; // Due to seperator (a) between the name and index and (2) a newline ('\n').
    return string_length;
  }
  return 0;
}
// Prints the data in the container
void mcl_bunch::print() {
  if(size_string > 0 && string != NULL) {
    printf("|");
    for(uint i=0; i < current_size_string; i++)
      putchar(string[i]);
    printf("|\n");
  } else printf("(empty)\n");
}


/**
   @Name: copy_line(..) -- Copies the lines from the arguments, assuming that only the last argument
   has a line end included.
   @Assumption: If arg_2 is empty, then arg_3 is empty as well.
   @Changed: 07.01.2011 by oekseth
*/
void mcl_bunch::copy_line(mcl_bunch &arg_1, mcl_bunch &arg_2, mcl_bunch &arg_3) {
  // TODO: Instead of removing the 'end_line' append the 'end_liene' after the line thing instead!
  const uint arg_1_size = arg_1.get_size_of_line(); 
  if(arg_1_size > 0) {
    const uint arg_2_size = arg_2.get_size_of_data();
    const uint arg_3_size =      arg_3.get_size_of_data();
    enlarge(arg_1_size + arg_2_size + arg_3_size); // ensures that enough memory is reserved
    // First argument:
    strncpy(string + current_size_string, arg_1.getLine(), arg_1_size);
    current_size_string += arg_1_size;
    if(arg_2_size > 0)  {
      uint offset_2 = 2; // the line ends at the two chars '$' and '\n', thereby an offset is needed
      char *data_2 = arg_2.getLineData();
      if(data_2[arg_2_size-1] != '\n')	offset_2 = 0;
      strncpy(string + current_size_string, data_2, arg_2_size- offset_2);
      current_size_string += arg_2_size - offset_2;
    }
    if(arg_3_size > 0) {
      uint offset_3 = 2; // the line ends at the two chars '$' and '\n', thereby an offset is needed
      char *data_3 = arg_3.getLineData();
      if(data_3[arg_3_size-1] != '\n')	offset_3 = 0;
      strncpy(string + current_size_string, data_3, arg_3_size - offset_3);
      current_size_string += arg_3_size-offset_3;
    }
    position_of_line_start = current_size_string;
    position_of_data_start = current_size_string;
  }
}


#ifdef USE_MPI
//! Writes the string to the file, and clears the memory
uint mcl_bunch::write_data_to_file(MPI_File  *file_out)
#else
uint mcl_bunch::write_data_to_file(FILE  *file_out)
#endif
 {
  set_end_of_block();
  if(size_string > 0 && string != NULL) {
    const int ret_val = mcl_bunch::write_string_to_file(file_out, string); 
    //    const int ret_val = fputs(string, file_out);
    if(ret_val == EOF) {
      fprintf(stderr, "!!\tDisk quota (most likely) exceeded, as return_value from writing %u elements to the file was '%d': Finishes the operation (ie, do not abort), but the result file will no be complete. Please either (a) change result path to a new location with bigger disk quota or remove files on your system, before you (b) rerun this software. If neither of these alternatives help, please contact the developer at his email oekseth@gmail.com . This error was generated at line %d in %s fund in method %s\n", size_string, ret_val, __LINE__, __FILE__, __FUNCTION__);
    }
    return (uint)ret_val;
  }
  return 0;
}

#ifdef USE_MPI
/**
   @brief Uses MPI-lib for the write operation:
   @remarks Is the interace for writing all the data to the result files (in the final phase of the ortAgogue algorithm).
**/
int mcl_bunch::write_string_to_file(MPI_File *file, char *string, uint length) {
  assert(file);
  assert(string);
  assert(length);
  MPI_Status status_return;
  MPI_File_write_shared(*file, string, length, MPI_CHAR, &status_return);
  //  MPI_File_iwrite_at() 
  // MPI_File_sync
  //  MPI_File_write(*file, string, length, MPI_CHAR, &status_return);
  //! A goal is to provide error information, for example when disk buffer is exceeded, with the concequence that files are not written.
  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
  //MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes);
  int count = 0; MPI_Get_count(&status_return, MPI_CHAR, &count);
  if(count != (int)length) {
    printf("!!\t[myrank=%d]\tGot count(%d) for data_length(%u), at line %d in file %s\n", myrank, count, length, __LINE__, __FILE__); 
    return EOF; 
  } else return length;
}

#else 

/**
   @brief Uses c-lib for the write operation:
   @remarks Is the interace for writing all the data to the result files (in the final phase of the ortAgogue algorithm).
**/
int mcl_bunch::write_string_to_file(FILE *file, char *string, uint current_size_string) {
  assert(file);
  if(false && current_size_string) ; // In order to hide the variable, making the itnerface generic:
  if(string) {
    static bool error_is_printed = false;
    //  fprintf(file, string);                                                                                                                                                    
    if(true) {
      const uint cnt_wrote = fwrite(string, 1, strlen(string), file);
      if(cnt_wrote != strlen(string)) {
	if(!error_is_printed) {
          fprintf(stderr,
                  "!!\t We might have exceeded the disk-size of the file system:\n"
		    "-\t either increase the available space of your disk, or use a more restricte filter when initiating orthAgogue.\n"
                  "-\t we infer this from your effort of writing \"%u\" chars, while only \"%u\" chars were written.\n"
		  "-\t To avoid cluttering of the output, we only print this message once for each core.\n"
                  "-\t To verify that errors is not \"fiction\", please compare \"%s\" with the output-files.\n"
                  "-\t If this error is not understood, please contact the developer at [oekseth@gmail.com].\n"
                  "Error generated at [%s]:%s:%d\n",  (uint)strlen(string), cnt_wrote,  string, __FUNCTION__, __FILE__, __LINE__);
          error_is_printed = true;
	}
        assert(cnt_wrote == strlen(string));
      }
      return (int)cnt_wrote;
    } else {
      return fputs(string, file);
    }
  } else return 0;
  /*
  if(string) {
    return write_string_to_file(file, string, (uint)strlen(string));
  } else return 0;
*/
}
#endif
void mcl_bunch::free_string_buffer(char *&string, uint &size_string, /*reset_internal_poitner*/ bool reset_internal_poitner) {
  delete [] string;    
  // TODO: Valgrind complains when using the tbb-internal-things: Why?
  //    tbb::tbb_allocator<char>().deallocate(string, size_string);
  string = NULL, size_string = 0;
  if(reset_internal_poitner) current_size_string = 0, position_of_data_start = 0, position_of_line_start = 0;
}

char *mcl_bunch::allocate_string_buffer(long long int new_size) {
  assert(new_size);
  char *new_string = new char[new_size+1];
  log_builder::test_memory_condition_and_if_not_abort(new_string!=NULL, __LINE__, __FILE__, __FUNCTION__);
  new_string[new_size] = '\0';
  return new_string;
  /*
  char *new_string = new char[new_size];
  log_builder::test_memory_condition_and_if_not_abort(new_string!=NULL, __LINE__, __FILE__, __FUNCTION__);    
  return new_string;
  */
  // TODO: Valgrind complains when using the tbb-internal-things: Why?
  //    return (char*)tbb::tbb_allocator<char>().allocate(sizeof(char)*new_size);
}

void mcl_bunch::free_mem() {
  if(string) {
    free_string_buffer(string, size_string, true);
  }
}
void mcl_bunch::free_class(mcl_bunch *&arg) {
  if(arg != NULL) {
    arg->free_mem();
    delete arg; arg = NULL;
  }
}

mcl_bunch::mcl_bunch() : size_string(0), current_size_string(0), position_of_line_start(0), position_of_data_start(0), 
			 string(NULL)
{enlarge();
  et("mcl_bunch_constructor_");
}



#ifdef assert_code
//! Compare the internal-object-string with argument given
bool mcl_bunch::compare_string(char *inp, const bool print_diff) {
  if(current_size_string != 0) {
    bool is_equal = true;
    uint size = current_size_string;
    if(strlen(inp) < current_size_string) 	size = strlen(inp);
    for(uint i =0; i<size; i++) {
      if(string[i] != inp[i]) {
	is_equal = false;
	if(print_diff) printf("'%c'(%d) != '%c' (%d)\n", string[i], string[i], inp[i], inp[i]);
      }
    }
    return is_equal;      
  } else if(inp == NULL || strlen(inp) == 0) {
    //if(print_diff) printf("input=NULL\n");
    return true;
  }
  else if(print_diff) printf("current_size_string(%u) not set.\n", current_size_string);
  if(print_diff) printf("!!\tstring not set\n");
  return false;
}

#endif
/**
   @Name: assert_class(..) -- Validates the code: The main validation of this code is exectuted in clas "mcl_format".
*/
void mcl_bunch::assert_class(const bool print_info) {
  const static char *class_name = "mcl_bunch";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  class mcl_bunch mcl = mcl_bunch();
  mcl.insert("ole", 5.120, 1);
  assert(mcl.compare_string("ole:5.120", true));
  mcl.set_line_end(true, false); 
  assert(mcl.compare_string("ole:5.120\t$\n", true));//  mcl.print();     
  mcl.set_end_of_block();
  assert(mcl.compare_string("ole:5.120\t$\n\0", true));//  mcl.print();   
  mcl.free_mem();

  mcl = mcl_bunch();
  mcl.insert(50, 5.0, 1); 
  assert(mcl.compare_string("50:5.000", true));//  mcl.print();  
  mcl.free_mem(); 
 
  mcl = mcl_bunch(); 
  mcl.insert(50, 5.0, 1);
  assert(mcl.compare_string("50:5.000", true));//  mcl.print();  
  mcl.free_mem();

  mcl = mcl_bunch();
  mcl.insert("ole", "ole", 50, 0);
  assert(!mcl.has_data());
  mcl.free_mem();

  mcl = mcl_bunch();
  mcl.insert("ole", "ole", 50, 1);
  assert(mcl.compare_string("ole\tole\t50.000", true));
  mcl.free_mem();

  mcl = mcl_bunch();
  mcl.insert(12, 12, 50, 1);
  //  mcl.print();   
  assert(mcl.compare_string("12\t12\t50.000", true));
  mcl.free_mem();
#endif  
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
  //  if(false)   print_constants();
}   

