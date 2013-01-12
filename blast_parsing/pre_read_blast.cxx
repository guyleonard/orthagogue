#include "pre_read_blast.h"
// -----------------------------------------------------------
// The interesting part of this class (private):
//! @return the file pointer.
FILE *pre_read_blast::get_file() {
  assert(blast_settings.FILE_NAME); // If name not set object used wrongly
  FILE *file = fopen(blast_settings.FILE_NAME, "rb");
  if(file) {
    return file;
  } else {
    printf( "!!\tUnable to open file named '%s' at line %d in class blast_parsing. Aborts.\n", blast_settings.FILE_NAME, __LINE__);
    fflush(stderr);
    exit(2);
  }
  return NULL;
}
 
//! @return true if exact positions will need to be calculated.
bool pre_read_blast::calculate_exact_positions(const loint file_size, const loint approx_size_each) {
  // (b) Reserves memory:
  file_start_pos = new loint[MACHINE_CNT];
  file_end_pos   = new loint[MACHINE_CNT];

  file_start_pos[0] = file_end_pos[0] = 0;
  assert(disk_buffer_size);
  const loint d_size = disk_buffer_size; 
  if(approx_size_each < d_size) {
    int myrank = 0;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
#endif
    if(myrank==0) { // Only write data for the root node, in order to avoid cluttering the termianl user otuput.
      char string[100]; memset(string, '\0', 100);
      sprintf(string, "Node designated id(%d) is not set doing the prosessing, probably due to the disk buffer size set to %u", myrank, disk_buffer_size);
      uint preferable_number_of_nodes = ((uint)file_size/(uint)disk_buffer_size);
      if(preferable_number_of_nodes) preferable_number_of_nodes-=1;
      fprintf(stderr, "!!\tAborts: If the disk buffer size(%d) was set manually, please read the manual, else reduce the number of nodes to %u nodes (file_size=%d/disk_buffer_size(%d).\n", disk_buffer_size, preferable_number_of_nodes, (int)file_size, (int)disk_buffer_size);
      log_builder::throw_warning(input_terminal, __LINE__, __FILE__, __FUNCTION__,string);
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif
    exit(1);  // Exits:
//     // The first machine_id reads all the data:
//     file_start_pos[0] = 0; file_end_pos[0] = file_size;
//     // Sets the rest of the macine id's to empty size avoiding them to read data:
//     for(int i = 1; i < (int)MACHINE_CNT; i++) {
//       file_start_pos[i] = 0; file_end_pos[i] = 0;
//     }
    return false;
  } else return true;
}

//! @return the file size.
loint pre_read_blast::get_file_size(char *file) {
  return log_builder::get_file_size(file, __LINE__, __FILE__);
}


// -----------------------------------------------------------
// Standard proecedures for getting the results:
//! Returns the starting position in the file for the machine id given.
loint pre_read_blast::get_start_pos(uint i) {
  // Tests some general properties that would be set if object is used correctly:
  assert(i < MACHINE_CNT);
  assert(file_start_pos);
  return file_start_pos[i];
}

//! Returns the end position in the file for the machine id given.
loint pre_read_blast::get_end_pos(uint i) {
  // Tests some general properties that would be set if object is used correctly:
  assert(i < MACHINE_CNT);
  assert(file_end_pos);
  return file_end_pos[i];
}

//! Returns number fo chars to read in the file for the machine id given.
loint pre_read_blast::get_read_length(uint i) {
  // Tests some general properties that would be set if object is used correctly:
  assert(i < MACHINE_CNT);
  assert(file_start_pos);
  assert(file_end_pos); 
  const loint length = file_end_pos[i] - file_start_pos[i];
  return length;
}

//! Prints the elmement given the index as param
void pre_read_blast::print_element(uint i) {
  if((i < MACHINE_CNT) && (file_start_pos) && (file_end_pos)) {
    printf("-\tThe file for macine_id[%u] is in [%llu, %llu]\n",
	   i, file_start_pos[i], file_end_pos[i]);
  }
}

// -----------------------------------------------------------
// The interesting part of this class (public):

/**
   @brief Produces a safe division of the blast file (into blocks).
   @remarks Procedure is:
   #- Tests some general properties set by the constructor, e.g. file exsistence.
   #- Reserves memory and checks if exact positions must be caluclated.
   #- Sets the approximate size for each machine id.
   #- Blocks are tests, ensuring there's a shift from one type of left-most protein to another. This avoids the need for the merging procedure in the mpi-implemntation having to consider overlapping regions.
   @author Ole Kristian Ekseth (oekseth)
   @date 02.01.2012 by oekseth (initial)
   @date 03.01.2012 by oekseth (asserts)
**/
void pre_read_blast::start_preprocessing(int myrank, char seperator) {    
  // Tests some general properties that would be set if object is used correctly:
  assert(MACHINE_CNT); // Number of machines > 0 for any purpose of this operator
  const loint file_size = get_file_size(blast_settings.FILE_NAME);
  if(MACHINE_CNT == 1) {
    if(!file_start_pos) file_start_pos = new loint[MACHINE_CNT];
    if(!file_end_pos)   file_end_pos   = new loint[MACHINE_CNT];
    file_start_pos[0] = 0;
    file_end_pos[0] =   file_size;
  } else {
    assert(file_size);           // Tests the files exsistence.
    loint approx_size_each = file_size/(MACHINE_CNT+1);
    // Reserves memory and checks if exact positions must be caluclated.    
    if(calculate_exact_positions(file_size, approx_size_each)) {
      approx_size_each = file_size/MACHINE_CNT; // Updates to improve error message:
      if(false) {
	// Comment: Code below included in order to examplify the goal of this code:
	printf("--\tThe indexes if the file had been split equal:\n");
	for(uint i = 0; i < MACHINE_CNT; i++) {
	  if(i) file_start_pos[i] = file_end_pos[i-1];
	  file_end_pos[i] = file_start_pos[i] + approx_size_each;
	  printf("#\tindex[%u] has start=%llu and end=%llu\n", i, file_start_pos[i], file_end_pos[i]);
	}
	printf("--\n");
      }
      // (d) Calculates the exact position in the file

      // TODO: Validate the routines for possible recuction of disk buffer size below!
      uint temp_disk_buffer_size = disk_buffer_size;
      if(disk_buffer_size > (1024*1024*20)) temp_disk_buffer_size = 1024*1024*20; // adjusts in order to remove a lag

      read_file_t *fileRead = read_file::init(temp_disk_buffer_size, blast_settings.FILE_NAME);
      log_builder::test_memory_condition_and_if_not_abort(fileRead!=NULL, __LINE__, __FILE__, __FUNCTION__);
      for(uint i = 0; i < MACHINE_CNT; i++) {
	taxon_list_t *listProteins = taxon_list::init();      
	log_builder::test_memory_condition_and_if_not_abort(listProteins!=NULL, __LINE__, __FILE__, __FUNCTION__);
	if(i!=0) { // Makes it consistent
	  file_start_pos[i] = file_end_pos[i-1];
	}
	// Moves the file to the assumed new break-point
	loint possible_end_location = file_start_pos[i] + approx_size_each;
	if(possible_end_location < file_size) {
	  if(fileRead->set_file_pointer_position(possible_end_location)) { // File was big enough for this operation
	    const loint jump_size = possible_end_location;
	    //	    const int ftell_position_in_file = fileRead->ftell_get_real_position_in_file();
	    assert(fileRead->ftell_get_real_position_in_file() == (int)jump_size); // They are eqaul

	    // Sets the end of the string
	    loint file_position = 0;
	    fileRead->set_overflow_buffer_to_emtpy();
	    do {
	      uint block_cnt = 0;
	      string_section_t *string = fileRead->read_file_blocks(seperator, block_cnt, false); 
	      // Tests if block has a shift from one type of left-most-protein to another
	      const loint chars_in_block_read =  string->more_than_one_leftmost_protein(listProteins, blast_settings);
	      file_position = jump_size + chars_in_block_read;
	      
	      // De-allocates the strin_section reserved:
	      string_section::close(string);
	    } while (file_position == 0);
	    file_end_pos[i] = file_position; // Updates the position.

#ifndef NDEBUG
	    if(!myrank) {
	      if(file_end_pos[i] > 10) {
		const int pos_of_newline = file_end_pos[i] -1;
		//! Asserts that a newline is the end:
		FILE *file_input = fopen(blast_settings.FILE_NAME, "r"); // assumes that reading in binary mode goes faster
		log_builder::test_file_reading_and_if_not_abort(file_input, blast_settings.FILE_NAME, __LINE__, __FILE__, __FUNCTION__);
		fseek(file_input, pos_of_newline, SEEK_SET); // positions the file
		char poss_newline = fgetc(file_input);
		if(poss_newline != EOF) {
		  if(poss_newline != '\n') {
		    log_builder::throw_warning(software_error, __LINE__, __FILE__, __FUNCTION__,"division of file into blocks erronous");
		    printf("!!\t%d\tDid not find thw newline ('\\n'=%d) at the expected place, instead we got:\t", (int)i, (int)'\n');
		    while(poss_newline != '\n' && poss_newline != EOF) {
		      poss_newline = fgetc(file_input);
		      printf("%c", poss_newline);
		      //		    printf("%c (%d),", poss_newline,  (int)poss_newline);
		    }
		    printf("\n");
		    assert(false);
		  }
		} 
		fclose(file_input);
	      }
	    }
#endif
	    // De-allocates the memory reserved:	  
	    if(listProteins) taxon_list::close(listProteins);
	  } else file_end_pos[i] = file_size;
	} else file_end_pos[i] = file_size;
	if(listProteins) taxon_list::close(listProteins);
	if(file_end_pos[i] == file_start_pos[i]) {
	  if(myrank==0) { // Only write data for the root node, in order to avoid cluttering the termianl user otuput.
	    char string[100]; memset(string, '\0', 100);
	    sprintf(string, "Node designated id(%d) is not set doing the prosessing, probably due to the disk buffer size set to %u", myrank, disk_buffer_size);
	    uint preferable_number_of_nodes = i;
	    if(preferable_number_of_nodes) preferable_number_of_nodes -=1;
	    for(uint o = 0; o < i+1; o++) {
	      print_element(o); // Prints the allocations.
	    }
	    fprintf(stderr, "!!\tAborts: If the disk buffer size(%d) was set manually, please read the manual, else reduce the number of nodes to %u nodes (file_size=%d/disk_buffer_size(%d).\n", disk_buffer_size, preferable_number_of_nodes, (int)file_size, (int)disk_buffer_size);
	    log_builder::throw_warning(input_terminal, __LINE__, __FILE__, __FUNCTION__,string);
	  }
#ifdef USE_MPI
	  MPI_Finalize();
#endif
	  exit(1);
	}
	//      printf("in pre_read_blast at line %d\n", __LINE__);       print_element(i);
      }

      // De-allocates the memory reserved:
      if(fileRead)     read_file::close(fileRead);
    }
  }

}


// -----------------------------------------------------------
// Initiators and de-allocations procedures:

//! Initiates an object of type pre_read_blast.
pre_read_blast *pre_read_blast::init(uint _disk_buffer_size, uint MACHINE_CNT, tsettings_input_t blast_settings) {
  return new pre_read_blast(_disk_buffer_size, MACHINE_CNT, blast_settings);
}  

//! De-allocates the memory reserved the object of type pre_read_blast.
void pre_read_blast::close(pre_read_blast *&obj) {
  if(obj){obj->free_memory(), delete obj, obj = NULL;}
}

//! De-allocates the memory reserved for this object.
void pre_read_blast::free_memory() {
  if(file_start_pos) {delete [] file_start_pos; file_start_pos = NULL;}
  if(file_end_pos) {delete [] file_end_pos; file_end_pos = NULL;}
}

/**
   @brief The constructor:
   @param <_disk_buffer_size> The number of chars to approximate read for each 'run'
   @param <_MACHINE_CNT> The number of machines (nodes) to use for the MPI operation.
   @param <f_name> The name of the file to evaluate.
   @author Ole Kristian Ekseth (oekseth)
*/
pre_read_blast::pre_read_blast(uint _disk_buffer_size, uint _MACHINE_CNT, tsettings_input_t _blast_settings)
  : disk_buffer_size(_disk_buffer_size), MACHINE_CNT(_MACHINE_CNT), blast_settings(_blast_settings), file_start_pos(NULL), file_end_pos(NULL) {}

// -----------------------------------------------------------
// Methods for testing:

void pre_read_blast::assert_private_parts() {
#ifdef assert_code

#endif
}

/**
   @brief The main test function for this class.
   @remarks This method includes:
   - Formalized tests for verification that a block starts at a line.
   - Valgrind used verifying memory usage.
   - Examples of how this class may be used.
   @author Ole Kristian Ekseth (oekseth)
   @date 03.01.2012 by oekseth.
**/
void pre_read_blast::assert_class(const bool print_info) {
  const static char *class_name = "pre_read_blast";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  // (1) Initiates an object: (a) opening the file, (b) building a list to store data in
  const int INDEX_IN_FILE_FOR_PROTEIN_NAME = 0;
  const int INDEX_IN_FILE_FOR_TAXON_NAME = 1;
  const bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE = false;
  const char SEPERATOR = '_';
  char *FILE_INPUT_NAME = "all.blast";
  const uint DEFAULT_NUMBER_OF_COLUMNS_IN_NAME = 6;
  tsettings_input_t blast_settings = tsettings_input(INDEX_IN_FILE_FOR_PROTEIN_NAME, INDEX_IN_FILE_FOR_TAXON_NAME, USE_LAST_BLAST_CLOMUN_AS_DISTANCE, SEPERATOR, FILE_INPUT_NAME,DEFAULT_NUMBER_OF_COLUMNS_IN_NAME);
  const uint MACHINE_CNT = 2; 
  const int disk_buffer_size = 1024;
  pre_read_blast *obj = pre_read_blast::init(disk_buffer_size, MACHINE_CNT, blast_settings);
  obj->assert_private_parts();
  // (2) Starts the divsion into blocks:
  obj->start_preprocessing(4, SEPERATOR);
  assert(disk_buffer_size == 1024); // Set as a mandatory criterium in order to get the below example working!
  // (3) Get the block positions for each machine (node). This exemplifies the usage of the data gained above
  if(false) {
    // Comment: This code included in order to show example of usage.
    for(uint i = 0; i < MACHINE_CNT; i++) {
      printf("\tMachine[%u] starts at %llu and ends at %llu, giving a length of %llu chars to read.\n", i,  
	     obj->get_start_pos(i), obj->get_end_pos(i), obj->get_read_length(i));    
    }
  }
  FILE *file = fopen(FILE_INPUT_NAME, "rb");
  const loint f_size = pre_read_blast::get_file_size(FILE_INPUT_NAME);
  char *data = new char[f_size];
  fread(data, f_size, 1, file);
  for(uint i = 0; i < MACHINE_CNT; i++) {
    const loint p_start = obj->get_start_pos(i);
    const loint p_end = obj->get_end_pos(i);
    //    obj->print_element(i);
    if(false) {
      printf("starts:\n----\n");
      for(loint k = p_start; k < p_end; k++) {
	putchar(data[k]);
      }
    }
    if(p_start > 0) assert(data[p_start-1] == '\n');
    if(p_end > 0 &&
       (p_end!= (f_size-1))
       ) assert(data[p_end-1] == '\n');
    if(i != 0) assert(obj->get_end_pos(i-1) == obj->get_start_pos(i));
  }
  delete [] data; data = NULL;
  fclose(file);
  // (4) Closes the object, releasing the memory:
  pre_read_blast::close(obj);
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

