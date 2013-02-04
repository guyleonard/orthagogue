#include "pipe_write.h"
pipe_write::pipe_write(log_builder_t *_log, taxa_t *listTaxa, int taxon_length, char *FILE_BINARY_LOCATION,
		       bool MODE_PAIRWISE_OUTPUT_ABC, bool MODE_PAIRWISE_OUTPUT_MCL,
		       bool PRINT_IN_ABC_FORMAT, bool PRINT_IN_MCL_FORMAT, mcl_t TYPE_OF_RESULTFILE_TO_STDOUT, bool _SORT_ABC_DATA) :
  tbb::filter(/*is_serial*/true), log(_log), SORT_ABC_DATA(_SORT_ABC_DATA)
#ifndef NDEBUG
  , file_chunk_received(NULL)
  , file_chunk_header(NULL)
#endif
{
#ifndef NDEBUG
  file_chunk_received = new list_file_chunk(MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT, SORT_ABC_DATA);
  file_chunk_header   = new list_file_chunk(MODE_PAIRWISE_OUTPUT_ABC, MODE_PAIRWISE_OUTPUT_MCL, PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT, SORT_ABC_DATA);
#endif
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes);
  assert(number_of_nodes);
  mpi_file_list = mcl_format::init_file_array(TYPE_OF_RESULTFILE_TO_STDOUT, FILE_BINARY_LOCATION, PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT); 
  int number_of_nodes = 1; MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes); 
  int myrank = 0;  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
  //! The 'root' node is given responsibility writing the header file
  if((number_of_nodes > 1 && myrank == 0) || (number_of_nodes==1)) {
    const uint actual_chars_inserted = mcl_format::write_file_headers(mpi_file_list, listTaxa, taxon_length, PRINT_IN_MCL_FORMAT); // Creates the headers, and writes them to the files
    const uint expected_chars_inserted = file_chunk_header->append_header_sizes(listTaxa, taxon_length);
    assert(actual_chars_inserted == expected_chars_inserted);
  }
#else
  out_file = mcl_format::init_file_array(TYPE_OF_RESULTFILE_TO_STDOUT, FILE_BINARY_LOCATION, PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT); 
  const uint actual_chars_inserted =  mcl_format::write_file_headers(out_file, listTaxa, taxon_length, PRINT_IN_MCL_FORMAT); // Creates the headers, and writes them to the files
#ifndef NDEBUG
    const uint expected_chars_inserted = file_chunk_header->append_header_sizes(listTaxa, taxon_length);
    assert(actual_chars_inserted == expected_chars_inserted);
#endif
#endif
  debug_dump_result = pipe_struct_result();

}

/**! Closes the files, adding the final trailings to the mcl files    */
void pipe_write::free_mem(bool SORT_ABC_DATA, char *FILE_BINARY_LOCATION, int n_threads) {
#ifdef USE_MPI 
  mcl_format::close_init_file_array(mpi_file_list, none, SORT_ABC_DATA, FILE_BINARY_LOCATION);
#else
  mcl_format::close_init_file_array(out_file, none, SORT_ABC_DATA, FILE_BINARY_LOCATION);
#endif
#ifndef NDEBUG
  FILE *f_log = log_builder::get_log_file_pointer(n_threads, "pipe_dump_write", __FILE__, __LINE__, true);
  if(f_log) {      
    fprintf(f_log, "--\t The dumping operation resulted in the strings having the following number of chars:\n");
    debug_dump_result.print_result(f_log);
    fclose(f_log);
  }
  bool verify_file_lengths_correctness = true;
#ifdef USE_MPI
  int number_of_nodes = 1; MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes); 
  int myrank = 0;  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
  //! Only the first node is to do the verification
  if(number_of_nodes > 1 && myrank != 0) {
    verify_file_lengths_correctness = false;
  }
#endif
  assert(file_chunk_received);
  assert(file_chunk_header);
  //! Compare the calculated file sizes with the actual file sizes:
  if(verify_file_lengths_correctness) {assert(file_chunk_received->compare_sizes_with_files(file_chunk_header, FILE_BINARY_LOCATION));}
  //! De-allocates the given data:
  if(file_chunk_received) {file_chunk_received->free_memory(); delete file_chunk_received; file_chunk_received=NULL;}
  if(file_chunk_header)   {file_chunk_header->free_memory();   delete file_chunk_header;   file_chunk_header=NULL;}
#endif
}

/**
   @brief The method of parallisation.
   @remarks Uses The "received object from pipe" is given the file pointer, using this to write data, before the received object being de-allocated.
**/
void* pipe_write::operator()(void* item) {
  static uint elements_inserted = 0;
  const clock_t tstart = times(NULL);
  mcl_format* mcl = static_cast<mcl_format*>(item);
  if(mcl != NULL) {
#ifndef NDEBUG
    //! Uses the actual string lengths for calculation
    mcl->log_update_result_cheme(debug_dump_result);
    debug_dump_result.append_meta_data(elements_inserted++, 0);
    //! Uses the calculated string length for calculation.
    assert(file_chunk_received);
    file_chunk_received->merge_lists(mcl->get_file_chunk());
    //! Compare the two lists, and verify they are the same:
    loint *two = NULL; debug_dump_result.get_list(two);
    assert(true == file_chunk_received->compare_with_list(two));
#endif

#ifdef USE_MPI 
    mcl->write_file(mpi_file_list); // writes the data to the file
#else
    mcl->write_file(out_file); // writes the data to the file
#endif
    mcl_format::free_class(mcl); // frees the memory reserved
  }
  const clock_t tend= times(NULL); log->append_measurement(write_resultfile, 2, tstart, tend);
  return NULL;
}

#ifdef assert_code
void pipe_write::assert_private_parts() {
  ;
}
#endif
void pipe_write::assert_class(const bool print_info) {
  const static char *class_name = "pipe_write";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  pipe_write temp();
  //  temp.free_mem();
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
