#include "mpi_read.h"
#ifdef USE_MPI

//! @return the file chunk:
char *mpi_read::get_file_chunk(loint &string_length ) {
  const loint chars_remaning = reading_file_end_position - reading_file_start_position;
  if(chars_remaning<=0) {string_length = 0; return NULL;}
  assert(block_length);
  string_length = block_length;
  if(chars_remaning < block_length) string_length = chars_remaning;
  reading_file_start_position += string_length; // Updates.
  char *string = new char[string_length+1]; string[string_length] = '\0'; // The end.
  MPI_File_read(fh, string,
		//number of objects to read=
		string_length, MPI_CHAR, MPI_STATUS_IGNORE);
  return string;
}
//! @return true if there is more data to read.
bool mpi_read::has_more_data_in_file_to_read() {
  const loint chars_remaning = reading_file_end_position - reading_file_start_position;
  if(chars_remaning<=0) {
    return false;
  }
  else return true;
}

/**
   @brief Reads the data
   @param <str>    The char buffer to put things into.
   @param <str_len>  The maximum number of elements to read.
   @return the number of elements read.
**/
loint mpi_read::read_chars(char *string, size_t string_length) {
  assert(string);
  assert(string_length);
  assert(block_length);
  const loint chars_remaning = reading_file_end_position - reading_file_start_position;
  if(chars_remaning<=0) {return 0;}
  loint chars_read = (loint)string_length;
  if(chars_remaning < (loint)string_length) chars_read = chars_remaning;
  //  memset(string+chars_remaning, '\0', string_length - chars_remaning);
  reading_file_start_position += chars_read; // Updates.
  MPI_File_read(fh, string,
		//number of objects to read=
		chars_read, MPI_CHAR, MPI_STATUS_IGNORE);
  if(reading_file_start_position >= reading_file_end_position) {
    const loint not_set = reading_file_start_position - reading_file_end_position;
    if(not_set) {
      memset(string+string_length-not_set, '\0', not_set);
    }
  }
  return chars_read;
}

  //! Deallocates the memory.
 void mpi_read::free_mem() {
    if(file_is_set) {
      MPI_File_close(&fh); 
      MPI_Type_free(&filetype);
      file_is_set = false;
    }
  }
//! The constructor:
mpi_read::mpi_read(MPI_Comm comm, int _myrank, char *file_name, loint _reading_file_start_position, loint _reading_file_end_position, loint _reading_file_length, loint _read_block_length) :
  myrank(_myrank), file_is_set(false), reading_file_start_position(_reading_file_start_position), reading_file_end_position(_reading_file_end_position), reading_file_length(_reading_file_length), block_length(_read_block_length)
{
  if(file_name) {// && reading_file_length && reading_file_end_position) {
    file_is_set = true; // Therefore the memory will later on be initialised.
    MPI_File_open(comm, file_name, MPI_MODE_RDONLY,  MPI_INFO_NULL, &fh);     
    //! Builds a file type retrieving the file reading from.
    MPI_Type_contiguous(block_length, MPI_CHAR, &filetype);
    MPI_Type_commit(&filetype); 
    MPI_Offset myrank_offset = reading_file_start_position;
    // TODO: Remove "native" and repalce with safer type!
    MPI_File_set_view(fh, myrank_offset, MPI_CHAR, filetype, "native",  MPI_INFO_NULL); 
  } else {
    log_builder::throw_warning(software_dependencies,__LINE__, __FILE__, __FUNCTION__,"Input not provided enabling reading of the file");
    assert(false);
  }
}

#endif
