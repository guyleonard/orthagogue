#include "bucket_pipe_binary.h"
//! @return true if data is set
bool bucket_pipe_binary::isNotEmpty() {return structData != NULL;}

//! De-allocates memory for the list_file_struct object.
void bucket_pipe_binary::free_structdata() {list_file_parse<rel>::close(structData, false);}

//! Deallocates the memory for this object.
void bucket_pipe_binary::free_mem(const uint taxon_length) {
  list_file_parse<rel>::close(structData, false);
  list_norm::close(arrNorm);
}

//! The constructor.
bucket_pipe_binary::bucket_pipe_binary() :
  arrNorm(NULL), structData(NULL)
{}

//! The constructor.
bucket_pipe_binary::bucket_pipe_binary(list_file_struct_t *&_structData) :
  arrNorm(NULL), structData(_structData)
{
  _structData = NULL;
}

//! The constructor.
bucket_pipe_binary::bucket_pipe_binary(list_norm_t *_arrNorm,  list_file_struct_t *&_structData) :
  arrNorm(_arrNorm), structData(_structData)
{
  _structData = NULL;
}
