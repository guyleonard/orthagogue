/**
   @file
   @brief Builds (writes) the result files, consisting of the strings given.
   @ingroup pipe_filters
   @author Ole Kristian Ekseth (oekseth)
**/
#ifndef pipe_write_h
#define pipe_write_h
/*
 * Copyright 2012 Ole Kristian Ekseth (oekseth@gmail.com)
 *
 * This file is part of orthAgogue.
 *
 * orthAgogue is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * orthAgogue is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with orthAgogue
. If not, see <http://www.gnu.org/licenses/>.
 */
#include "log_builder.h"
#include "types.h"
#include "libs.h"
#include "taxa.h"
#include "file_parse.h"
#include "build_string.h"
#include "bucket_norm.h"
#include "mcl_format.h"
#include "enum_mcl.h"
#include "pipe_struct_result.h"
#include "list_file_chunk.h"
/**
   @class pipe_write
   @brief Builds (writes) the result files, consisting of the strings given.
   @ingroup pipe_filters
   @remarks Processes in serial the block of strings sent from class 'pipe_struct'.
   @author Ole Kristian Ekseth (oekseth)
   @date 30.09.2009 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 13.08.2012 by oekseth (support for MPI).
**/
class pipe_write: public tbb::filter {
 private:
  log_builder_t *log;
  FILE** out_file;
  pipe_struct_result debug_dump_result;
#ifdef USE_MPI
  MPI_File **mpi_file_list;
  int myrank;
  int number_of_nodes;
#endif
#ifndef NDEBUG
  /**
     @brief Hold the number of chars included for each of the files written. 
     @remarks Accessed via procedures found in class mcl_format
  **/
  list_file_chunk *file_chunk_received; 
  list_file_chunk *file_chunk_header; 
#endif
 public:
  /**! Closes the files, adding the final trailings to the mcl files. De-allocates memory reserved.    */
  void free_mem(bool SORT_ABC_DATA, char *FILE_BINARY_LOCATION, int n_threads);
  /**
     @brief The method of parallisation.
     @remarks Uses The "received object from pipe" is given the file pointer, using this to write data, before the received object being de-allocated.
  **/
  /*override*/void* operator()(void* item);
  //! The constructor:  
  pipe_write(log_builder_t *_log, taxa_t *listTaxa, int taxon_length, char *FILE_BINARY_LOCATION, 
	     bool MODE_PAIRWISE_OUTPUT_ABC, bool MODE_PAIRWISE_OUTPUT_MCL,
	     bool PRINT_IN_ABC_FORMAT, bool PRINT_IN_MCL_FORMAT, mcl_t TYPE_OF_RESULTFILE_TO_STDOUT);
  
#ifdef assert_code
  //! The test function for the private parts of this class.
  void assert_private_parts();
#endif
  //! The main test function for this class  
  static void assert_class(const bool print_info);

};

#endif

