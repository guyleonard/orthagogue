/**
   @file
   @brief Merges containers building a filtered set of orthologs- and inparalogs.
   @ingroup pipe_filters
   @author Ole Kristian Ekseth (oekseth)
**/
#ifndef pipe_merge_h
#define pipe_merge_h
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
#include "taxa.h"
#include "list_file_parse.h"
#include "bucket_pipe_binary.h"
#include "pipe_t.h"
/**
   @class pipe_merge
   @brief Merges containers building a filtered set of orthologs- and inparalogs.
   @ingroup pipe_filters
   @remarks Inputs an object of type 'bucket_pipe_binary', containing blocks of information to be merged.
   - Merges containers of type list_file_struct and type **norm.
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
   @date 25.12.2011 by oekseth (cleanup).
**/
class pipe_merge: public tbb::filter {
 private:
  taxa_t *listTaxa;
  log_builder_t *log;
  const uint taxon_length;
  const bool USE_EVERYREL_AS_ARRNORM_BASIS;
  const pipe_t PIPE_TYPE;
  list_norm_t *arrNorm;
  list_file_struct_t *fileStruct;
  //  list_file_struct_t *listStructData;
  loint debug_size_fileStruct_at_init;
  //! Allocates- and intiates memory for the 'norm_t'structure
  list_norm_t *init_arrNorm();
  /**! Updates the global list of the normalization values  */
  void mergeArrNorm(list_norm_t *norm);
 public:
  //! Sets the *list_norm object with the param given.
  void set_arrNorm(list_norm_t *norm) {
    arrNorm = norm;
  }
  //! @return the *list_norm object
 list_norm_t *get_arrNorm(){return arrNorm;}
  //! Sets the *list_file_struct object with the param given.
  void setFileStruct(list_file_struct_t* arg) {fileStruct = arg;}
  /**
     @brief When the list_file_struct object is returned, it's the end-of-life for pipe_merge, and a log file is generated.
     @return the processed *list_file_struct object.
  **/
  list_file_struct_t *getFileStruct(int n_threads, bool MODE_PAIRWISE_OUTPUT_ABC, bool MODE_PAIRWISE_OUTPUT_MCL, char *FILE_BINARY_LOCATION);

  //! The method of parallisation.
  /*overrride*/ void* operator()(void* item);
  //! The constructor.
  pipe_merge(uint _taxon_length, const bool _use_everyrel_as_arrnorm_basis,
	     const pipe_t type, log_builder_t *_log,
	     list_file_struct_t *&_listStructData, taxa_t *_listTaxa);
#ifdef assert_code
  //! The test function for the private parts of this class.
  void assert_private_parts();
#endif
  //! The main test function for this class  
  static void assert_class(const bool print_info);
};
#endif
