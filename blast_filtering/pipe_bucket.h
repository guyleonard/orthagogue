#ifndef pipe_bucket_h
#define pipe_bucket_h
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
#include "taxon_pair.h"
/**
   @file
   @class pipe_bucket
   @brief Produces buckets of numbered items to parse.
   @ingroup pipe_filters
   @details Goal is to make the threads in the rest of the pipe effective:
   @todo Should be considered removed, replaced by preestimated blocks used by each of the consquative threads
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/
class pipe_bucket: public tbb::filter {
 private:
  list_file_parse_t *listParseData; 
  list_file_struct_t *listStructData; 
  taxa_t *listTaxa;
  log_builder_t *log;
  const bool is_inpa;
  /**
     @brief Produces buckets of numbered items to parse. Goal is to make the threads in the rest of the pipe effective:
     @todo Should be considered removed, replaced by preestimated blocks used by each
     of the consquative threads.
     @date 21.12.2010 by oekseth
  */
  struct taxon_pair* init_taxon_p(const bool inpa_ops, uint taxon_start, uint taxon_end,
				  uint taxon_length, uint n_threads);  
  //! Holds the blocks to be used for the parsing operation.
  struct taxon_pair *listPair;
  uint list_pair_pos;
 public:
  //! The constructor.
  pipe_bucket(const bool inpa_ops, uint taxon_start, uint taxon_end,
	      uint taxon_length, log_builder_t *_log,
	      list_file_parse_t *&_listParseData,
	      list_file_struct_t *&_listStructData, taxa_t *_listTaxa, uint n_threads
	      );

  //! The method of parallisation.
  void* operator()(void* item);
  //! De-allocates the memory of this object.
  void free_data() {free(listPair);}
  //! De-allocates the memory of this object.
  void free_mem() {free_data();}
#ifdef assert_code
  //! The test function for the private parts of this class.
  void assert_private_parts();
#endif
  //! The main test function for this class  
  static void assert_class(const bool print_info);
};

#endif
