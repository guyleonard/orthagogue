#ifndef pipe_norm_h
#define pipe_norm_h
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
#include "../configure.h"
#include "types.h"
#include "tbb_libs.h"
#include "log_builder.h"
#include "taxa.h"
#include "pipe_t.h"
#include "bucket_norm.h"
/**
   @file
   @class pipe_norm
   @brief Updates an norm_t** object, if it's not done so in a previous phase.
   @ingroup pipe_filters
   @todo If the  norm_t** object is already up to data, this step is a waste of time. Consider using an altnerative approach.
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)
**/
class pipe_norm: public tbb::filter {
 private:
  log_builder_t *log;
  const uint taxon_length;
  const bool USE_EVERYREL_AS_ARRNORM_BASIS;
  const pipe_t PIPE_TYPE;
  //! Updates the global list of the normalization values
  void mergeArrNorm(list_norm_t *norm);
 public:
  //! The object containing the basis for the normalization procedure.
  list_norm_t *arrNorm;  
  //! The method of parallisation.
  void* operator()(void* item);
  //! The constructor:
  pipe_norm(uint _taxon_length, const bool _use_everyrel_as_arrnorm_basis,
	    const pipe_t type, log_builder_t *_log);

  //! The test function for the private parts of this class.
  void assert_private_parts();
  //! The main test function for this class  
  static void assert_class(const bool print_info);
};
#endif
