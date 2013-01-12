#ifndef bucket_pipe_binary_h
#define bucket_pipe_binary_h
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
#include "types.h"
#include "list_file_parse.h"
#include "list_norm.h"
/**
   @file
   @struct bucket_pipe_binary
   @brief A transport container
   @ingroup filter_transport_container
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by oekseth (init).
   @date 25.12.2011 by oekseth (cleanup).
**/
struct bucket_pipe_binary {
  //! @return the **norm object holding the basis for the normalization procedure.
  list_norm_t *arrNorm;
  //! Holds the filtered blast data.
  list_file_struct_t *structData;
  //! @return true if data is set
  bool isNotEmpty();
  //! De-allocates memory for the list_file_struct object.
  void free_structdata();  
  //  //! Deallocates memory for the **norm object given the length of it.
  //  void free_arrNorm(const uint taxon_length);
  //! Deallocates the memory for this object.
  void free_mem(const uint taxon_length);

  //! The constructor.
  bucket_pipe_binary();
  //! The constructor.
  bucket_pipe_binary(list_file_struct_t *&_structData);
  //! The constructor.
  bucket_pipe_binary(list_norm_t *_arrNorm,  list_file_struct_t *&_structData);
};

/**
   @brief A transport container
   @ingroup filter_transport_container
   @author Ole Kristian Ekseth (oekseth)
   @date 25.12.2011 by oekseth (cleanup).
**/
typedef struct bucket_pipe_binary bucket_pipe_binary_t;

#endif
