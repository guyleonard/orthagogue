#ifndef bucket_norm_h
#define bucket_norm_h
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
#include "list_norm.h"
#include "enum_write_list_t.h"
/**
   @file
   @struct bucket_norm
   @brief Data container holding changes in the normative array
   @author Ole Kristian Ekseth (oekseth)
   @date 02.11.2011 by oekseth (initial)
**/
struct bucket_norm {
  //! Holds the basis for the normalization values
  list_norm_t *arrNorm;  
  //! Deallocates the memory
  void free_mem(const uint taxon_length);
  //! Initializes the class
  static bucket_norm *init(list_norm_t *_arrNorm);

  //! Constructs the class given the input params.
  bucket_norm(list_norm_t *_arrNorm);

  //! Asserts the class given with test functions.
  static void assert_class(const bool print_info);
};
/**
   @brief Data container holding changes in the normative array
   @author Ole Kristian Ekseth (oekseth)
   @date 02.11.2011 by oekseth (initial)
**/
typedef struct bucket_norm bucket_norm_t;
#endif
