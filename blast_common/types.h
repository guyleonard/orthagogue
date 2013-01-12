#ifndef types_h
#define types_h
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
/**
  @file
  @brief Defines commonly used types in the project.
  @ingroup common
  @author Ole Kristian Ekseth (oekseth)
  @date 21.12.2010 by Ole Kristian Ekseth (init)
  @date 28.12.2011 by oekseth (clean-up)
 */
/**  @ingroup common **/
/** \{ **/
#include "libs.h"
#include "o_rel.h"

//! Depricated, but some methods still uses it holding the memmory address.
typedef long int mem_loc; 

//typedef unsigned short int overlap_t; //! Defines the maximum size of the overlap.

//! For convenience, giving the system easier job casting error if value is negative.
typedef unsigned int uint; 

/**
  @brief Defines the maximum size of the overlap.
  @todo Changed to 'short int' as valgrind complains on 'usigned short int'.
  @remarks If changed from other type, MPI structures using this must be updated appropriatly, as the complainer does not have means detecting such inconsistencies!.
*/
typedef uint overlap_t; 

//! Holding the internal address of a data in a file, requires a a large number. 
typedef unsigned long long int loint; 
//typedef long long unsigned int loint; 
//! Holding the internal address of a data in a file, requires a a large number. 
typedef long int lint; 
/** \} **/
#endif
