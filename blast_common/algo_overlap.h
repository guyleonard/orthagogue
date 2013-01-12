#ifndef algo_overlap_h
#define algo_overlap_h
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
/**
   @file
   @brief Specific methods describing the overlap filtering procedure
   @author Ole Kristian Ekseth (oekseth)
   @date 02.11.2011 by oekseth (initial)
**/
//!    aboveOverlapLimit -- returns true if data is to be used.
bool aboveOverlapLimit(overlap_t overlap_left_left, overlap_t overlap_left_right,
		       overlap_t overlap_right_right, overlap_t overlap_right_left, short int AMINO_LIMIT, bool PRINT_OVERLAP_VALUES_ABOVE);

/**
   aboveOverlapLimit_improved -- returns true if data is to be used.
*/
bool aboveOverlapLimit_improved(overlap_t overlap_left_left, overlap_t overlap_left_right,
				overlap_t overlap_right_right, overlap_t overlap_right_left,   overlap_t recip_left_right_left, overlap_t recip_right_left_right, short int AMINO_LIMIT, bool PRINT_OVERLAP_VALUES_ABOVE);

  //! The main test function for this class  
void assert_overlap_module(const bool print_info);
//  static void assert_class(const bool print_info);
#endif
