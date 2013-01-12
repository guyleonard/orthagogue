#ifndef tsettings_input_h
#define tsettings_output_h
/**
   @file
   @brief Holds the properties for the output file.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
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
   @struct tsettings_output
   @brief Holds the properties for the output file.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
struct tsettings_output {
  // TODO: For future inmplementation in order to get return values.
  // ....
  //! The constructor
  tsettings_output();
};

/**
   @typedef tsettings_output_t
   @brief Holds the properties for the output file.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
typedef struct tsettings_output tsettings_output_t;
#endif
