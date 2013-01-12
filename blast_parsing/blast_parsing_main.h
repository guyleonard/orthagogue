#ifndef parsing_main_h
#define parsing_main_h
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
#include "blast_parsing.h"
#include "cmd_list.h"
#include "bp_container.h"
/**
   @file
   @brief The launcher of the parsing library.
   @ingroup parsing pipe_parsing parsing_container transfer_container parsing_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 30.12.2011 by oekseth (initial)
**/


// ---------------------(methods)------------------------------
/**
   @brief This file is regarded as the launcher of the software, implying using the librariy  'blast_parsing'. 
   @remarks In order to understand the above library, this file might be of use the the developer.
   @author Ole Kristian Ekseth (oekseth)
   @date 27.12.2011 by oekseth (cleanup)
**/
int main(int argc, char *argv[]);
#endif
