#ifndef blast_extractors_h
#define blast_extractors_h
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
#include "parse.h"
#include "tsettings_input.h"
#include "../configure.h"
/**
   @file
   @brief Functions for parsing the block of chars given as input from a previous run.
   @ingroup parsing_ops
   @remarks Uses object of type Parse storing th results in
   @author Ole Kristian Ekseth (oekseth)
**/
//! \{    @ingroup parsing_ops

/**
   @class blast_extractors
   @brief Functions for parsing the block of chars given as input from a previous run.
   @ingroup parsing_ops
   @remarks Uses object of type Parse storing th results in
   @author Ole Kristian Ekseth (oekseth)
**/
class blast_extractors {
 public:

/**
   @brief Prints the segment given.
   @param <start_p> The Position in the memory string to start the reading.
   @param <end> The Position in the memory string to end the reading.
   @author Ole Kristian Ekseth (oekseth)
 **/
  static void print_segment(char *start_p, char *end);

/**
   @brief Prints each char seperately.
   @param <start_p> The Position in the memory string to start the reading.
   @param <end> The Position in the memory string to end the reading.
   @remarks Sometimes special chars affects the data. This methods primary usage is for visual validation that the input corresponds the the programmers expecation.
   @author Ole Kristian Ekseth (oekseth)
 **/
  static void print_segment_stepwise(char *start_p, char *end);

  /**
     @brief Gets the labels for a given protein, using a line in the blast file.
     @param <first_column> True if it's the leftmost protein in the given blast-file line
     @param <p> The Parse object to set the names into.
     @param <pos_column_start> The memory position to start the work from.
     @param <logical_end> The last address where data may be find: used to avoid errors.
     @param <b> The tsettings_input_t object iot parse correctly.
     yy   @return The position after the data were the given column were located.
  **/
  static char *getIDColumn(const bool first_column, Parse &p, char *pos_column_start, char *logical_end, tsettings_input_t b);

  /**
     @brief Gets the labels for a protein pair, using a line in the blast file.
     @param <p> The Parse object to set the names into.
     @param <pos_column_start> The memory position to start the work from.
     @param <logical_end> The last address where data may be find: used to avoid errors.
     @param <blast_settings> The tsettings_input_t object iot parse correctly.
     @return True if data were read correctly.
  **/
  static bool getParseHeaders(Parse &p, char *&pos_column_start, char *logical_end, tsettings_input_t blast_settings);

  /**
     @brief Methods for assertion of this class.
     @remarks Included for further use.
  **/
  static void assert_class(const bool print_info) {
    const static char *class_name = "blast_extractors";
    if(print_info) printf("--\tStarts testing class %s\n", class_name);
    // TODO: include some code here.
    if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
  }
};

/**
   @typedef blast_extractors_t
   @brief Functions for parsing the block of chars given as input from a previous run.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
**/
typedef class blast_extractors blast_extractors_t;
//! \}
#endif
