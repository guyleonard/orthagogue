#ifndef cmd_list_h  
#define cmd_list_h
/**
   @file
   @class cmd_list
   @brief Defines the interface between arguments and variables. Used by the programmer directly.
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
 * along with orthAgogue. If not, see <http://www.gnu.org/licenses/>.
 */
#include "../configure.h"
#include "../blast_common/types.h"
#include "../blast_common/libs.h"
#include "../blast_common/constants.h" 
//#include "types.h"
//#include "constants.h" 
#include "cmd_constants.h"
#include "cmd_argument.h"
#include <time.h> 
#include <stdio.h>
/**
   @class cmd_list
   @brief Interface between arguments and variables. Used by the programmer directly.
   @remarks This is an easy implementation, due to the observation that
   reading the command line input is not a critical part of the software running time.
   The following class should therefore not be used in parts requireing a high performance
   with regard to string parsing. In such a situation, a suffix tree would be advisable.
   @ingroup terminal
   @author Ole Kristian Ekseth (oekseth)
   @date 27.12.2010 by oekseth (cleanup)
**/
class cmd_list {
 private:
  class cmd_argument *list; // the list holding each set argument
  uint list_pos; // the next postion to assgn to the list
  uint list_size; // the size of the lsit
  // Variables updated during parsing:
  bool print_argument_after_init; // if set to true, prints the argument list
  bool print_man_page; // if set to true, prints the man page
  bool print_help; // if set to true, prints the help list
  char *DEFAULT_OPTION_NAME;
  uint DEFAULT_OPTION_NAME_COUNT;
  char *FILE_INPUT_NAME;
  /**
     @brief Iterates through the argument list, setting the variable if found  
     @return true if argument was found
  **/
  bool insert_argument(char *&short_id, char *&long_id, char *string, bool &remove_previous_references);
 public:
  //! Prints the arguments to several 'outs':
  void print_arguments(FILE *output_help, FILE *man_page_file);
  /**! Binds the input arguments with variables: */
  void set_arguments(char **array, int array_size, bool create_man_page);
  //! Returns the size of the argument list
  uint getListSize(){return list_pos;}  
  //! Set the argument list to NULLB
  void clean_cmd_argument() {delete [] list;list_pos=0, list_size =0;} 
  //! Called when arguments are to be added
  void add_cmd_argument(class cmd_argument arg); 
  /**! Builds a man page header for the software  */
  FILE *create_man_page_header();
  /**! Builds the man page tail, and ends the operation of the man page building  */
  void create_man_page_tail(FILE *&file);
  //! De-allocates the memory reserved
  void free_mem(){} // Included for future use
  //! De-allocates memory for the object given as param.
  static void close(cmd_list *&obj) {if(obj){obj->free_mem(), delete obj, obj=NULL;}}
  //! Constructor
 cmd_list(char *_DEFAULT_OPTION_NAME, uint &_DEFAULT_OPTION_NAME_COUNT, char *_FILE_INPUT_NAME) :
  list(NULL)
    , list_pos(0), list_size(0), DEFAULT_OPTION_NAME(_DEFAULT_OPTION_NAME), DEFAULT_OPTION_NAME_COUNT(_DEFAULT_OPTION_NAME_COUNT),
    FILE_INPUT_NAME(_FILE_INPUT_NAME)
    {
      list = new cmd_argument[5];
      list_pos = 0, list_size = 5;
      print_argument_after_init = false; // if set to true, prints the argument list_list::print_argument_after_init = false; // if set to true, prints the argument list
      print_man_page = false; // if set to true, prints the man page
      print_help = false; // if set to true, prints the help list
      if(false)   print_cmd_constants();
    }
  //! The main test function for this class  
  static void assert_class(const bool show_positive_results);
};

#endif
