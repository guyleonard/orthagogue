#ifndef cmd_argument_h
#define cmd_argument_h
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
#include "../blast_common/types.h"
#include "../blast_common/libs.h"
#include "../blast_common/constants.h" 
//#include "types.h"
//#include "constants.h" 
#include "cmd_constants.h"
#include "enum_argument_types.h"
/**
   @file
   @class cmd_argument
   @brief Makes internal variables accessable for the user. 
   @ingroup terminal
   @remarks Provides a data wrapper for each of the possible arguments on the command line
   @author Ole Kristian Ekseth (oekseth)
   @date 27.12.2010 by oekseth
**/
class cmd_argument {
 private:
  char *description; // the textual description of the variable
  char *identifier_short; // the primary identifier. E.g. "-f <file_name>"
  char *identifier_long; // the secondary identifier. E.g. "--file_name <file_name>"
  ARGUMENT_TYPES type; // holds informat9ion on how to validate the input argument
  int *argument_i; // Holds the reference to he argument given.
  bool *argument_b; // Holds the reference to he argument given.
  char **argument_cc; // Holds the reference to he argument given.
  char *argument_cc_s; // Holds the reference to he argument given.
  float *argument_f; // Holds the reference to he argument given.
  bool data_is_set; // set to true if data is set: speeds up the search
  char *sub_section_name; // the name of the subsection whom the argument is to be printed below: set to NULL when the argument is printed
  //! Returns the correct usage of the enums, correspondning to the rules of validations
  char *get_correct_usage();
  //! Prints error message if the parsing did not correspond to the specification, found in the function 'get_correct_usage()'
  void print_wrong_value(char *identifier, char *value);
  //! Inserts the value if it conforms to the input type set, else prints an error.
  void insert_value(char *identifier, char *value);
  //! Returns true if argument has lphanumeric vlaues, assuming that the uinput has a '\0' at its end
  bool is_alpha_num(char *value);
  /**! Returns true if argument has numeric values, assuming that the uinput has a '\0' at its end  */
  bool is_numeric(char *value, bool decimals);
  /** Inititates the arguments for the constructor:  */
  void init_arguments(ARGUMENT_TYPES type, void *_arg);
 public:
  //! Returns the short identifier
  char *getShortIdentifier() {return   identifier_short;} // the primary identifier. E.g. "-f <file_name>"
  /**
     @brief Inserts an argument if the text is found euqal to the input
     @param <_identifier_short> The identifier reqcognised by the prefix '-'
     @param <_identifier_long> The identifier reqcognised by the prefix '--'
     @param <value> The value set "after" the identifier.
     @param <arg_a_boolean> Set to true if the identifier represents a boolean.
     @return true if inserted, i.e. if the string was found.
  */
  bool insert_argument(char *_identifier_short, char *_identifier_long, char *value, bool &arg_a_boolean);
  /**! Returns true if both the short- and long identifier is not equal*/
  bool is_not_equal(class cmd_argument arg);
  //! Frees the memory rservered
  void free_mem();
  //! Prints the value set
  void print_value_set();
  /** Prints the description of the arguments if */
  char *print_cmd_argument(char *input_name, FILE *output_help, FILE *file); 
  //! Constructor (1)
  cmd_argument(char *_description, char *_identifier_short, char *_identifier_long,
	       ARGUMENT_TYPES _type, void *_arg, char *DEFAULT_OPTION_NAME, uint &DEFAULT_OPTION_NAME_COUNT);
  //! Constructor (2)
  cmd_argument(char *_description, char *_identifier_short, char *_identifier_long,
	       ARGUMENT_TYPES _type, void *_arg, char *_sub_section_name);
  //! Constructor (3)
  cmd_argument();    
#ifdef assert_code
  //! Tests the code with regard to the isenrtion of arguments, and theri validation
  static void assert_insert_argument();
#endif
  //! The main test function for this class  
  static void assert_class(const bool print_info);
};


#endif
