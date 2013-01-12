#ifndef mcl_bunch_h
#define mcl_bunch_h
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
#include "../configure.h"
#include "taxa.h"
#include <set>
using namespace std;
/**
   @file
   @class mcl_bunch
   @brief Encapsulats a string, with specific methods adding information.
   @ingroup output_format
   @remark Format according to output rules (mostly defined by the MCL-input-criteria).
   @author Ole Kristian Ekseth (oekseth)
   @date 08.01.2011 by oekseth (init).
   @date 01.09.2011 by oekseth (assert-functions).
   @date 25.12.2011 by oekseth (cleanup).
**/
class mcl_bunch {  
  const static char SEPERATOR_AFTER_ROW_START = MCL_SEPERATOR_AFTER_ROW_START;
  const static char SEPERATOR_IN_ROWS = MCL_SEPERATOR_IN_ROWS;
 private:
  uint size_string; // the reserved size of the string
  uint current_size_string; // the length/size of the string containing data
  uint position_of_line_start; // holds the position of the line start of the last started line
  uint position_of_data_start; // holds the position of where the data start of the last started line
  char *string; // the data container
  //! Enlarges the char array, and updates the size variabel
  void enlarge();
  /**! Enlarges the string  */
  void enlarge(uint size_to_be_inserted);
  //! Initally written for visual asseritng the input with regard to changing the format
  void validate();
  //! Inserts a char at the end of the data block in memory
  void set_char(char c);
  //! Adds the adjusted distance for the inparalog:    
  void insert_sim_score(float sim, float div_factor);
  //! Inserts the name as an integer.
  void set_line_start(uint world_index_in);
  //! Inserts the name: @return true if name not NULL
  bool set_line_start(char *name);
 public:
  /**
     @return true if all values were intialised
     @remarks
     - Should be de-activated during a live run.
     - Inspects a limited set of objects in this class
  **/
  bool valgrind_try_provoking_uninitialised_values() {
    bool error_found = false;
    if(position_of_data_start > size_string) return false;
    if(position_of_line_start > size_string) return false;
    //! Visits all the members in order to test that they are set.
    for(uint i = 0; i < size_string; i++) {
      if(isascii(string[i]) == false) {error_found=true;}
    }
    return error_found;
  }

  //! @return the eof of the format string.
  char *get_eof() {
    if(string) return   string + current_size_string;
    else return NULL;
  }
  /**
     @return the total number of chars used rerpesenting the header.
     @remarks
     - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
     - Makes is possible using more effective writes when MPI is activated.
  **/
  static uint get_size_of_header(char *label);
  /**
     @return the total number of chars used rerpesenting the header.
     @remarks
     - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
     - Makes is possible using more effective writes when MPI is activated.
  **/
  static uint get_size_of_header(uint index);
  /**
     @return the length of pair to be inserted.
     @remarks
     - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
     - Makes is possible using more effective writes when MPI is activated.
  **/
  static uint get_size_of_inserted_pair(uint world_index_in, uint world_index_out, float sim_score);
  /**
     @return the length of an mcl-row-element
     @remarks
     - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
     - Makes is possible using more effective writes when MPI is activated.
  **/
  static uint get_size_of_inserted_pair(uint world_index_out, float sim_score);


  /**
     @return the length of pair to be inserted.
     @remarks
     - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
     - Makes is possible using more effective writes when MPI is activated.
  **/
  static uint get_size_of_inserted_pair(char *name_in, char *name_out, float sim_score);

  /**
     @return the length of an mcl-row-element
     @remarks
     - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
     - Makes is possible using more effective writes when MPI is activated.
  **/
  static uint get_size_of_inserted_pair(char *name_out, float sim_score);

  //! Sets the header
  void set_header(char *name);
  //! Sets the header
  void set_header(uint world_index);
  /**! Prints the current status of the string; useful for dbugging purposes,
     getting a correspondence between expactaitns and known facts about whats inserted.   */
  void print_data();
  /**
     @remarks Called when no more data is to be added to the current line. Adds line end if data set
     @return the numbe of elements inserted.
  */
  uint set_line_end(const bool use_dollar_sign, const bool use_test);
  //! Called when no more data is to be added to this block
  void set_end_of_block() {set_char('\0'); }
  /**
     @brief Sets the name_index-thing
     @return the number of chars used representing it
  **/
  uint set_name_index(char *name, uint world_index_in);
  //! @return the number of chars used representing the given string
  static uint get_size_name_index(char *name, uint world_index_in);

  //! Inserts a pairwise relation using integer data
  void insert(uint world_index_in, uint world_index_out, float sim_score, float div_factor);
  //! Inserts a pairwise relation using names
  void insert(char *name_in, char *name_out, float sim_score, float div_factor);
  //! Inserts, using input of type arrKey[index_out], str_len);
  void insert(char *name, float sim_score, float div_factor);
  //! Inserts using the div-factor, i.e.const float div_factor_out = arrAvgNorm[taxon_in][taxon_out];
  void insert(uint world_index_in, float sim_score, float div_factor);
  //! Used in order to know where to copy data from
  uint get_line_start() {return position_of_line_start;}
  //! @return the last row inserted
  char *getLine() {return string + position_of_line_start;}
  //! @return the complete data set
  char *getLineData() {return string + position_of_data_start;}
  //! @return the data
  char *getData() {return string;}
  /**! @return the size of the line. Needed to store this string of data   */
  uint get_size_of_line();
  //! Asserts that the current string resides inside the line
  void et(char *s);
  /**! @return the size of the data in a line. Needed to store this string of data   */
  uint get_size_of_data() {   return (current_size_string - position_of_data_start);}
  //! @return the length of the data
  loint get_current_size_string(){return current_size_string;}
  //! Copies the lines from the arguments, assuming that only the last argument has a line end included  
  void copy_line(mcl_bunch &bunch_inpa, mcl_bunch &bunch_pair_orth, mcl_bunch &bunch_orth);
  //! @return true if a string is set
  bool has_data() {return current_size_string != 0;}
  //! Prints the data in the container
  void print();
  /**
     @return The internal variable 'string' consisting.
     @param <length> Sets the length of the string returned.
   **/
  char *get_data(uint &length) {length = current_size_string; return string;}
  /**
     @return The internal variable 'string' consisting.
   **/
  char *get_data() {return string;}
  /**
     @return The data for 'string+current_size_string' .
   **/
  char *get_end_of_buffer() {if(string!=NULL)return string + current_size_string; return NULL;}
#ifdef USE_MPI
  /**
     @brief Writes the string to the file, and clears the memory
     @return the number of chars written to the file
  */
  uint write_data_to_file(MPI_File  *file_out);
#else
  /**
     @brief Writes the string to the file, and clears the memory
     @return the number of chars written to the file
  */
  uint write_data_to_file(FILE  *file_out);
#endif

#ifdef USE_MPI
  /**
     @brief Uses MPI-lib for the write operation:
     @remarks Is the interace for writing all the data to the result files (in the final phase of the ortAgogue algorithm).
  **/
  static int write_string_to_file(MPI_File *file, char *string) {
    if(string) {
      return write_string_to_file(file, string, (uint)strlen(string));
    } else return 0;
  }
  /**
     @brief Uses MPI-lib for the write operation:
     @remarks Is the interace for writing all the data to the result files (in the final phase of the ortAgogue algorithm).
  **/
  static int write_string_to_file(MPI_File *file, char *string, uint current_size_string);
#else
  /**
     @brief Uses c-lib for the write operation:
     @remarks Is the interace for writing all the data to the result files (in the final phase of the ortAgogue algorithm).
  **/
  static int write_string_to_file(FILE *file, char *string, uint current_size_string);
  /**
     @brief Uses c-lib for the write operation:
     @remarks Is the interace for writing all the data to the result files (in the final phase of the ortAgogue algorithm).
  **/
  static int write_string_to_file(FILE *file, char *string) {
    return mcl_bunch::write_string_to_file(file, string, 0);
  }
#endif


  //! Frees the string given.
  void free_string_buffer(char *&string, uint &size_string, /*reset_internal_poitner*/ bool reset_internal_poitner);
  //! @return an allocated string vien the new size as input.
  char *allocate_string_buffer(long long int new_size);
  //! Frees memory.
  void free_mem();
  //! @return a new object of the class.
  static class mcl_bunch *init_class() {return new mcl_bunch();}
  //! Deallocates the memory of the object given.
  static void free_class(mcl_bunch *&arg);
  mcl_bunch();
#ifdef assert_code
  bool show_positive_results; // if set, shows the tests run in an texutal representation
  uint test_cnt; // if set, counts the tests
  //! Compare the internal-object-string with argument given
  bool compare_string(char *inp, const bool print_diff);
#endif
  //! The main test function for this class  
  static void assert_class(const bool show_positive_results);
};
#endif
/**
   @brief Encapsulats a string, with specific methods adding information.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
typedef class mcl_bunch mcl_bunch_t;
