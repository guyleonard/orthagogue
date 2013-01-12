#ifndef mcl_format_h
#define mcl_format_h
/**
   @file
   @brief Produces a set of data string for mcl-like output
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
**/
#include "mcl_bunch.h"
#include "enum_mcl.h"
#include "mcl_print_settings.h"
#include "list_file_chunk.h"
#include "../blast_filtering/pipe_struct_result.h"
/**
   @class mcl_format
   @brief Produces a set of data string for mcl-like output
   @ingroup output_format
   @remark: Uses enums of type "mcl_t", with the corresponding variables
   holding access data for this class.
   - Updates (adding of new methods) of this class requires uppdates of this structure
   as well, as those are closely bundles together.
   (Therefore the enum is located above).
   - In order to do print outs of this class, static methods are provided.
   @author Ole Kristian Ekseth (oekseth)
   @date 08.01.2010 by oekseth (initial)
   @date 26.08.2010 by oekseth (asserts)
   @date 25.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
   @date 25.12.2011 by oekseth (cleanup).
**/
class mcl_format {
 private:
  int taxon_length;
  taxa_t *listTaxa;
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT;
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT;
  bool MODE_PAIRWISE_OUTPUT_ABC; // for the abc file: if set, the data out pairwise in stead of as in a row
  bool MODE_PAIRWISE_OUTPUT_MCL; // for the mcl file: if set, the data out pairwise in stead of as in a row
  // The following variables decides what output is to be printed
  bool PRINT_IN_ABC_FORMAT;
  bool PRINT_IN_MCL_FORMAT;
  bool SORT_ABC_DATA; // if true, sorts teh abc files before outprint
  mcl_bunch_t **string;
  set<uint> set_orth_inpa; // holds the list of orth_inpa for each row.
  class list_file_chunk *file_chunk;
 public:
  //! @return the file of file chunks.
  list_file_chunk *get_file_chunk() {return file_chunk;}
 private:  
  /**
     @return true if all values were intialised
     @remarks
     - Should be de-activated during a live run.
     - Inspects a limited set of objects in this class
  **/
  bool valgrind_try_provoking_uninitialised_values() {
    bool error_found = false;
    for(uint i = 0; i < mcl_t_size; i++) {
      if(string[i]->valgrind_try_provoking_uninitialised_values()) {
	error_found = true;
      }
    }
    return error_found;
  }
  //! Sets the values for the header using properties of this object.
  void set_debug_header_size(uint world_index, char *name, uint type_string, uint type_number);
  //! Sets the values for the header using properties of this object.
  void set_debug_pair_size(uint world_index_out, char *label_in, char *label_out, float sim_score, uint type_string, uint type_number);
  //! Sets the values for the header using properties of this object.
  void set_debug_pair_size(uint world_index_out, char *label_out, float sim_score, uint type_string, uint type_number);


  //! Initiatlizes the mcl-bunch-classes:
  void init_string();
  //! Frees the memory reserved for the strings
  void free_string();
 public:
  //! Updates the result cheme for use in the log:
  void log_update_result_cheme(pipe_struct_result &result);
  //! Frees the memory reserved for the strings
  void free_mem() {free_string();}
  //! Prints the current status
  void print_size_of_data(); 
  //! Inserts the headers for inparalogs.
  void set_header_inpa(uint world_index, char *name); 
  //! Inserts inparalogs: Does not check the uniqueness. (The data source to be used has unique values)
  void insert_inpa(uint world_index, char *name, float sim_score, float div_factor);
  //! Inserts inparalogs for pairwise output: Does not check the uniqueness. (The data source to be used has unique values)
  void insert_inpa(uint world_index_in, uint world_index_out, char *name_in, char *name_out, float sim_score, float _div_factor);
  //! Sets the name_index-thing
  void set_name_index(char *name, uint world_index);
  //!     Sets the  header for ortho pairs.
  void set_header_pair_ortho(char *name, uint world_index);
  /**
     @brief Inserts ortho pairs: Does not check the uniqueness. (The data source to be used has unique values)
     @return '1' if an eleemnt was inserted, else '0'.
     @remarks The return value is used in debug mode in class pipe_struct verifying that the number of inserted realtions correspons to our expectations.
  **/
  void insert_pair_ortho(uint world_index_in,uint world_index_out,char *name_in,char *name_out, float sim_score, float _div_factor);
  /**
     @brief Inserts orthologs into both the 'name' and 'number' strings: Checks uniqueness.
     @return '1' if an eleemnt was inserted, else '0'.
     @remarks The return value is used in debug mode in class pipe_struct verifying that the number of inserted realtions correspons to our expectations.
  **/
  uint insert_ortho_inpa(char *name, uint world_index, float sim_score, float _div_factor);
  //! Inserts orthologs into both the 'name' and 'number' strings: Checks uniqueness.
  uint insert_ortho_inpa(uint world_index_in,uint world_index_out,char *name_in, char *name_out, float sim_score, float div_factor);
  //! Inserts header for orth_inpa
  void set_header_ortho_inpa(char *name, uint world_index);
  //! @return the current size of the line
  int get_size_of_line(mcl_t type) {if(string!=NULL) return (int)string[type]->get_size_of_line(); return -1;}
  //! @return the current size of the data
  int get_size_of_data(mcl_t type) {if(string!=NULL) return (int)string[type]->get_current_size_string(); return -1;}
  /**! Sets the line end, and updates the 'all' and 'all_number' strings */
  void set_line_end();
  //! @return the eof of the current string using integer representation:
  char *get_ortholog_current_eof(const bool type_integer) {
    if(string) {
      if(type_integer) {if(string[pair_orth_number]) return string[pair_orth_number]->get_eof();}
      else if(string[pair_orth]) {return string[pair_orth]->get_eof();}
    }
    return NULL; // If non of the above matched, we return this.
  }

  /**
     @return the file names
     @remarks The naming convention follows:
     - Names with integers (numbers in this context) corresponds to mcl type, thereby .mcl'
     - Those files with strings (those other in this contect), the naming covention '.abc' is used
  **/
  static char *get_file_name(uint id);
#ifdef USE_MPI
  /**! Opens the file given the identifier, describing the purpose of the file   */
  static MPI_File *file_open(char *identifier, char *FILE_BINARY_LOCATION);
  /**! @return a list of open files */
  static MPI_File ** init_file_array(mcl_t type, char *FILE_BINARY_LOCATION, bool PRINT_IN_ABC_FORMAT, bool PRINT_IN_MCL_FORMAT);
#else
  /**! Opens the file given the identifier, describing the purpose of the file   */
  static FILE *file_open(char *identifier, char *FILE_BINARY_LOCATION);
  /**! @return a list of open files */
  static FILE ** init_file_array(mcl_t type, char *FILE_BINARY_LOCATION, bool PRINT_IN_ABC_FORMAT, bool PRINT_IN_MCL_FORMAT);
#endif
#ifdef USE_MPI
  /**! Closes a list of open files */
  static void close_init_file_array(MPI_File **&file, mcl_t type, bool SORT_ABC_DATA, char *FILE_BINARY_LOCATION);
  /**
     @brief Creates the headers, and writes them to the files 
     @return a checksum corresponding to the total number of chars inserted
  */
  static uint write_file_headers(MPI_File **file, taxa_t *listTaxa, int taxon_length, bool PRINT_IN_MCL_FORMAT);
  /**! Writes the data to the files, given the data set. */
  void write_file(MPI_File **file);
#else
  /**! Closes a list of open files */
  static void close_init_file_array(FILE **&file, mcl_t type, bool SORT_ABC_DATA, char *FILE_BINARY_LOCATION);
  /**
     @brief Creates the headers, and writes them to the files
     @return a checksum corresponding to the total number of chars inserted
  */
  static uint write_file_headers(FILE **file, taxa_t *listTaxa, int taxon_length, bool PRINT_IN_MCL_FORMAT);
  /**! Writes the data to the files, given the data set. */
  void write_file(FILE **file);
#endif
  //! @return a new objectof this class.
  static mcl_format *init(mcl_print_settings_t _print_settings) {  return new mcl_format(_print_settings);}
  //! @return a new objectof this class.
  static mcl_format *init(int taxon_length, taxa_t *listTaxa, bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT, 
			  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT, bool MODE_PAIRWISE_OUTPUT_ABC,  
			  bool MODE_PAIRWISE_OUTPUT_MCL, bool PRINT_IN_ABC_FORMAT, bool PRINT_IN_MCL_FORMAT, 
			  bool SORT_ABC_DATA);

  //! Destructor for the class given an object as input
  static void free_class(mcl_format *arg) {arg->free_string(); delete arg; arg = NULL;}
  //! Class constructor:
  mcl_format(  int _taxon_length,
	       taxa_t *_listTaxa, 
	       bool _DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT, 
	       bool _DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT, 
	       bool _MODE_PAIRWISE_OUTPUT_ABC,  // for the abc file: if set, the data out pairwise in stead of as in a row
	       bool _MODE_PAIRWISE_OUTPUT_MCL,  // for the mcl file: if set, the data out pairwise in stead of as in a row
	       // The following variables decides what output is to be printed
	       bool _PRINT_IN_ABC_FORMAT, 
	       bool _PRINT_IN_MCL_FORMAT, 
	       bool _SORT_ABC_DATA  // if true, sorts teh abc files before outprint)
	       );
  //! Class constructor:
  mcl_format(mcl_print_settings_t print_settings);

#ifdef assert_code
  //! Prints the elements found in the 'orth_inpa list
  void print_ortho_inpa_line_list();
  //! Prints the content of the string given its type.
  void print(mcl_t type) {string[type]->print();}
  //! Compares the strings
  bool compare_string(char *inp, const bool print_diff, mcl_t type);
#endif

  //! The main test function for this class  
  static void assert_class(const bool show_positive_results);
}; 
#endif
/**
   @brief Produces a set of data string for mcl-like output
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
typedef class mcl_format mcl_format_t;
