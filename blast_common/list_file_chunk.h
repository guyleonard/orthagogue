#ifndef list_file_chunk_h
#define list_file_chunk_h
/**
   @file
   @brief Holds estimations for the file chunks sized.
   @remarks To reduce complexity when validating result.
   - For simplification assumes that the standard format is used:
   - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
   - Makes is possible using more effective writes when MPI is activated.
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
/**
   @class list_file_chunk
   @brief Holds estimations for the file chunks sized.
   @remarks To reduce complexity when validating result.
   - For simplification assumes that the standard format is used:
   - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
   - Makes is possible using more effective writes when MPI is activated.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
class list_file_chunk {
 private:
  bool MODE_PAIRWISE_OUTPUT_ABC; // for the abc file: if set, the data out pairwise in stead of as in a row
  bool MODE_PAIRWISE_OUTPUT_MCL; // for the mcl file: if set, the data out pairwise in stead of as in a row
  // The following variables decides what output is to be printed
  bool PRINT_IN_ABC_FORMAT;
  bool PRINT_IN_MCL_FORMAT;
  uint *list; // Hold the number of chars included for each of the files written.
 public:
  //! @return a pointer to the list:
  uint *get_list() {return list;}
  //! @return the size (length) given the index;
  uint get_length(uint index) {return get_size(index);}
  //! @return the size (length) given the index;
  uint get_size(uint index);
  //! Compare lists using the first list as basis, and assuming a length defined in file enum_mcl.h
  bool compare_with_list(loint *two);
 private:
  //! @return the file size, given the index:
  uint get_file_size(uint index, char *FILE_BINARY_LOCATION);
 public:
  /**
     @brief 
     @param <file_chunk_header> Combines the internal list with this parameter to produce the size for comparison.
     @param <FILE_BINARY_LOCATION> The folder where the temporary files are stored.
     @return true if all sizes corresponds to expected file sizes.
  **/
  bool compare_sizes_with_files(list_file_chunk *file_chunk_header, char *FILE_BINARY_LOCATION);
  //! @return an initiated version of the array of containers holding the number of chars read.
  uint *init_list();
  //! Inserts an element into the list given.
  void insert_size(uint index, uint value);
  //! Resets the value of an element into the list given.
  void reset_size(uint index, uint value);
  /**
     @brief Merges two lists.
     @param <src_obj>  The source of the data. Leaves it un-tuched.
     @remarks Does not copy data if src not set.
  **/
  void merge_lists(list_file_chunk *src_obj);
  //! @return the expected sice of the mci header for an mci file.
  static uint get_expected_size_of_header(taxa_t *listTaxa, int taxon_length);
  /**
     @brief appends the header sizes to the list.
     @return a checksum corresponding to the total number of chars inserted
  **/
  uint append_header_sizes(taxa_t *listTaxa, int taxon_length);
  //! @brief sets the size consumption for the header, given the settings:
  void get_size_of_header(char *label, uint world_index, uint &size_string, uint &size_number);
  //! Sets the values for the header using properties of this object.
  void insert_sizes_header(uint world_index, char *name, uint type_string, uint type_number);
  /**
     @return the length of pair to be inserted.
     @remarks
     - For simplification assumes that the standard format is used:
     - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
     - Makes is possible using more effective writes when MPI is activated.
  **/
  void get_size_of_inserted_pair(uint world_index_out, char *label_in, char *label_out, float sim_score, uint &size_string, uint &size_number);

  /**
     @return the length of pair to be inserted.
     @remarks
     - For simplification assumes that the standard format is used:
     - Used for making the lengths of the file calculatable before- and during- and after the operation, validating the process.
     - Makes is possible using more effective writes when MPI is activated.
  **/
  void get_size_of_inserted_pair(uint world_index_out, char *label_out, float sim_score, uint &size_string, uint &size_number);
  //! Sets the values for the header using properties of this object.
  void insert_sizes_for_a_pair(uint world_index_out, char *label_in, char *label_out, float sim_score, uint type_string, uint type_number);

  //! Sets the values for the header using properties of this object.
  void insert_sizes_for_a_pair(uint world_index_out, char *label_out, float sim_score, uint type_string, uint type_number);
  //! De-allocates the memory reserved.
  void free_memory();

  //! The constructor:
  list_file_chunk();
  //! The constructor:
  list_file_chunk(
		  bool _MODE_PAIRWISE_OUTPUT_ABC,  // for the abc file: if set, the data out pairwise in stead of as in a row
		  bool _MODE_PAIRWISE_OUTPUT_MCL,  // for the mcl file: if set, the data out pairwise in stead of as in a row
		  // The following variables decides what output is to be printed
		  bool _PRINT_IN_ABC_FORMAT, 
		  bool _PRINT_IN_MCL_FORMAT
		  );

};
#endif
