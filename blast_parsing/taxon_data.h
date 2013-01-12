#ifndef taxon_data_h
#define taxon_data_h
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
#include <algorithm>
#include "log_builder.h"
using namespace std;
/**
   @brief The default size of the buffer for class taxon_data.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 16.09.2011 by oekseth (asserts)  
**/
static const uint START_SIZE_BUFFER = 3;  

/**
   @file
   @class taxon_data
   @brief Holds a list of names (i.e. sometimes called 'labels' or 'strings').
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)  
**/
class taxon_data {  
 public:
  //! The total number of elements in the buffer
  uint prots_reserved; 
  //! The number of proteins whos room for
  uint prots_used; 
 private:
  char **buffer;
  char *taxon_name;
  char *prev_protein_name; // Some proteins have multiple defention: this variable tries to avoid the duplicaiton of keys
  //! Initiates the 'taxon_name' to at set of '\0' chars.
  //  void init_taxon();
 public:
  //! @return true if object is set.
  bool has_data(){return ( (buffer!=NULL) && (taxon_name != NULL));}
  /**
     @brief Inserts the given object into this.
     @param <obj> The object to get variables from, and remove references from.
  **/
  void steal_object(taxon_data &obj) {
    assert(!prots_reserved);
    assert(!buffer);
    assert(!taxon_name);
    assert(!prev_protein_name);
    prots_reserved = obj.get_prots_reserved();
    prots_used = obj.get_prots_used();
    buffer = obj.get_buffer();
    taxon_name = obj.get_taxon_name();
    prev_protein_name = obj.get_prev_protein_name();
    obj.remove_references();
  }

  /**
     @brief Removes the references, without deallocating the memory.
     @remarks Used in conjunction with stealing an object.
  **/
  void remove_references() {
    prots_reserved = 0;
    prots_used     = 0;
    buffer =            NULL;
    taxon_name =        NULL;
    prev_protein_name = NULL;
  }
  //! Sets data for 'this'
  void set_object(uint reserved, uint used, char **&_buffer, char *&_taxon_name, char *&prev) {
    prots_reserved = reserved, prots_used = used, buffer = _buffer, taxon_name = _taxon_name, prev_protein_name = prev;
    //    strcpy(taxon_name, _taxon_name);    
  }
  //! Prints data about this object.
  void print_data(uint index);
  //! @return the name given the param.
  char *get_name(uint index);
  //! @return the number of proteins marked in this object.
  uint get_prots_used(){return prots_used;}
  //! @return the number of proteins marked in this object.
  uint get_length() {return prots_used;}
  //! @return the number of proteins reserved for this object.
  uint get_prots_reserved(){return prots_reserved;}
  //! @return the buffer.
  char **get_buffer() {return buffer;}
  //! @return the taxon name.
  char *get_taxon_name() {return taxon_name;}
  //! @return the previous protein name found for this object.
  char *get_prev_protein_name() {return prev_protein_name;}
  
  //! @return true if input name is the same as the taxon-name of this object.
  bool is_equal(char *name);
  //!  @return True if the strings given are equal.
  bool compare_strings(char *name_1, char *name_2) {
    if(name_1 && name_2) {
      return (0 == strcmp(name_1, name_2));
    } else { // If both are not set, they are equal.
      if((name_1 == NULL) && (name_2 == NULL)) return true;
      else return false; // only one of them is set ie, they are not equal.
    }
  }

  //! Sets the taxon name
  void setTaxonName(char *taxon);
  //! Sets the taxon name
  void set_taxon_name(char *taxon) {setTaxonName(taxon);}
  //! Sets the parm with value of this object.
  void get_taxon_name(char *&name, loint &size_arg);
  //! @return the first inserted protein
  char *getFirstInsertedProtein(){return buffer[0];}
  //! @return the last inserted protein: If empty, returns NULL
  char *getLastInsertedProtein();
  //! @return the number of proteins reserved with regard to the storage (list)
  uint getProteinsReserved(){return prots_reserved;}

  //! @return The size of the proteins
  uint getProteinsUsed(){return prots_used;}

  //! @return The buffer of the protein names
  char **getBuffer(uint &size) {
    printf("___ prots_used(%d) and buffer_is_set(%d) at line %d in file %s\n", size, buffer!=NULL, __LINE__, __FILE__);
    size = prots_used; return buffer;
  }

  //! Prints the buffer
  void printf_buffer();
  /**
     @brief The 'main method' of this structure: 
     @param <name> Inserts it into the the list.
     @return The pointer to the name given as input argument.
     @remark If buffer is full, enalarges the list.
     @date 01.12.2010
  **/
  char* insert_protein(char *name);
  //! Inserts the label at specific index
  void insert_label(char *name, uint index);
  /**
     @brief Deletes the buffer holding the indexes
     @deprecated For backward compability,
  */
  void delete_buffer();

  //! Merges the argument 'taxon_data' with 'this'
  void mergeData(taxon_data data, bool overlap);
  
  /**
     Enlarges the buffer if the data to add is larger than the free space available
     @param <min_size> The minimum size the number of elements the buffer must consist of.
  **/
  void enlarge_buffer(loint min_size);
  //! @return true if the identifier is the same as the argument (given as input).
  bool is_taxon(char *arr_taxon_name);
  /**
     @param <size> The size of this object not initialized.
     @return A list of this object
  */
  static taxon_data *init_empty(uint size);
  /**
     @param <size> The size of this object initilazed.
     @return A list of this object
  */
  static taxon_data *init(uint size);
  /**
     @brief Closes the object, de-allocating the memory.
     @param <obj> The list of objects to close
  */
  static void close(taxon_data *&obj) {
    // TODO: why is "obj.delete_buffer()" not called??
    delete [] obj; obj = NULL;
  }
  /**
     @brief Closes the list of objects, de-allocating the memory.
     @param <list> The list of objects to close
     @param <size> The size of this object initilazed.
  */
  static void close(taxon_data *&list, uint size) {
    if(list) {
      for(uint i = 0; i < size; i++) list[i].delete_buffer();
      delete [] list, list = NULL;
    }
  }
  /**
     @brief Puts the data into the taxon_data object, and sets the 'this' to empty.
     @remarks After this operation is called, deallocating this object is safe, ie, not conflicting with the recently added data.
  **/
  void take_and_reset_data(taxon_data &obj) {    
    // Sends 'my' data to the argument:
    obj.set_object(prots_reserved, prots_used, buffer, taxon_name, prev_protein_name); 
    // Resets own data:
    prots_reserved=0, prots_used=0, buffer = NULL, taxon_name=NULL, prev_protein_name=NULL;     
  }
  //! Copies protein name from 'src' to 'dest'.
  static void copy_protein_name(char *&dest, char *src);
  //! Copies taxon name from 'src' to 'dest'.
  static void copy_taxon_name(char *&dest, char *src);

#ifdef USE_MPI
  /**
     @brief Sets the label of proteins
     @param <protein_length> The number of proteins to use as input.
     @param <listOfProteins> The 1d-list holding the protein names.
     @param <listOfLengths>  The list holding the length of each protein.
     @return The number of chars used for the names inserted.
     @remarks Important is to use the correct offsets when calling these functions!
  **/
  uint set_labels_using_1d_object(uint protein_length, char *listOfProteins, uint *listOfLengths);

  /**
     @brief Compares in order to verify the equalness.
     @param <obj> The object to compare with this object.
     @param <verbouse> Exaplains all the differences found.
     @return true if they are equal.
  **/
  bool compare_object(taxon_data &obj, bool verbouse);

  /**
     @brief Inserts the taxon-name into the list.
     @remarks In order to easyly be sendt with the mpi-procedure, each name is ended with an '\0' symbol.
  **/
  void insert_the_name_into_argument(char *&taxonNamesArr, uint &taxonNamesArr_size, uint &taxonNamesArr_used);
  /**
     @brief Asserts that the the list is correctly transformed.
     @remarks Only active when DBEBUG macro varaible is not defined.
  **/
  void assert_1d_list(uint bufferLength_startpos, uint &start_position_in_buffer, char *bufferTmp, uint *bufferLength);

  //! Transform "**buffer" into "*bufferTmp" and "int *bufferLength". Is validated internally. Used for MPI
  //  void transform_protein_label_buffer_into_1d_list(char *&bufferTmp, uint &bufferTmp_size, uint *&bufferLength,uint &bufferLength_size);
  void transform_protein_label_buffer_into_1d_list(char *&bufferTmp, uint &bufferTmp_char_used, uint &bufferTmp_size,uint *&bufferLength,uint &bufferLength_size, uint bufferLength_startpos);
#endif

  //! Deletes the buffer holding the indexes
  void free_mem() {delete_buffer();}
  //! The constructor: by default it does not reserve memory
  taxon_data();
  //! @brief The constructor. @param <reserve_size> The length of the data to reserve memory for.
  taxon_data(loint reserve_size);

  //! The main test function for this class  
  static void assert_class(const bool show_positive_results);
#ifdef assert_code
 private:  // The following code contains test function for this class:  
  //! Testing the merging of two data strctures
  void assert_mergeData();
  //! Asserting that the taxon name is set correct, using the functions for inserting, and retrieval
  void assert_setTaxonName(taxon_data &list);

  //! Testing that the isnertin of proteins are correct:
  void assert_insert_protein();
#endif
};

/**
   @brief Holds a list of names (i.e. sometimes called 'labels' or 'strings').
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
   @date 16.09.2011 by oekseth (asserts)  
**/
typedef taxon_data taxon_data_t;

#endif
