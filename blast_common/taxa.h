/**
   @file
   @brief Holds data for each protein given a specific taxon.
   @ingroup blastfile_container
**/
#ifndef taxa_h
#define taxa_h
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
#include "types.h"
#include "libs.h"
//#include "tbb_libs.h"
#include "types.h"
#include "../configure.h"
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h> // some systems require it
#include <sys/stat.h>
#include <sys/termios.h> // for winsize
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "log_builder.h"
/**
   @class taxa
   @brief Holds data for each protein given a specific taxon.
   @ingroup blastfile_container
   @remark Supporting functions are included, among such:
   # Builds a linked list for storing ortholog relations, with regard to each protein.
   # Assures symmetri of the (above) linked lists.
   # Produces the mcl- and abc-output for the given ortholog relation.
   @author Ole Kristian Ekseth (oekseth).
   @date 25.08.2011 by oekseth
   @date 25.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
   @date 25.12.2011 by oekseth (cleanup).
**/
class taxa{ //! Holds data used for each protein
 private:
  // Below variables included for convenience: As the size of it is marginally compared to the whole consumption, it does not affect the overall memory signature.
  char  SEPERATOR; // the seperator to be used during the parsing
  int INDEX_IN_FILE_FOR_PROTEIN_NAME; // The index in the blast file where the protein name is found
  int INDEX_IN_FILE_FOR_TAXON_NAME; // The index in the blast file where the taxon name is found

 public:
  //! The default size of the arrOrtho
  static const uint BASE_ARR_ORTHO_STEP = 10; 
  //! variabel set in constructor, representing the total length of this taxon
  int total_cnt; 
  //! the index where those proteins who belongs to this taxa shall start at
  int rel_start; 
  //! the index fo where each protein of this taxa must be within
  int rel_end; 
  //! The name of the taxon this object represents.
  char *name; //[SIZE_TAXO];
  //! list of protein labels.
  char **arrKey; 
 private:
  //! Holds the list of the overlap values for the proteins.
  uint *arrOverlap;
 public:
  
  //!  @return the list of overlap values for the given taxon:
  uint *getArrOverlap(){return arrOverlap;}
  
  //! Holds the list of the inparalog limits for the proteins.
  float *arrInpaLimit;

  //!  @return the list of inparalog limits for the given taxon:
  float *getArrInpaLimit(){return arrInpaLimit;}
  //!  @return the inparalog limit for the given taxon:
  float get_arrInpaLimit(int key) {
    if(!arrInpaLimit) return 0;
    assert(key < total_cnt);
    return arrInpaLimit[key];
  }
  /**! @return the mcl header, having a very long string as the input */
  char *mcl_get_all_header_proteins(char *string);  
  //! Sets the taxon name, and set internal position of the string to its end
  //! @remarks: increments the pointer.
  void setTaxonName(char *&string, const bool increment_pointer);
 public:
  //! Inserts the overlap list without allocting memory
  void set_unsafe_overlap_list(uint *&over);
  /**
     Copies the overlap list into this
     @remarks Assumes that if the list in this object is set, it's set tot the proper length.
  **/
  void copy_arr_overlap_into_list(uint *over);

  /**
     Copies the inpaLimit list into this
     @remarks Assumes that if the list in this object is set, it's set tot the proper length.
  **/
  void copy_arr_inpalimit_into_list(float *inpalist);

  /**
     @return The total overlap score
     @remarks: Useful function for validating after an update
  **/
  uint get_total_overlap_score();
  /**
     @return The total overlap score
     @remarks: Useful function for validating after an update
  **/
  static uint get_total_overlap_score(taxa *&listTaxa, uint index_start, uint index_end_pluss_one);
  //! @return the upper estimate for the object size
  uint get_object_size(uint taxon_length);
  //! @return an approximation for the memory consumption.
  uint get_object_memory_size();
   /**
     @brief Prints the memory consumption for the object.
     @remarks For anaylsing the memory fingerprint.
  **/
  static void print_getTotal_memoryConsumption_for_object(taxa *listTaxa, uint taxon_length, FILE *f);

  /**
     @brief Sets the protein label into the list, and updates the internal position of the string to its end
     @param <local_index> If index is out of bound, the input string is set to NULL.
     @param <string> Set to NULL if index not found.
  **/
  void setProteinName(uint local_index, char *&string);
  //! @returns the complete protein name, to be used in the mcl output
  char *getCompleteProteinName(uint local_index);
  //! @returns the relative start position, using internal variable 'rel_start'.
  int getRelativeStartingIndex() {return rel_start;}
  //! @returns the relative end-position, using internal variable 'rel_end'.
  int getRelativeEndIndex() {return rel_end;}
  //! @returns the absolute end-position.
  int getAbsoluteEndIndex() {return rel_end+total_cnt;}
  //! @returns the absolute start-position.
  int getAbsoluteStartIndex() {return rel_start+total_cnt;}
  //! Takes a key from teh world index and prints the chars representing it
  void printProtein(int real_key, const bool print_taxon, char seperator) {
    printProtein(stdout, real_key, print_taxon, seperator);
  }
  //! Takes a key from teh world index and prints the chars representing it
  void printProtein(FILE *f, int real_key, const bool print_taxon, char seperator);
  //! @returns a local key to be used in this subset of the data set
  uint getLocalIndex(uint world_key);
  //! @returns the unique key tbu in the whole data set
  uint getWorldIndex(uint local_key);
  //! @returns the length of the protein space (total number of proteins)
  uint getLength() {return (uint)total_cnt;}
  //! @returns true if the local key is in this range
  bool is_a_local_key(int protein_index) {
    return ((protein_index >-1) && (protein_index < total_cnt));
  }
  //! @returns true if the world key given as input belongs to this taxon
  bool isRelation(int world_key);
  //! Initiates the *arrInpaLimit
  void initArrInpaLimit(bool only_if_empty);
  //! @return true if protein is accepted:
  bool aboveInparalogLimit(uint protein_id, const float averaged_similarity, const bool INPARALOG_OPERATION);
  /**
     @return An arrOrtho row built.
     @remark Used in the process of bulidng a two-way-relation (setting the pairs)
     not malloced before at this point (after the iteration of the orthologs,
     but before the iteration of the inparalogs).
  */
  uint* initOrthoRow(uint size);
  /**
     @brief Updates in 'sim_score' arbuemnt is higher than the
     previous value set at position 'protein_in'
  */
  void insertInpaLimit(const uint protein_in, float sim_score);
  /**
     @param <real_key_out> the key of the world outer index
     @param <real_key_in> the key of the world inner index
     @return True if the 'arrOrtho' relation exsists.
  */
  bool is_ArrOrthoRelation(uint real_key_out, uint real_key_in);
  /**   Builds a two way relation for a given row  */
  void buildTwoWayArrOrtho_row(const uint protein_id, taxa *listTaxa, const uint taxon_length);
  /**
     @brief Builds a two way relation
     @warning Not thread safe.
  **/
  void buildTwoWayArrOrtho(const uint taxon_length, taxa *listTaxa);
  /** Deallocates 'arrInpaLimit'   */
  void delete_initArrInpaLimit();
  /**!  Frees 'arrInpaLimit' Used before collecting the data to an output file
   */
  void delete_initArrOrtho();
  /**
     @param <arrOrtho> the ortho array
     @param <protein_in> the index to get the ortholog-matches for.
     @param <&size> Sets the size of the array returned
     @return A list (array) of the proteins 'protein_in' has ortholog pairs with.
  */
  uint *getOrthoRow(uint **arrOrtho, uint protein_in, uint &size);  
  /** Print the list of orthologs based on the array given as input    */
  void printArrOrthoList(uint **arrOrtho, const uint start, uint end, const bool use_names, taxa *listTaxa, uint taxon_length);
  /**
     @brief Frees the memory allocated.
     @remark Is to be done at the end of the program.
  */
  void delete_buffer();
  //! Prints the overlap values.
  void printOverlap(const bool use_name);
  //! Prints the variables set in this object.
  void print_variables();
  //! Prints the variables set in this object.
  void print_variables(FILE *F);
  //! @returns the overlap using the local index
  uint getOverlap(uint local_index);
  //! @returns the inparalog limit using the local index
  float getArrInpaLimit(uint local_index);
  //! Prints the limits (i.e. the similarity scores) whom the inpralogs must be above
  void printArrInpaLimit(bool use_names);
  //! Prints hte key name pair for this taxon
  void printArrKeyList(int taxon_index);  
  //! @returns the taxon id based on the inputs given
  static uint getTaxonId(uint world_key, taxa *listTaxa, uint taxon_length);
  //! @returns the taxon name of this object.
  char *getName() {return &name[0];}
  //! @returns the name of the protein, given as a world index
  char *getProteinNameWorld(uint world_ind_in);
  /**
     @return the name of the protein, given as a local index
  */
  char *getProteinName(uint local_ind_in) {return arrKey[local_ind_in];}
  /**
     @return the list of protein labels.
  */
  char **getArrKey() {return arrKey;}
  /**     @return the protein label for the given argument id.  */
  char *getArrKey(uint protein_id);
  /**
     @return the protein label for the given argument id.
     @param <world_protein_id> The index given using the global id scheme.
  */
  char *getArrGlobalKey(int world_protein_id);

  /**
     @brief Sets the name of the taxon.
     @remark To be executed in the first (of two) parsing operations through the input blast file.
  */
  void setName(char *_name);
  //! Frees the memory
  void free_mem() {delete_buffer();}
  /**
     @param <size> The number of objects of type taxa to provide.
     @return an initialized list of objects with the length provided
  **/
  static taxa *init(int size) {return new taxa[size]();}

  /**
     @brief De-allocates only the pointer to the class wrapper, and not the content itself
     @remarks Used after parallel processing.
  **/
  static void close_only_pointer(taxa *&obj) {
    if(obj) {delete [] obj; obj = NULL;}
  }
  //! Deallocates the memory for the given object set.
  static void delete_taxa(int &taxon_size, taxa *&list);
  //! Deallocates the memory for the given object set.
  static void close(taxa *&list, int &taxon_size) {
    delete_taxa(taxon_size, list);
  }

  /**
     @brief Initiates the **arrOrtho
     @remark Is first needed in the 'ortho_set_bin', therefore
     not malloced before at this point (after the iteration of the orthologs,
     but before the iteration of the inparalogs).
     @param <def_size> if set to true, intialize the rows: if not set, then the rows are noe intiated
  */
  void initArrOrtho(bool def_size);
  /**
     @brief Initiates a row of the '**arrOrtho
     @remark Is first needed in the 'ortho_set_bin', therefore
     not malloced before at this point (after the iteration of the orthologs,
     but before the iteration of the inparalogs).
  */
  void initArrOrthoRow(uint protein_in, uint size);
    /**
     @brief Generate (writes) a log file describing the details of the memory signature for taxa collected during parsing
     @remarks Useful for analyzing- optimizing algorithmic correctness- and memory allocation procedures.
  **/
  static void generate_memory_allocation_overview(char *file_name, taxa *listTaxa, int taxon_length);

  //! Writes the information about the class to the given given param.
  void print_class_info(FILE *out);  
  //! Writes the information about the object to the given given param.
  static void print_class_info(FILE *out, taxa *listTaxa, int taxon_length);
  //! The Constructor.
  taxa();
  //! The Constructor.
  taxa(char *_name, int _total_cnt, int _rel_start, int _rel_end, char **_arrKey, char sep, int index_protein, int index_taxon);
  //! Purpose of this code is to provide a default length of each taxon to be used
  static taxa * init_test_array(uint taxon_size, uint protein_size);
  //! @return a char pointer holding the taxa names, given a taxa object.
  static char **get_taxa_names(taxa *obj, int taxon_length);
  //! Sets the proteins above the given threshold.
  static void get_cnt_taxons_above_limit(taxa *obj, int taxon_length, int protein_limit,  loint &cnt_proteins_temp, loint &cnt_taxa_temp);
  /**
     @brief Used validating/comparing results, ie, before- and after MPI transfer.
     @return the total score
  **/
  static float get_total_score_for_arrInpaLimit(taxa *&listTaxa, uint taxon_length);
#ifdef USE_MPI
/**
   @brief Sends- and receives the list of arrInpaLimit* objects accross nodes.
**/
  static void mpi_send_and_receive_arrInpaLimit(taxa *&listTaxa, uint taxon_length);
/**
   @brief Sends- and receives the list of norm_t* objects accross nodes.
**/
  static void mpi_send_and_receive_arrOverlap(MPI_Comm comm, int myrank, int number_of_nodes, taxa *&listTaxa, uint index_start, uint index_end_pluss_one);
#endif


#ifdef assert_code
 private:
  //! Tests methods for bulding of reciprocal ortholog relations.
  void assert_arrOrthoFunctions();
 public:
#endif
  //! The main test function for this class
  static void assert_class(const bool show_positive_results);
 };

/**
   @brief Holds data for each protein given a specific taxon.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth).
   @date 25.12.2011 by oekseth (cleanup).
**/
typedef class taxa taxa_t;
#endif
