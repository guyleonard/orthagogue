#ifndef blast_parsing_h 
#define blast_parsing_h
/**
   @file
   @brief Handling a blast file, dumping them into retrievable containers.
   @ingroup parsing
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
#include "../configure.h"
#include "log_builder.h"
#include "norm_t.h"
#include "list_norm.h"
#include "pipe_parse_file.h"  
#include "pipe_parse_parse.h"
#include "pipe_parse_merge.h"   
#include "cmd_list.h"
#include "bp_container.h"
#include "tsettings_input.h"
/**
   @class blast_parsing
   @brief Handling a blast file, dumping them into retrievable containers.
   @ingroup parsing
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by Ole Kristian Ekseth (init)
   @date 18.08.2011 by Ole Kristian Ekseth (Cleaning.)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/
class blast_parsing {
 private:
  // Below are config-settings for this scope of the processing.
  // --------------------
  //! If set, prints the data used as basis for the normalization procedure.
  bool DEBUG_NORM; 
  //! If set, prints the data used for normalization of the similarity scores.
  bool PRINT_NORMALIXATION_BASIS; 
  //! If set, uses all realtions in the input as basis for the normalization preocedure. Alternative, only those values filtered out, is to be used.
  bool USE_EVERYREL_AS_ARRNORM_BASIS; 
  bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE; // If set, uses the last blast bolumn as input for data, instead of the second last.
  bool USE_BEST_BLAST_PAIR_SCORE; // If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
  //! The numbers of cpu's to build the paralisation for: only approximate value is neccassary to be set.
  int CPU_TOT; 
  int INDEX_IN_FILE_FOR_PROTEIN_NAME; // The index in the blast file where the protein name is found
  int INDEX_IN_FILE_FOR_TAXON_NAME; // The index in the blast file where the taxon name is found
  lint MAX_PARSE_BUFFER_SIZE; // TODO: move to configure-file.
  //! An initial guess on how many columns there are in the blast file with regard to the name-label.
  uint DEFAULT_NUMBER_OF_COLUMNS_IN_NAME;
  char SEPERATOR; // The seperator to be used during the parsing
  //! The path to store the folders whom the binary data will be located at
  char *FILE_BINARY_LOCATION; 
  char *FILE_INPUT_NAME; // the name (the whole path) of the input file
  const static bool use_improved_overlap_algo = true;
  // Below are data containers for this scope of the processing.
  list_norm_t *arrNorm; // Container storing the normative values
  taxa_t *listTaxa; // Container storing info about each taxonfound.
  int taxon_length; // The number of taxa in the collection.
  list_file_parse_t *listParseData; // Holds the blast file in a modified form.
  float max_input_value; // The maximum similairty score found in the blast file used.
  uint LIMIT_MINIMUM_NUMBER_OF_PROTEINS_FOR_EACH_TAXA;
  uint disk_buffer_size;
  bool DEBUG_print_pairs_in_file_parse_log_file;
 public:
/** 
    @brief Updates the command line interface with values to be set by the user:
    @param <cmd> The input list to be made available from the terminal console.
    @parm <first_pass> In order to get the correct order of the fields, i.e. "OUTPUT" before "OPERATIONAL".
 **/
  void init_values(cmd_list *cmd, const bool first_pass=true);
  
  /**
     @brief Gives the user settings for the blast parsing process.
     @param <obj> The object of class tsettings_input to update the settings found in this object with.
     @remarks The filtering process will depend upon some of these settings, therby relevant.
   **/
  void get_input_settings(tsettings_input_t &obj);
  //! @return True if the normalization process shall be debugged during the running of the software.
  bool get_debug_norm(){return DEBUG_NORM;}
  //! @return True if the basis for the normalization shall be printed.
  bool normalization_basis(){return PRINT_NORMALIXATION_BASIS;}
  //! @return True if every values in hte parsing shall be used as basis for the normalization process.
  bool everyrel_as_arrnorm_basis(){return USE_EVERYREL_AS_ARRNORM_BASIS;}
  //! @return The locaiton of where to store temporary- end resulting data from the filtering.
  char *get_file_binary_location(){return FILE_BINARY_LOCATION;}
  //! @return The name of the blast file to be parsed, defined by the user.
  char *get_file_input_name(){return FILE_INPUT_NAME;}
  //! @return The seperator to be used in the parsing of the blast file, defined by the user.
  char get_seperator(){return SEPERATOR;}
  //! @return The Number of cpu's set by the user.
  int get_cpu_tot(){return CPU_TOT;}
  //! @return The **norm object holding information about the total number- and sum of relations for each taxon-pair.
  list_norm_t *get_arrNorm(){return arrNorm;}
  /**
     @param <taxon_l> Sets this variable to the total number of taxa found in the collection.
     @return The collection of the taxa found with additional meta-data.
   **/
  taxa_t *get_listTaxa(int &taxon_l){taxon_l = taxon_length;return listTaxa;}
  //! Returns the highest similairty score found in the blast file parsed.
  float get_max_input_value(){return max_input_value;}
  //! Returns the total number of taxa found in the collection.
  int get_taxon_length(){return taxon_length;}
  //! Writes general info about the class:
  void print_class_info();
  //! Deallocates the hash. @remarks Ends the life of the hash labeled proteins.
  void free_hashProteins(){
    if(listParseData) listParseData->free_hashProteins();
  }

  //! Deallocates the data for the taxa objects.
  void free_listTaxa() {
    if(listTaxa) {
      for(int i =0; i< taxon_length; i++) listTaxa[i].free_mem();
      delete [] listTaxa, listTaxa = NULL;
    }
  }
  //! Deallocates the data for the list_file_parse objects.
  void free_list_file_parse() {
    if(listParseData) {
      listParseData->free_memory(true);
      delete listParseData;
      //      for(int i =0; i< taxon_length; i++) listTaxa[i].free_mem();
      //      delete [] listTaxa, listTaxa = NULL;
    }
  }
  //! De-allocates memory for the list_norm object.
  void free_arrNorm(){list_norm::close(arrNorm);}
  //! De-allocates memory for this object.
  void free_memory();
  //! De-allocates memory for the object given as input.
  static void close(blast_parsing *&obj, const bool delete_internals);
  //! Initiates the list for parsing input arguments fromthe terminal:
  static cmd_list* init_cmd_list(char *DEFAULT_OPTION_NAME, uint &DEFAULT_OPTION_NAME_COUNT, char *FILE_INPUT_NAME);
  //! Maps the internal variables to the input given from the terminal.
  static void build_cmd_list(cmd_list *cmd, int argc, char *argv[]);
  /**
     Executes the main operation for the parsing:
     @param <comm> The communicator to use if the macro variable "USE_MPI" is activated.
     @param <log> The log object to store measurements in.
     @param <bp>  The object containing the results of the blast parsing.
     @param <number_of_nodes> If MPI is not used, it would imply using only one node.
     @return An object containing an structured version of the blast file.
  **/
#ifdef USE_MPI
  void start_parsing(MPI_Comm comm, log_builder_t *log, bp_container_t &bp, int number_of_nodes);
#else
  /**
     Executes the main operation for the parsing:
     @param <log> The log object to store measurements in.
     @param <bp>  The object containing the results of the blast parsing.
     @param <number_of_nodes> If MPI is not used, it would imply using only one node.
     @return An object containing an structured version of the blast file.
  **/
  void start_parsing(log_builder_t *log, bp_container_t &bp, int number_of_nodes);
#endif

#ifdef USE_MPI
  /**
     Executes the main operation for the parsing:
     @param <comm> The communicator to use if the macro variable "USE_MPI" is activated.
     @param <log> The log object to store measurements in.
     @param <number_of_nodes> If MPI is not used, it would imply using only one node.
     @return An object containing an structured version of the blast file.
  **/
  list_file_parse_t* start_parsing(MPI_Comm comm, log_builder_t *log, int number_of_nodes);
#else
  /**
     Executes the main operation for the parsing:
     @param <log> The log object to store measurements in.
     @param <number_of_nodes> If MPI is not used, it would imply using only one node.
     @return An object containing an structured version of the blast file.
  **/
  list_file_parse_t* start_parsing(log_builder_t *log, int number_of_nodes);
#endif
  /**
     @brief the constructor.
     @param <cmd> The object to add variables accessible from the terminal.
  **/
  blast_parsing(cmd_list *cmd);  
  //! The test function for this class.
  static void assert_class(bool print_info);
};

/**
   @brief Handling a blast file, dumping them into retrievable containers.
   @ingroup parsing
   @author Ole Kristian Ekseth (oekseth)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
**/
typedef class blast_parsing blast_parsing_t;
#endif
