/**
   @file
   @brief Code producing a filtered output of the input.
   @ingroup filtering
**/
#ifndef blast_filtering_h
#define blast_filtering_h
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
#include "cmd_list.h"
#include "../configure.h"
#include "log_builder.h"
#include "taxa.h"
#include "id_simil_list.h"
#include "list_file_parse.h"
#include "enum_mcl.h"
#include "../blast_parsing/bp_container.h"
#include "cmd_list.h"
// The pipes:
#include "pipe_bucket.h"
#include "pipe_ortholog_inparalog.h"
#include "pipe_merge.h"
#include "pipe_norm.h"
#include "pipe_write.h"
#include "pipe_struct.h"
#include "../blast_parsing/tsettings_input.h"
#include "meta_blast.h"
using namespace std;
/**
   @class blast_filtering
   @brief Produces a filtered output of the input.
   @ingroup filtering
   @remarks The purpose of this wrapper module is handling a the filtering process given an input of a taxa_t- and list_file_parse containers, producing a filtered output.
   @author Ole Kristian Ekseth (oekseth)
   @date 21.12.2010 by Ole Kristian Ekseth (init)
   @date 18.08.2011 by Ole Kristian Ekseth (Cleaning.)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
   @date 31.12.2011 by oekseth (clean-up)
**/
class blast_filtering {
 private:
  // Variables below are either set by the constructor or by the cmd-line param given as input
  bool OUTPUT_PIPE_MCI_ALL;
  bool DEBUG_NORM; 
  bool PRINT_NORMALIXATION_BASIS; 
  bool USE_EVERYREL_AS_ARRNORM_BASIS;
  char *FILE_BINARY_LOCATION; 
  //  char SEPERATOR;
  int CPU_TOT;
  // Below are actual data that must be filled on forehand iot do some usefull tasks.
  list_norm_t *arrNorm;
  taxa_t *listTaxa;
  int taxon_length;
  list_file_parse_t *listParseData;
  float max_input_value;

  // Below constants used without modification, but with an option for further scope.
  static const bool MODE_PAIRWISE_OUTPUT_ABC = true; // for the abc file: if set, the data out pairwise in stead of as in a row
  static const bool MODE_PAIRWISE_OUTPUT_MCL = false; // for the mcl file: if set, the data out pairwise in stead of as in a row

  // Below are config-settings for this scope of the processing.
  bool DEBUG_PRINT_DISCARDED_PAIRS;
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT; // If unset, uses the raw similarity score for the abc-files.
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT; // If unset, uses the raw similarity score for the mci-files.
  bool MODE_PAIRWISE_OUTPUT, MODE_INTEGER_OUTPUT; // if set, the data is written in a pairwise- in stead ofwise format.
  bool PRINT_IN_ABC_FORMAT, PRINT_IN_MCL_FORMAT; // If unset, the *.<file> (i.e. the file with <type> as protein names) discardes
  bool PRINT_OVERLAP_VALUES_ABOVE; // If set, dumps the discarded value based on the overlap criteria.
  bool RESTRICTED_DEFENITION; // If set the co-orthologs (relations between inparalogs of orthologs) are not in the results.
  bool SORT_ABC_DATA; // If set, sorts the abc files before printing them. Implemented ineffecitvely.
  short int AMINO_LIMIT; // The threshold value for the overlap.
  float MIN_SIMILARITY_LIMIT; // How close the proteins must be in order not to be ignored .
  static const bool use_improved_overlap_algo = true;
  // Below are data containers for this scope of the processing.
  list_file_struct_t *listStructData;
  id_simil_list listOrtho; 
  //! Holds properties about the blast file:
  meta_blast_t blastInfo;
  //  class meta_blast blastInfo;

  // Functions for filtering of inparalogs- and orthologs:
  void exec_ortho_operation(int n_thread, uint taxon_start, uint taxon_end, log_builder_t *log);
  void exec_operation(int n_thread, int taxon_start, int taxon_end, bool is_inparalog, log_builder_t *log);
  void exec_ortho(int n_thread, log_builder_t *log);
  void exec_inpa(int n_thread, log_builder_t *log);
  void exec_filtering(int n_thread, log_builder_t *log, const bool is_inpa);
  void exec_co_orth(int n_thread, log_builder_t *log, stack_rel *&stackRel);
  void write_result_file(int n_thread, log_builder_t *log, stack_rel *&stackRel, pipe_write &write, float **&arr_avgNorm);
 public:
  //! Prints information describing this class.
  void print_class_info();
/** 
    @brief Updates the command line interface with values to be set by the user:
    @param <cmd> The input list to be made available from the terminal console.
    @parm <first_pass> In order to get the correct order of the fields, i.e. "OUTPUT" before "OPERATIONAL".
 **/
  void init_values(cmd_list *cmd, const bool first_pass=true);
  //! Initiates the list for parsing input arguments fromthe terminal:
  static cmd_list* init_cmd_list(char *DEFAULT_OPTION_NAME, uint &DEFAULT_OPTION_NAME_COUNT, char *FILE_INPUT_NAME);
  //! Maps the internal variables to the input given from the terminal.
  static void build_cmd_list(cmd_list *cmd, int argc, char *argv[]);
  /**
     @brief Sets the user settings from the blast parsing process.
     @param <obj> The object of class tsettings_input containing the data to use.
     @remarks The filtering process will depend upon some of these settings, therby relevant.
   **/
  void set_values(tsettings_input obj);
  //! Initializes values after those set in the blast parsing (dound in library blast_parsing):
   void set_values(bool _DEBUG_NORM, bool _PRINT_NORMALIXATION_BASIS, bool _USE_EVERYREL_AS_ARRNORM_BASIS, char *_FILE_BINARY_LOCATION/*,  char _SEPERATOR*/, int _CPU_TOT, bool USE_LAST_BLAST_CLOMUN_AS_DISTANCE);

  /**
     @brief Executes the main operation for the filtering:
     @param <log> The log object to store measurements in.
     @param <bp>  The object containing the basis whom the filtering will work on.
  **/
  void start_filtering(log_builder_t *log, bp_container_t &bp);
  //! De-allocates memory for this object.
  void free_memory();
  //! De-allocates memory for the object given as input.
  static void close(blast_filtering *&obj);
  //! The constructor.
  blast_filtering(cmd_list *cmd);

  //! The assert method:
  static void assert_class(bool print_info);
};
/**
   @brief Produces a filtered output of the input.
   @ingroup filtering
   @author Ole Kristian Ekseth (oekseth)
   @date 31.12.2011 by oekseth (clean-up)
**/
typedef class blast_filtering blast_filtering_t;
#endif
