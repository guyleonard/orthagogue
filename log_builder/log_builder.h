#ifndef log_builder_h 
#define log_builder_h
/** 
   @file
   @brief Builds- and updates a log structure.
   @ingroup log_builder
   @author Ole Kristian Ekseth (oekseth)
   @date 17.12.2011 by oekseth (initial)
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
//#include "types.h"
#include "../blast_common/types.h"
#include "../blast_common/libs.h"
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
#include <typeinfo>       // operator typeid
#include <exception>      // std::exception

#include "enum_logid.h"
#include "cmd_list.h"
#include "../blast_parsing/tsettings_input.h"


/**
   @brief Identifies the different warning types
   @enum warnings
   @ingroup log_builder
   @author Ole Kristian Ekseth (oekseth)
**/
enum warnings {
  input_terminal,
  blastp_syntax,
  pointer,
  software_dependencies,
  memory_constraint,
  log_file,
  software_error // When a test does not pass
};

/**
   @brief Identifies the different warning types
   @ logid
   @ingroup log_builder
   @author Ole Kristian Ekseth (oekseth)
**/
typedef enum warnings warnings_t;


/**
   @class log_builder
   @brief Builds- and updates a log structure.
   @ingroup log_builder
   @remark The goal is to provide data for easy reconstruction- and evaluation
   of the process enabling scientific reability of time consumption given.
   @author Ole Kristian Ekseth (oekseth)
   @date 17.12.2011 by oekseth (initial)
**/
class log_builder {
 private:   // Variables added in order to ease the dump of data to a file
  // Using tables instead of variable-lsit, with enums as index-accesspoints:
  double list_time[logid_size];
  double list_pipes[logid_size][3];
  clock_t clock_time[logid_size];
  clock_t clock_pipes[logid_size][3];
  loint list_memory[logid_size];
  // Important variables for the purpose of describing the system
  // TODO: Consider removing the variables below and replace them with e.g. tsettings_input.
  char file_name[200];
  char seperator;
  lint FILE_BLOCK_READ_SIZE;
  int cpu_cnt;
  bool silent_run; // If set to true, does a silent run.
  bool verbose_as_it_goes; // If set to true, prints information as it is set.
  char **terminal_input;
  uint terminal_input_cnt;
  //! Initiates the above data containers to emptyness (or for the time the current).
  void init_containers();
  static void err_sys(char *err);
 public:

  /**
     @brief Perform comparisons of floats.
     @param <value1> The first value to use for comparison.
     @param <value2> The second value to use for comparison.
     @param <size_of_original_components> The value used calculating the offset-coefficient, required to avoid error due to rounding (when the values given as parameters to this function were generated).
     @param <line> The line this function was called from
     @param <file> The file this function was called from.
     @param <description> Short description of the test performed.
     @remarks 
     # Background: Extensive tests written August 08 2012 by oekseth due to problems concluded to originate in rounding of floats, bringing life to this function.
     # Include here for easy update of mulitple tests:
     @todo Consider moving this function to an own class handling such things.
  **/
  static void compare_floats(float value1, float value2, loint size_of_original_components, const int line, const char *file, const char *description); 
  /**
     @brief Catch errors in memory location.
     @return true if the error was of type std::bad_alloc.
  **/
  static bool catch_memory_exeception(const long int allocated_memory_size, const char *function, const char *file, const int line, const bool print_warning = false); 


  //! Test if condidtion is passed.
  static void test_memory_condition_and_if_not_abort(bool condition, const int line, const char *file, const char *function) {
    /*
      log_builder::test_memory_condition_and_if_not_abort(listTaxa!=NULL, __LINE__, __FILE__, __FUNCTION__);
    */
    if(!condition) {
      fprintf(stderr, "!!\tWas not able allocating data due to memory constraints. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", line, file, function);
      exit(2);
    }
  }

  
  //! Throws a standarised error, and abort.
  static void throw_error_and_abort(const char *message, const int line, const char *file, const char *function) {
    fprintf(stderr, "!!\t%s. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", message, line, file, function);
    exit(2);
  }

  //! Throws a standarised error for file reading, and abort.
  static void test_file_reading_and_if_not_abort(FILE *file_pointer, char *file_name, const int line, const char *file, const char *function) {
    if(!file_pointer) {
      fprintf(stderr, "!!\t Unable to open file %s. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", file_name, line, file, function);
      exit(2);
    }
  }

  //! Throws a specific error, given the "warnings_t" type, using the message as basis for specificity.
  static void throw_warning(warnings_t type, const int line, const char *file, const char *function,const char *message) {
    /*
      log_builder::throw_warning(software_dependencies, __LINE__, __FILE__, __FUNCTION__,"");
      char str[100]; memset(str, '\0', 100);
      log_builder::throw_warning(software_error, __LINE__, __FILE__, __FUNCTION__,str);
    */
    if(type == input_terminal) fprintf(stderr, "usage-warning!!\t%s. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", message, line, file, function);
    // TODO: Consider putting in the macro variables below!
#ifndef NDEBUG    
    else if(type == blastp_syntax) fprintf(stderr, "blastp-data-discarded\t%s. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", message, line, file, function);
#endif
    else if(type == pointer) fprintf(stderr, "pointer-error!!\t%s. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", message, line, file, function);
    else if(type == software_dependencies) { 
      fprintf(stderr, "software_dependency-error!!\t%s. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", message, line, file, function);
    } else if(type == memory_constraint) fprintf(stderr, "memory_constraint-error!!\t%s. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", message, line, file, function);
    else if(type == log_file) fprintf(stderr, "log-file-writing-error!!\tFound it difficult opening '%s'. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", message, line, file, function);
    else if(type == software_error) fprintf(stderr, "software-test-did-not-pass!!\tFound it difficult performing '%s'. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", message, line, file, function);
  }

  //! @return the file size:
  static loint get_file_size(char *file, const int line_caller, const char *file_caller);

  //! Takes a pause in the program; useful for analysing- and debugging.
  static void pause_program(const int line, const char *file, const char *function) {
    /*
    log_builder::pause_program(__LINE__, __FILE__, __FUNCTION__); 
    */
    char x = ' ';
    printf("Press any key to continue...at line %d in file %s for method %s\n", line, file, function);    
    x = getchar();
    printf("The key you entered is: %c \n", x);
  }
  //! Sets the terminal input line: If set, the data is included in the textual logs.
  void set_program_inits(uint cnt, char **array) {terminal_input=array, terminal_input_cnt=cnt;}
  //! Sets the status of the verbouse-parameter.
  void set_verbose_as_it_goes(bool verbouse) {verbose_as_it_goes = verbouse;}
  //! Starts the log measurement of a given pipe
  void start_measurement(logid_t type_id, uint pipe_id);
  //! Appends information for the  log measurement of a given pipe.
  void append_measurement(logid_t type_id,uint pip_id,const clock_t start,const clock_t end);
  //! Ends the log measurement of a given pipe
  void end_measurement(logid_t type_id, uint pipe_id);  
  //! Starts the log measurement for a given general id
  void start_measurement(logid_t type_id);
  //! Ends the log measurement for a given general id
  void end_measurement(logid_t type_id);
  /**
     @brief Creates the directory to store the log-data in.
     @param <size> The size of the log-folder-name
     @return true if the path is ok;
  **/     
  static bool create_log_folder(loint &size) {
      bool path_ok = true;
#ifdef LOG_FOLDER_NAME
    size+=sizeof(FUNC_LOG_FOLDER_NAME);
    struct stat st;
    if(stat(FUNC_LOG_FOLDER_NAME,&st) != 0) { // The path does not already exists
      if(-1 == mkdir(FUNC_LOG_FOLDER_NAME, S_IRWXU)) {path_ok = false;}
    }
    return path_ok;
#else
    return true;
#endif
  }
  
  //! @return the file pointer to a log-file-location
  static FILE *get_log_file_pointer(char *file_description, uint cpu_cnt, char *caller_file, int caller_line) {
    char sep = '_';
    // The log file, and contents of it:
    if(sizeof(FUNC_LOG_FOLDER_NAME) > 0) sep = '/';
    char str[200]; memset(str, '\0', 200);
#ifdef USE_MPI
    //! Some mpi variables:
    int myrank = 0;  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    // Gets specific data about the reading.
    int number_of_nodes = 1; MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes); 
    if(number_of_nodes > 1) {
      sprintf(str, "%s%c%s.log.%u.%d", FUNC_LOG_FOLDER_NAME, sep, "list_file_parse", cpu_cnt, myrank);
    } else sprintf(str, "%s%c%s.log.%u", FUNC_LOG_FOLDER_NAME, sep, "list_file_parse", cpu_cnt);
#else
    sprintf(str, "%s%c%s.log.%u", FUNC_LOG_FOLDER_NAME, sep, "list_file_parse", cpu_cnt);
#endif
    FILE *f_log = fopen(str, "w");
    if(!f_log) {fprintf(stderr, "!!\tUnable to open file %s at %s:%d\n", str, caller_file, caller_line);}
    return f_log;
  }
 private:
  bool is_alpha_num(char *value);   //! Returns true if argument has alphanumeric values,
  FILE *open_file(char *str, bool is_text);
  FILE **build_log_path();
  char *get_date_string();
  /**
     @brief PRoduces a row for use with gnuplot
     @param <file> The stream to write data to.
     @param <add_data> If set to false adds the column headers.
     @remarks Seperated out as own method iot ensure consistency.
  **/
  void produce_row(FILE *file, bool add_data);
  //! Outputs the results of the time measurements.
  void print_general_status(FILE *file);
  //! Below code to hide variables fom apperating as 'unused'
  void print_log_names();
 public:
  /**
     @brief Updates with the settings given by the user.
     @param <sinp> Defines the settings used during the parsing.
  **/
  void update_settings(tsettings_input_t sinp); 
  /**
     @return a file pointer to the log file.
  **/
  static FILE * get_log_file_pointer(uint n_threads, char *description, char *code_file, int code_file_line, const bool create_new_file);

  /**
     @return a file pointer to the log file.
  **/
  static FILE * get_log_file_pointer(uint n_threads, char *description, char *code_file, int code_file_line, const bool create_new_file, int myrank);

  //! Releases the memory.
  void free_mem();
  //! Initializes an object of class log_builder
  static log_builder *init() {return new log_builder();}
  //! Releases the memory of an object given.
  static void close(log_builder *&obj);
  //! Produces the log file given the measurements taken.
  void write_log_file();
  //! The constructor.
  log_builder();
  //! The constructor.
  log_builder(char *f_name, char _sep, lint b_r_size, lint di_size, int _cpu_cnt);    
  /**
     @brief The constructor.
     @param <sinp> Defines the settings used during the parsing.
   **/
  log_builder(tsettings_input_t sinp);    
  /**
     @brief The constructor.
     @param <cmd> The object to add variables accessible from the terminal.
   **/
  log_builder(cmd_list *cmd);    
  //! The main test function for this class  
  static void assert_class(const bool print_info);
};

/**
   @brief Builds- and updates a log structure.
   @ingroup log_builder
   of the process enabling scientific reability of time consumption given.
   @author Ole Kristian Ekseth (oekseth)
   @date 17.12.2011 by oekseth (initial)
**/
typedef class log_builder log_builder_t;
#endif
