#include "log_builder.h"

static const char *pipe_name[3] = {
  "get", "extract", "merge"
};
static const char *logid_name[logid_size] = {
  "read_first", "collect_read", "read_second", 
  "collect_final", "collect_lfp", "collect_overlap",
  "ortholog_init", "orthologs", "inparalogs", "co_orthologs",
  "write_result",
  "complete",
  "blast_parsing", "blast_filtering"
};


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
void log_builder::compare_floats(float value1, float value2, loint size_of_original_components, const int line, const char *file, const char *description) {
  //! Usage:
  /*
    log_builder::compare_floats(, , , __LINE__, __FILE__, __FUNCTION__);
  */
#ifndef NDEBUG
  int myrank=0; 
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    
#endif
  // TODO: The below 'coefficient' must be creaful evaluated, in order to get high precision, wihout the cost of thorwing an error when the error is due to the im-precision of the floating-point-numbers:
  //! Downscales the floating points to avoid error when the difference is due to rounding:
  float coefficient = 10/(10+size_of_original_components); // TODO: Must validate this coefficient
  if(coefficient > 0.0001) coefficient = 0.001; // TODO: A local adjustment.
  //! Compares, using the downscaling:
  const int adjust_1 = (int)(coefficient*value1);
  const int adjust_2 = (int)(coefficient*value2);
  if(adjust_1 != adjust_2) {
    fprintf(stderr, "!!\t[myrank=%d]\tFailed for %s, as %f != %f (wth coefficient=%f => (%d =?= %d)), at %s:%d\n",
	    myrank, description, value1, value2, coefficient, adjust_1, adjust_2, file, line);
    assert((coefficient*value1) == (coefficient*value2));
  }
#endif

}

//! Initiates the above data containers to emptyness (or for the time the current).
void log_builder::init_containers() {
  for(uint i = 0; i < logid_size; i++) {
    list_memory[i] = 0, list_time[i] = 0; clock_time[i] = clock_t();
    for(uint o = 0; o < 3; o++) list_pipes[i][o] = 0, clock_pipes[i][o] = clock_t();;      
  }
}

void log_builder::err_sys(char *err) {
  fprintf(stderr, "!!\n!!\tAborts due to %s at line %d in class log_builder.\n", err, __LINE__);
  exit(2);
}

void log_builder::start_measurement(logid_t type_id, uint pipe_id) {
  assert(pipe_id < 3);
  struct tms tmsstart;
  if((clock_pipes[type_id][pipe_id] = times(&tmsstart)) == -1) // starting values
    err_sys("times error");
}

void log_builder::append_measurement(logid_t type_id, uint pipe_id, const clock_t start, const clock_t end) {
  if(false) print_log_names(); // In order to hide the "not used-warning
  long clktck = 0;
  if((clktck =  sysconf(_SC_CLK_TCK)) < 0) 	err_sys("sysconf error");  
  clock_t t_real = end - start;
  list_pipes[type_id][pipe_id] += t_real/(double)clktck;
  /*if((type_id==read_first) && (pipe_id = 0) &&  list_pipes[type_id][pipe_id] != 0) printf("in log_builder list_pipes[%s][%s] = %f\n",
    logid_name[type_id], pipe_name[pipe_id], list_pipes[type_id][pipe_id]);*/
}

void log_builder::end_measurement(logid_t type_id, uint pipe_id) {
  assert(pipe_id < 3);
  struct tms tmsend;
  clock_t end; long clktck = 0;
  if(list_pipes[type_id][pipe_id] != 0) { // Then increments:

  } else {
    if((end = times(&tmsend)) == -1) err_sys("times error"); // starting values
    if((clktck =  sysconf(_SC_CLK_TCK)) < 0) 	err_sys("sysconf error");
    clock_t t_real = end - clock_pipes[type_id][pipe_id];
    list_pipes[type_id][pipe_id] = t_real/(double)clktck;
  }
}

void log_builder::start_measurement(logid_t type_id) {
  if(false) print_log_names(); // In order to hide the "not used-warning
  list_memory[type_id] = (loint)sbrk(0);
  struct tms tmsstart;
  if((clock_time[type_id] = times(&tmsstart)) == -1) // starting values
    err_sys("times error");
}

void log_builder::end_measurement(logid_t type_id) {
  assert((uint)	 type_id < (uint)logid_size);
  struct tms tmsend; clock_t end; long clktck = 0;
  if((end = times(&tmsend)) == -1) err_sys("times error"); // starting values
  if((clktck =  sysconf(_SC_CLK_TCK)) < 0) 	err_sys("sysconf error");
  clock_t t_real = end - clock_time[type_id];
  list_time[type_id] = t_real/(double)clktck;
  // Sets the memory usage:
  const loint mem_pos_end = (loint)sbrk(0);
  if(mem_pos_end > list_memory[type_id]) 
    list_memory[type_id] = mem_pos_end - list_memory[type_id];
  else list_memory[type_id] = 0; // For small values, this may happen.

  if(verbose_as_it_goes) {
    fprintf(stdout, "-%c%-25s%5s", ' ', logid_name[type_id], " used");
    fprintf(stdout, "%10.2f seconds: ", list_time[type_id]);
    fprintf(stdout, "\n");
  }
}

//! Returns true if argument has alphanumeric values,
// assuming that the uinput has a '\0' at its end.
bool log_builder::is_alpha_num(char *value) {
  if(value != NULL) {
    bool is_alpha_num = true;
    char *start = value;
    char *end = strrchr(value, '\0');    
    if(end != NULL)
      while(start != end) {
	if(!(isalnum(*start) || ispunct(*start))) return false;
	start++;
      }
    //      if(!isalnum(*start++)) return false; 
    else {
      is_alpha_num=false; // do not have an alpha numeric value at its end
      //      printf("!!\tdid not find ending char\n");
    }
    if(is_alpha_num) return true; // at this point, the tests are succesful
    else {
      fprintf(stderr, "!!\tError: String '%s' were not in a pure alpha numeric mode: Email the devloper (found in the documentation) sending a copy of this warning to the devloper.\n", value);
      return false;
    }
  } else return false;
}

FILE *log_builder::open_file(char *str, bool is_text) {
  if(str != NULL) {
    FILE *f;
    if(is_text) f = fopen(str, "w");
    else f = fopen(str, "wb");
    if(f!= NULL)  return f;
    else {fprintf(stderr, "!!\tUnable to open file in class log_builder of type %d and name '%s'\n", is_text, str); assert(false);}
  }
  return NULL;
}

FILE **log_builder::build_log_path() {
#ifdef BUILD_LOG_FILES
  loint size = 100;
  bool path_ok = create_log_folder(size);
  if(path_ok) {
    { // 	 grant the owner read, write and execute authority
      const uint str_size = 3;
      char **str = new char*[str_size];
      for(uint i =0;i<str_size;i++) {str[i] = new char[size]; memset(str[i], '\0', size);}
      //  const static char *res_file = "results_turbo.txt";
      char res_file[100]; memset(res_file, '\0', 100); 
      int myrank=0; bool mpi = false;
#ifdef USE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      mpi = true;
#endif
      if(mpi) sprintf(res_file, "results_orthAgogue.txt.rank_%d", myrank);
      else    sprintf(res_file, "results_orthAgogue.txt");
#ifdef LOG_FOLDER_NAME
      char sep = '_';
      if(sizeof(FUNC_LOG_FOLDER_NAME) > 0) sep = '/';
      //! The textural result files:
      if(mpi) sprintf(str[0], "%s%corthAgogue.log.%d.rank_%d", FUNC_LOG_FOLDER_NAME, sep, cpu_cnt, myrank);
      else    sprintf(str[0], "%s%corthAgogue.log.%d", FUNC_LOG_FOLDER_NAME, sep, cpu_cnt);
      //! The binary result files:
      if(mpi) sprintf(str[1], "%s%corthAgogue.bin.%d.rank_%d", FUNC_LOG_FOLDER_NAME, sep, cpu_cnt, myrank);
      else    sprintf(str[1], "%s%corthAgogue.bin.%d", FUNC_LOG_FOLDER_NAME, sep, cpu_cnt);
      //! The concatenated result for for all nodes:
      sprintf(str[2], "%s%c%s", FUNC_LOG_FOLDER_NAME, sep, res_file);
#else
      if(mpi) {
	sprintf(str[0], "orthAgogue.log.%d.rank_%d", cpu_cnt);
	sprintf(str[1], "orthAgogue.bin.%d.rank_%d", cpu_cnt);
      } else {
	sprintf(str[0], "orthAgogue.log.%d", cpu_cnt);
	sprintf(str[1], "orthAgogue.bin.%d", cpu_cnt);
      }
      sprintf(str[2], "%s", res_file);
#endif
      FILE ** f_list = new FILE*[3];
      f_list[0] = open_file(str[0], true);
      f_list[1] = open_file(str[1], false);
 
      bool result_file_is_empty = false; 
      struct stat stat_new;       // Adds a header if the file is not in the system:
      if(stat(str[2], &stat_new) != 0) result_file_is_empty = true;
      f_list[2] = fopen(str[2], "a+"); // Creates a new file it it doesn't exist, and appends to the end of it.
      if(f_list[2]) {
	if(result_file_is_empty) { // The file is not set: makes a header 
	  produce_row(f_list[2], false);
	}
      }
      for(uint i = 0; i <str_size; i++) {delete [] str[i];} delete [] str; str = NULL;
      return f_list;
    } /*else {
	fprintf(stderr, "!!\tIn class log_builder, unable to build a log directory at path: %s\n", log_folder_name); fflush(stderr);
	assert(false);
	}*/
  } else {
    fprintf(stderr, "!!\tIn class log_builder, do not regard the directory given as a set of alphanumerical chars, therefore disregards the path given: %s\n", FUNC_LOG_FOLDER_NAME); fflush(stderr);
    assert(false);
  }
  return NULL;
#endif
  return NULL;
}


/**
   @brief PRoduces a row for use with gnuplot
   @param <file> The stream to write data to.
   @param <add_data> If set to false adds the column headers.
   @remarks Seperated out as own method iot ensure consistency.
**/
void log_builder::produce_row(FILE *file, bool add_data) {
  if(add_data) {
    // Adds a row of information to the log file:
    fprintf(file, "%d\t", cpu_cnt);
    for(uint i =0; i <logid_size; i++) {
      fprintf(file, "%10.2f\t", list_time[i]);
    }    
    fprintf(file, "#\n");
  } else {
    fprintf(file, "#\t");
    for(uint i =0; i <logid_size; i++) {
      fprintf(file, "%s\t", logid_name[i]);
    }
    fprintf(file, "#\n");
  }
}
char *log_builder::get_date_string() {
  char *buffer = new char[100]; memset(buffer, '\0', 100);
  time_t curtime = time(NULL); // Get the current time.  
  struct tm *loctime = localtime(&curtime); // Convert it to local time representation.
  //    fputs (asctime (loctime), stdout);  // Print out the date and time in the standard format.  
  strftime(buffer, 100, "%A, %B %d %I:%M %p", loctime);
  //strftime (buffer, SIZE, "The time is %I:%M %p.\n", loctime);
  return buffer;
}

//! Outputs the results of the time measurements.
void log_builder::print_general_status(FILE *file) {
  char *str = get_date_string();
  fprintf(file, "The settings for the orthAgogue run at %s were: ",str);
  delete [] str; str = NULL;
  fprintf(file, "file_name(%s), FILE_BLOCK_READ_SIZE(%dKB), #cpus(%d).", file_name,  (int)FILE_BLOCK_READ_SIZE/(1024), (int)cpu_cnt);
  if(terminal_input) {
    fprintf(file, "\n-\t The software was started using the following terminal line: \"");
    for(uint i = 0; i < terminal_input_cnt; i++)
      fprintf(file, "%s ", terminal_input[i]);
    fprintf(file, "\".\n");
  }
  fprintf(file, "- Given these configurations, the measurements became:\n");
  // The time-settings:
  //  printf("...logid_size=%d\n", (int)logid_size);
  for(uint i = 0; i < logid_size; i++) {
    fprintf(file, "-%c%-25s%5s", ' ', logid_name[i], " used");
    //      if(list_time[i]) 
    fprintf(file, "%10.2f seconds: ", list_time[i]);
    for(uint o = 0; o < 3; o++) {
      if(list_pipes[i][o]) 	fprintf(file, "%10s:%-4.4f ", pipe_name[o], list_pipes[i][o]);
      else                      	fprintf(file, "%9s%-9c ", " ", '_');  
    }
    if(list_memory[i] > 1024)
      fprintf(file, " using %lluMB of memory.", list_memory[i]/(1024*1024));      
    fprintf(file, "\n");
  }
}

//! Below code to hide variables fom apperating as 'unused'
void log_builder::print_log_names() {
  for(uint i = 0; i < logid_size; i++) {
    printf("logid_name[%d] = %s:", i, logid_name[i]);
    for(uint o = 0; o < 3; o++) {
      printf(" pipe_name[%d] = %s ", o, pipe_name[o]);
    }
    printf("\n");
  }
}
void log_builder::write_log_file() {
#ifdef BUILD_LOG_FILES
  FILE **f_list = build_log_path(); 
  if(f_list[1]) {
    fwrite(this, sizeof(this), 1, f_list[1]); // Writes the binary file
    fclose(f_list[1]);
  }
  if(f_list[0]) {// Then does the tiresome process of writing the text-file:
    if(!silent_run) print_general_status(stdout);
    print_general_status(f_list[0]);
    fclose(f_list[0]);
  } else fprintf(stderr, "!!\tunable to open the file at line %d in class log_builder.\n", __LINE__);
  if(f_list[2]) {
    produce_row(f_list[2], true);
    fclose(f_list[2]);
  }
  if(f_list)  {delete [] f_list, f_list = NULL;}
#endif
}

log_builder::log_builder() :
  seperator('_'), FILE_BLOCK_READ_SIZE(0), 
  cpu_cnt(0), silent_run(true), verbose_as_it_goes(false), terminal_input(NULL), terminal_input_cnt(0)
{
  memset(file_name, '\0', 200); 
  init_containers();
}


log_builder::log_builder(char *f_name, char _sep, lint b_r_size, lint di_size, int _cpu_cnt) :
  seperator(_sep), FILE_BLOCK_READ_SIZE(b_r_size),
  cpu_cnt(_cpu_cnt), silent_run(true), verbose_as_it_goes(false),  terminal_input(NULL), terminal_input_cnt(0)
{
  if(false) printf("di_size=%d", (int)di_size);
  memset(file_name, '\0', 200);
  if(f_name)  memcpy(file_name, f_name, sizeof(f_name)+1);
  init_containers();
}

/**
   @brief The constructor.
   @param <sinp> Defines the settings used during the parsing.
**/
log_builder::log_builder(tsettings_input_t sinp) : 
  seperator(sinp.SEPERATOR), FILE_BLOCK_READ_SIZE(BLOCK_FILE_READ_SIZE),
  cpu_cnt(sinp.CPU_TOT), silent_run(true), verbose_as_it_goes(false), terminal_input(NULL), terminal_input_cnt(0)
{
  memset(file_name, '\0', 200);
  if(sinp.FILE_NAME)  memcpy(file_name, sinp.FILE_NAME, sizeof(sinp.FILE_NAME)+1);
  init_containers();
}

/**
   @brief The constructor.
   @param <cmd> The object to add variables accessible from the terminal.
**/
log_builder::log_builder(cmd_list *cmd) : 
  seperator('_'), FILE_BLOCK_READ_SIZE(0), 
  cpu_cnt(0), silent_run(true), verbose_as_it_goes(false), terminal_input(NULL), terminal_input_cnt(0)
{
  memset(file_name, '\0', 200); 
  init_containers();
#ifndef NDEBUG
  class cmd_argument cl2; 
  cl2 = cmd_argument("Show the log  after program termination", "L", "show_log", BOOLEAN, &silent_run, "DEBUG");
  cmd->add_cmd_argument(cl2);  // USE_EVERYREL_AS_ARRNORM_BASIS
#endif
}

/**
   @brief Updates with the settings given by the user.
   @param <sinp> Defines the settings used during the parsing.
**/
void log_builder::update_settings(tsettings_input_t sinp) {
  seperator = sinp.SEPERATOR,  cpu_cnt = sinp.CPU_TOT;
  if(sinp.FILE_NAME)  memcpy(file_name, sinp.FILE_NAME, sizeof(sinp.FILE_NAME)+1);
}

/**
   @return a file pointer to the log file.
   @remarks If in NDEBUG mode return NULL, discarding this option.
**/
FILE * log_builder::get_log_file_pointer(uint n_threads, char *description, char *code_file, int code_file_line, const bool create_new_file) {
#ifndef NDEBUG
  return get_log_file_pointer(n_threads, description, code_file, code_file_line, create_new_file, INT_MAX);
#else
  return NULL;
#endif
}

//! @return the file size:
loint log_builder::get_file_size(char *file, const int line_caller, const char *file_caller) {
  if(file) {
    struct stat sb;
    if (stat(file, &sb) == -1) {
      char str[100]; sprintf(str, "Unable to open file %s, reqeusted from object found at file %s, line %d", file, file_caller, line_caller);
      throw_error_and_abort(str, __LINE__, __FILE__, __FUNCTION__);
    }
    const loint size = (loint)sb.st_size;
    return size;
  } else {
    char str[100]; sprintf(str, "File %s not defined, reqeusted from object found at file %s, line %d", file, file_caller, line_caller);
    throw_error_and_abort(str, __LINE__, __FILE__, __FUNCTION__);
    return 0;
  }
}

/**
   @return a file pointer to the log file.
**/
FILE * log_builder::get_log_file_pointer(uint n_threads, char *description, char *code_file, int code_file_line, const bool create_new_file, int myrank) {
  assert(description);
  char sep = '_';
  // The log file, and contents of it:
  bool is_folder = false;
  char log_path[100];memset(log_path, '\0', 100);
  //sprintf(log_path, "");
  char str[200]; memset(str, '\0', 200);
#ifdef LOG_FOLDER_NAME
#ifdef FUNC_LOG_FOLDER_NAME
  if(sizeof(FUNC_LOG_FOLDER_NAME) > 0) {
    sep = '/';
    is_folder = true;
  }
  sprintf(log_path, "%s", FUNC_LOG_FOLDER_NAME);
  if(is_folder) {
    struct stat stat_folder;       // Adds a header if the file is not in the system:
    if(stat(log_path, &stat_folder) == 0) {
      if(-1 == mkdir(log_path, S_IRWXU)) {
	// TODO: Consider doing something here:
	// -- Problem is that this error may be cast even though the directory has been created. (11.04.2012 by oekseth).
      }
    }
  } 
#endif
#endif
  if(INT_MAX == myrank) sprintf(str, "%s%c%s.log.%u", FUNC_LOG_FOLDER_NAME, sep, description, n_threads);
  else   sprintf(str, "%s%c%s.log.%u.%d", FUNC_LOG_FOLDER_NAME, sep, description, n_threads, myrank);
  FILE *f_log = NULL;
  if(create_new_file) f_log = fopen(str, "w");
  else f_log = fopen(str, "a");
  if(f_log) return f_log;
  else {
    log_builder::throw_warning(log_file, __LINE__, __FILE__, __FUNCTION__, str);
    return NULL;
  }
}
//! Frees memory
void log_builder::free_mem() {
  ;  //  if(log_folder_name) {delete [] log_folder_name, log_folder_name = NULL;}
}
void log_builder::close(log_builder *&obj) {
  if(obj) {obj->free_mem(), delete obj, obj = NULL;}
}

void log_builder::assert_class(const bool print_info) {
  const static char *class_name = "log_builder";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code      
  // TODO: Consider inserting something.
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}
