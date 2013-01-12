#include "blast_parsing_main.h"
/**
   @brief This file is regarded as the launcher of the software, implying using the librariy  'blast_parsing'. 
   @remarks In order to understand the above library, this file might be of use the the developer.
   @author Ole Kristian Ekseth (oekseth)
   @date 27.12.2011 by oekseth (cleanup)
**/
int main(int argc, char *argv[]) {
  blast_parsing::assert_class(true); 
  if(true) {
    char *DEFAULT_OPTION_NAME = "OPTIONS"; // Used when no option is specified
    uint DEFAULT_OPTION_NAME_COUNT = 0; // counts the number of times the 'options' is used
    char *FILE_INPUT_NAME ="temp.raritet.txt";
    // Constructs the object taking input from the terminal:
    cmd_list *cmd = blast_parsing::init_cmd_list(DEFAULT_OPTION_NAME, DEFAULT_OPTION_NAME_COUNT, FILE_INPUT_NAME);
    // Initializes the library we are interested in:
    blast_parsing * blast = new blast_parsing(cmd);
    // Maps internal variables with the terminal input:
    blast_parsing::build_cmd_list(cmd, argc, argv);
    blast->print_class_info();
    if(FILE_INPUT_NAME) {
      // Gets the user configuration parameters:
      const char SEPERATOR = blast->get_seperator();
      const int CPU_TOT = blast->get_cpu_tot();

      // Starts the logg tool gaining spesific measurements of time usage:
      FILE_INPUT_NAME = blast->get_file_input_name();
      log_builder *log = new log_builder(FILE_INPUT_NAME, SEPERATOR, BLOCK_FILE_READ_SIZE, 1024*1024, CPU_TOT);
      log->start_measurement(complete_runningtime);

      // Executes the parsing of blast file into efficient structures:
      bp_container_t bp;
#ifdef USE_MPI
      blast->start_parsing(MPI_COMM_WORLD, log, bp, 1);  
#else
      blast->start_parsing(log, bp, 1);  
#endif

      // Completes the logging operation:
      log->end_measurement(complete_runningtime);
      log->write_log_file();
    
      // Frees allocated memory:
      log_builder::close(log);
    }
    // Frees allocated memory:
    //    blast_parsing::close(blast);
    cmd_list::close(cmd);
  }

}
