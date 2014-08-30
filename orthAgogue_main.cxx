#include "orthAgogue_main.h"
  
/* 
   \brief The Ececutable of the software
   @Name:  This file is regarded as the launcher of the software, implying using the two libraries 'blast_parsing' and 'blast_filtering'. In order to understand the aboe libraries, this file might be of use the the developer.
   @Author: Ole Kristian Ekseth (oekseth)
   @Last_Major_Update: 27.12.2011 by oekseth (cleanup)
**/
int main_operation_turbo(int array_cnt, char **array) { 
  int number_of_nodes = 1; // If MPI is not used, it would imply using only one node.
#ifdef USE_MPI
  MPI_Init(&array_cnt, &array);
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_nodes); 
#endif
  // Inputs the parameters from the command line termianl execution of this software:
  char *DEFAULT_OPTION_NAME = "OPTIONS"; // Used when no option is specified
  uint DEFAULT_OPTION_NAME_COUNT = 0; // counts the number of times the 'options' is used
  char *FILE_INPUT_NAME ="temp.raritet.txt";

  // Constructs the object taking input from the terminal:
  cmd_list *cmd = blast_parsing::init_cmd_list(DEFAULT_OPTION_NAME, DEFAULT_OPTION_NAME_COUNT, FILE_INPUT_NAME);

  // Initializes the libraries we are interested in:
  blast_parsing *blast = new blast_parsing(cmd);
  blast_filtering *filter = new blast_filtering(cmd);
  blast->init_values(cmd, false);
  filter->init_values(cmd, false);
  blast->init_values(cmd, 2);
  log_builder *log = new log_builder(cmd);

  // Maps internal variables with the terminal input:
  blast_parsing::build_cmd_list(cmd, array_cnt, array);

  // Gets the user configuration parameters:
  tsettings_input_t settings_input;
  blast->get_input_settings(settings_input);

  // Updates with the settings given by the user.
  log->update_settings(settings_input);
  // Starts the logg tool gaining spesific measurements of time usage:
  log->start_measurement(complete_runningtime);

  // Executes the parsing of blast file into efficient structures:
  bp_container_t bp;
#ifdef USE_MPI
  blast->start_parsing(MPI_COMM_WORLD, log, bp, number_of_nodes);  
#else
  blast->start_parsing(log, bp, number_of_nodes);  
#endif
  // Gets the updated containers- and supporting values:
  filter->set_values(settings_input); 
  if(DEBUG) {
    printf("in orthAgogue_main at line %d the blast_memory consists of:\n", __LINE__);
    assert(bp.listParseData);
    bp.listParseData->print_info_size_matrix();
  }

  cmd_list::close(cmd);
  blast_parsing::close(blast, false);
  
  if(!bp.listTaxa) 
    {
#ifdef USE_MPI
    /**!
       If no data were generated, we assume the input was too small for the number of nodes given:
    **/

    char string[100]; memset(string, '\0', 100);
    sprintf(string, "Node designated id(%d) is not set doing the prosessing, probably due to the disk buffer size set", myrank);
    log_builder::throw_warning(input_terminal, __LINE__, __FILE__, __FUNCTION__,string);
    MPI_Finalize();
    exit(1);
#else
    log_builder::throw_warning(input_terminal, __LINE__, __FILE__, __FUNCTION__, "No data generated");    
#endif
  } else {
    // Exectutes thes filtering using the efficients structures gained in the step above
    filter->start_filtering(log, bp);  

    // Completes the logging operation:
    log->end_measurement(complete_runningtime);
    log->set_program_inits(array_cnt, array);
    log->write_log_file();
    // Frees allocated memory:
    if(blast) {
      blast->free_hashProteins();
      blast->free_list_file_parse();
      //blast_parsing::close(blast);
      blast->free_arrNorm();
      blast->free_listTaxa();
      delete blast, blast=NULL;
    }
    blast_filtering::close(filter);
    cmd_list::close(cmd);
    log_builder::close(log);
  }
#ifdef USE_MPI
     MPI_Finalize();
#endif
  return 0;
}
int main (int argc, char *argv[])
{
  main_operation_turbo(argc, argv);
  return 0;
}
