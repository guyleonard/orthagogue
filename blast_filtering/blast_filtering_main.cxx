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
#include "blast_filtering_main.h"


/**
   ---------------------------------------------------------------------------------------
   @Name: main.cxx -- This file is regarded as the launcher of the software, implying using the librariy  'blast_filtering'.
   -- In order to understand the above library, this file might be of use the the developer.
   @Author: Ole Kristian Ekseth (oekseth)
   @Last_Major_Update: 27.12.2011 by oekseth (cleanup)
   ---------------------------------------------------------------------------------------
**/
int main(int argc, char *argv[]) {
  //blast_filtering::assert_class(true);
  /*
    #ifdef assert_code
    const bool print_info = true;
    const static char *class_name = "main";
    if(print_info) printf("--\tStarts testing class %s\n", class_name);
    #include "pipe_struct.h"
    #include "pipe_bucket.h"
    #include "pipe_merge.h"
    #include "pipe_norm.h"
    #include "pipe_write.h"
#include "build_string.h"

#include "pipe_ortholog_inparalog.h"
pipe_binary::assert_class(print_info);
pipe_struct::assert_class(print_info);
pipe_bucket::assert_class(print_info);
pipe_merge::assert_class(print_info);
pipe_norm::assert_class(print_info);
pipe_write::assert_class(print_info);
#endif
  */
  char *DEFAULT_OPTION_NAME = "OPTIONS"; // Used when no option is specified
  uint DEFAULT_OPTION_NAME_COUNT = 0; // counts the number of times the 'options' is used
  char *FILE_INPUT_NAME ="temp.raritet.txt";
  char SEPERATOR = '_'; const int CPU_TOT = 1;
  // Constructs the object taking input from the terminal:
  cmd_list *cmd =blast_filtering::init_cmd_list(DEFAULT_OPTION_NAME, DEFAULT_OPTION_NAME_COUNT, FILE_INPUT_NAME);
  // Initializes the library we are interested in:
  blast_filtering * filter = new blast_filtering(cmd);
  log_builder *log =new log_builder(FILE_INPUT_NAME, SEPERATOR, BLOCK_FILE_READ_SIZE, 1024*1024, CPU_TOT);
  // Maps internal variables with the terminal input:
  blast_filtering::build_cmd_list(cmd, argc, argv);
  filter->print_class_info();

  // Exectutes thes filtering using the efficients structures gained in the step above
  bp_container_t bp; 
  filter->start_filtering(log, bp);  // bp not set. Aborts

  // Completes the logging operation:
  log->end_measurement(complete_runningtime);
  log->write_log_file();

  // Frees allocated memory:
  blast_filtering::close(filter);
  cmd_list::close(cmd);
  log_builder::close(log);
  return 0;
}
