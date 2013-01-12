#include "log_builder_main.h"

/**
   @brief Launches a test script using the log_builder library.
   @author Ole Kristian Ekseth (oekseth)
   @date 28.12.2011 by oekseth (cleanup)
**/
int main () {
  const char SEPERATOR = '_';
  const int CPU_TOT = 1;
  char * FILE_INPUT_NAME = "../all.blast";
    // Starts the logg tool gaining spesific measurements of time usage:
  log_builder *log = new log_builder(FILE_INPUT_NAME, SEPERATOR, BLOCK_FILE_READ_SIZE, 1024*1024, CPU_TOT);
  log->start_measurement(complete_runningtime);
  // Completes the logging operation:
  log->end_measurement(complete_runningtime);
  log->write_log_file();
  log_builder::close(log);
  return 0;
}
