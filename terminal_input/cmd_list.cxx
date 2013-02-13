#include "cmd_list.h"

/**! Builds a man page for the software
   @Reference: http://www.fnal.gov/docs/products/ups/ReferenceManual/html/manpages.html
   @Changed: 17.02.2011 by oekseth
*/ 
FILE * cmd_list::create_man_page_header() {
  if (true) {//if(argument_cnt > 0) { 
    if(SOFTWARE_NAME_MAN_PAGE != NULL) {
      //      printf("builds the manpages at line %d in cmd_list\n", __LINE__);
      SOFTWARE_NAME_MAN_PAGE = "~/temp.txt";
      // TODO: Ingores the local setting in the consttstants.h' file: Find a better way of doing i..
      char *home_dir = new char[70]; for(uint i=0; i<70;i++)home_dir[i]='\0';
#ifdef PROJECT_BINARY_DIR
      sprintf(home_dir, "%s/%s", FUNC_PROJECT_BINARY_DIR, "temp_man_page.txt");
      SOFTWARE_NAME_MAN_PAGE = home_dir;
#else
      sprintf(home_dir, "%s/%s", getenv("HOME"), "temp_man_page.txt");
      SOFTWARE_NAME_MAN_PAGE = home_dir;
#endif

      FILE *file = fopen(SOFTWARE_NAME_MAN_PAGE, "w");  
      if(file != NULL) {
	// Setting the data in the man page
	time_t rawtime = time_t();
	struct tm * timeinfo;
	timeinfo = localtime (&rawtime);
	if(false) printf("temp=%d\n", (int)(timeinfo->tm_sec));
	time_t t = time(0);
	struct tm* lt = localtime(&t);
	fprintf(file, ".TH %s 1 %02d/%02d/%04d\n\n", PROGRAM_NAME, lt->tm_mday, lt->tm_mon + 1, lt->tm_year + 1900);
	fprintf(file, ".SH NAME\n\n%s - %s\n\n", "orthAgogue", SOFTWARE_DESCRIPTION);
	fprintf(file, ".SH SYNOPSIS\n");	
	// TODO: Build a function in 'cmd_argument.h' to get the options on a line, like in 'man man'
	fprintf(file, ".B %s -i\n.I <file_name> [options] \n", PROGRAM_NAME); 

	if(SOFTWARE_REQUIREMENTS != NULL) 
	  fprintf(file, ".SH REQUIREMENTS\n\n%s\n\n", SOFTWARE_REQUIREMENTS);     
	if(SOFTWARE_AVAILABILITY != NULL) 
	  fprintf(file, ".SH AVAILABILITY\n\n%s\n\n", SOFTWARE_AVAILABILITY);
	if(SOFTWARE_DESCRIPTION_LONG != NULL)
	  fprintf(file, ".SH DESCRIPTION\n\n%s\n\n", SOFTWARE_DESCRIPTION_LONG);
	return file;
      } else fprintf(stderr, "!!\tFile name '%s' not possible to open: The path you have ssigned migh be wrong. Contact %s.\n", SOFTWARE_NAME_MAN_PAGE, DEVELOPER_EMAIL);
    } else fprintf(stderr, "!!\tFile name '%s' not possible to open. Contact %s.\n", SOFTWARE_NAME_MAN_PAGE, DEVELOPER_EMAIL);
  } else fprintf(stderr, "!!\tNo arguments provided. Function used wrong. Contact %s if this option is wishable!\n", DEVELOPER_EMAIL);
  return NULL;
}
/**! Builds the man page tail, and ends the operation of the man page building
   @Changed: 18.02.2011 by oekseth
*/
void cmd_list::create_man_page_tail(FILE *&file) {  
  if(file != NULL) {
    fprintf(file, ".SH EXAMPLES\n%s\n%s\n\n", SOFTWARE_EXAMPLES_INTRO, SOFTWARE_EXAMPLES);
   
    fprintf(file, ".SH COPYRIGHT\n\n%s\n\n", SOFTWARE_COPYRIGHT);
    fprintf(file, ".SH AUTHOR \n\n %s (email: %s).\n %s.\n", DEVELOPER_NAME, DEVELOPER_EMAIL, DEVELOPER_COMPANY);
    fprintf(file, ".SH SEE ALSO\n\n%s\n", SOFTWARE_SEE_ALSO);
    fclose(file);
    char cmd[300]; for(uint i=0; i<300;i++)cmd[i] = '\0';
    sprintf(cmd, "nroff -man %s > ~/orthAgogue.1", SOFTWARE_NAME_MAN_PAGE);//PROGRAM_NAME);
    assert(-1 != system(cmd));
    remove(SOFTWARE_NAME_MAN_PAGE); // removes the temporar file
    delete [] SOFTWARE_NAME_MAN_PAGE; SOFTWARE_NAME_MAN_PAGE = NULL;
  }
}
/** Prints the arguments
    @TODO: Create a man-page version of this
    @Changed: 11.03.2011 by oekseth. Edited the standard help hwen no arguments are provided
*/
void cmd_list::print_arguments(FILE *output_help, FILE *file) {
  char *sub_header = DEFAULT_OPTION_NAME;
  bool inserted_element = true;
  if(DEFAULT_OPTION_NAME_COUNT > 0) {
    if(output_help != NULL) fprintf(output_help, " %s PARAMETERS\n", DEFAULT_OPTION_NAME);
    if(file!= NULL) {
      fprintf(file, ".SH %s PARAMETERS\n\n", DEFAULT_OPTION_NAME);
      fprintf(file, ".TP 5\n\n"); 
    }
  } else sub_header = NULL;
  while(inserted_element == true) {
    inserted_element = false; // resets
    for(uint x =0; x<list_pos; x++) {
      if(sub_header == NULL) {
	sub_header = list[x].print_cmd_argument(sub_header, output_help, file);
	if(sub_header != NULL) { // a new sub header found
	  inserted_element = true; // an element was inserted
	}
      } else if(NULL != list[x].print_cmd_argument(sub_header, output_help, file)) {   
	inserted_element = true; // an element was inserted
      }
    }
    if(inserted_element == false) {
      if(sub_header != NULL) {
	sub_header = NULL; inserted_element = true; // make a new run
      } else 	inserted_element = false;
    }
  }
}

//! Print an error message when argument is not found.
void print_error_msg_cmd_not_found(char *&short_id, char *&long_id, char *string) {
    printf("!!\tNo data found for <%s>, using parameter <", string);
    if(long_id !=NULL) printf("%s", long_id);
    else if(short_id !=NULL) printf("%s", short_id);
    printf(">. Ignores.\n");
}
/**! Iterates through the argument list, setting the variable if found
   @Changed: 28.12.2010 by oekseth
*/
bool cmd_list::insert_argument(char *&short_id, char *&long_id, char *string, bool &remove_previous_references) {
  bool found = false;
  for(uint x =0; x<list_pos; x++) {
    if(!found){
      if(list[x].insert_argument(short_id, long_id, string, remove_previous_references)){
	found = true;
	DEFAULT_OPTION_NAME_COUNT++; 
      }
    }
  }	
  if(found) {long_id = NULL, short_id = NULL;} // resets
  return found;
}

/**! Binds the input arguments with variables
   @Changed: 17.02.2011 by oekseth
*/
void cmd_list::set_arguments(char **array, int array_size, bool create_man_page) {
  PROGRAM_NAME = array[0]; // if the name has changed from teh initail setting, its updated here for later reference
  FILE *file_man_page = NULL;
  if(create_man_page) file_man_page = create_man_page_header();
  if(array_size > 1) {
    char *long_id = NULL; char* short_id = NULL;
    bool remove_previous_references = false;
    for(uint i =1; i< (uint)array_size; i++) { // starts at one iot jump over the name of the program      
      if(array[i][0] == '-') { // start of short identifier
	if(strlen(array[i]) > 1) {
	  if(array[i][1] == '-') { // start of long identifier
	    //	    if(long_id != NULL) 
	    {
	      if(insert_argument(short_id, long_id, "\0", remove_previous_references)) {;// cmd_found=true;
	      }
	    }
	    long_id = array[i]+2; // avoid the '--'
	    //    id_updated = true;
	  } else {
	    //	    if(short_id != NULL) 
	    {
	      if(insert_argument(short_id, long_id, "\0", remove_previous_references)) {;// cmd_found=true;
	      }
	    }
	    short_id = array[i]+1; // avoid the '-'
	    //	    id_updated = true;
	  }
	}
      } else { // a value to set for the id
	if(insert_argument(short_id, long_id, array[i], remove_previous_references)) {;// cmd_found=true;
	}
      }
    }
    if(long_id != NULL) {if(insert_argument(short_id, long_id, "\0", remove_previous_references)) {// cmd_found=true;
      }}
    if(short_id != NULL) {if(insert_argument(short_id, long_id, "\0", remove_previous_references)) {// cmd_found=true;
      }}
    //    if(!cmd_found) {print_error_msg_cmd_not_found(short_id, long_id, "\0");}
  }
  if (create_man_page || (array_size <2) || FILE_INPUT_NAME == NULL)  { // file name not given or zero arguments given. PRints default message:   
    if(array_size < 2) fprintf(stdout, "!!\tBlast file not specified.\n\nUSAGE:\northAgogue -i <file_path> [-t][-p][-s][-O][-P][-S][-A][-w][-b][-c][-u|-e][-o][-C][m]\n");
    //    if(array_size < 2) fprintf(stdout, "!!\tBlast file not specified.\n\nUSAGE:\northAgogue -i <file_path> [-p-t-s-O-P-sABC-A-us-bho-mp-c-e-o-nii-L-nf-N-pd]\n");
    if(array_size < 2)    print_arguments(stdout, file_man_page);

    fprintf(stdout, "\nEXAMPLES\n1. orthagogue -i myblast.out -e 6 -O myoutdir # excludes protein pairs with e-values above 1e-06 and saves results in 'myoutdir'.\n");
    fprintf(stdout, "2. orthagogue -i myblast.out -e 6 -o 50 -O myoutdir # the same as above but excluding protein pairs with the overlap below 50%%.\n");
    fprintf(stdout, "3. orthagogue -i myblast.out -u -o 50 -O myoutdir # the same as above but without filtering by e-values and using BLAST scores instead of e-values in order to resolve HSPs with the '0.0' e-value; the required e-value cutoff should be set while running BLAST.\n");
    fprintf(stdout, "4. orthagogue -i myblast.out -b -e 6 -O myoutdir # the same as 1 but using only the best HSP for any pair of proteins (OrthoMCL emulation).\n");
    if(array_size < 2)   fprintf(stdout, "\nSee also https://code.google.com/p/orthagogue/\n\nThe software was developed by O.K. Ekseth under supervison of Dr. V.N. Mironov. Questions to be forwarded to orthagogue-issue-tracker@googlegroups.com\n");
// oekseth@gmail.com\n");
    //
    if(create_man_page) create_man_page_tail(file_man_page);
    if(array_size < 2) exit(2); // Aborts the execution
      
  } else if(create_man_page) {
    print_arguments(stdout, file_man_page);
    //    print_arguments(NULL, file_man_page);
    create_man_page_tail(file_man_page);
  }
  
}

//! Called when arguments are to be added
void cmd_list::add_cmd_argument(class cmd_argument arg) {
  bool not_unique = false;
  uint not_unique_index = UINT_MAX;
  for(uint x =0; x<list_pos; x++)
    if(!(list[x].is_not_equal(arg))) {not_unique = true; not_unique_index = x;}
  if(not_unique)  {
    fprintf(stderr, "!!\tThe id of the short identifier is not unique: Ingores \""); 
    arg.print_cmd_argument(NULL, stdout, NULL);// printf("\"!\n");
    list[not_unique_index].print_cmd_argument(NULL, stdout, NULL); printf("\"!\n");
    assert(false);
  } else {
    if(list_pos == list_size) { // enlarges the list
      const uint new_size = 2*list_size;
      class cmd_argument *temp = new cmd_argument[new_size]();
      memcpy(temp, list, sizeof(cmd_argument)*list_size);
      list_size = new_size;
      delete [] list;
      list = temp;
    }
    list[list_pos] = arg;
    list_pos++;
  }
}

//! The main test function for this class  
void cmd_list::assert_class(const bool print_info) {
  const static char *class_name = "cmd_list";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code 
  uint cnt = 0;   class cmd_list obj = cmd_list("STATUS", cnt, "orthAgogue_man");
  class cmd_argument cmd;
  static const uint list_size = 3;
  //! Holdres for variables we test:
  bool list_bool[list_size];   float list_float[list_size];   uint list_uint[list_size];
  //! Builds a set of parameters:
  uint i = 0; cmd = cmd_argument("bool", "b0", "bool0", BOOLEAN, &list_bool[i], "TEST");
  obj.add_cmd_argument(cmd); list_bool[i] = false;
  i = 1; cmd = cmd_argument("bool", "b1", "bool1", BOOLEAN, &list_bool[i], "TEST");
  obj.add_cmd_argument(cmd); list_bool[i] = false;
  i = 2; cmd = cmd_argument("bool", "b2", "bool2", BOOLEAN, &list_bool[i], "TEST");
  obj.add_cmd_argument(cmd); list_bool[i] = false;

  i = 0; cmd = cmd_argument("float", "f0", "float0", FLOAT, &list_float[i], "TEST");
  obj.add_cmd_argument(cmd); list_float[i] = 0; //(float)i;
  i = 1; cmd = cmd_argument("float", "f1", "float1", FLOAT, &list_float[i], "TEST");
  obj.add_cmd_argument(cmd); list_float[i] = 0; //(float)i;
  i = 2; cmd = cmd_argument("float", "f2", "float2", FLOAT, &list_float[i], "TEST");
  obj.add_cmd_argument(cmd); list_float[i] = 0; //(float)i;

  i = 0; cmd = cmd_argument("uint", "u0", "uint0", UINT_NOT_NULL, &list_uint[i], "TEST");
  obj.add_cmd_argument(cmd); list_uint[i] = 0; //(uint)i+1;
  i = 1; cmd = cmd_argument("uint", "u1", "uint1", UINT_NOT_NULL, &list_uint[i], "TEST");
  obj.add_cmd_argument(cmd); list_uint[i] = 0; //(uint)i+1;
  i = 2; cmd = cmd_argument("uint", "u2", "uint2", UINT_NOT_NULL, &list_uint[i], "TEST");
  obj.add_cmd_argument(cmd); list_uint[i] = 0; //(uint)i+1;


  //! Starts the parsing of it
  //! Note: if the below size is worngly set, errors are caused!
  static const uint array_1_size = 16;
  //! Note: The <program_name> must be the first arugment.
  char *array_1[] = {"<program_name>", "-b0", "-b1","-b2", "-u0", "1", "-u1", "2", "-u2", "3", "-f0", "0", "-f1", "1", "-f2", "2"};
  obj.set_arguments(array_1, array_1_size, false);  

  for(uint i = 0; i < list_size; i++) {
    assert(list_uint[i] == (i+1));
  }
  for(uint i = 0; i < list_size; i++) {
    assert(list_bool[i] == true);
  }

  for(uint i = 0; i < list_size; i++) {
    assert(list_float[i] == (float)i);
  }
  obj.free_mem();
#endif
  if(false)   print_cmd_constants();
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}


