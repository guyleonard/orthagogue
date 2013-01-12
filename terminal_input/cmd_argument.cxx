#include "cmd_argument.h"
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

//! Returns the correct usage of the enums, correspondning to the rules of validations
char *cmd_argument::get_correct_usage() {
  if(false) { /// currently deactivated.
    if(type == INTEGER)            return "Input expected is an integer: e.g. \"12\" or \"-12\".";
    else if(type == BOOLEAN)       return "";//Input expected is a et, activates it (i.e. sets to true if its false, alternative sets it to true)";
    else if(type == UINT_NOT_NULL) return "Input expected is an integer higher than 0: e.g. \"12\".";
    else if(type == CHAR_LIST)     return "Input expected is a name, e.g. \"all.blast\".";
    else if(type == CHAR_SINGLE)   return "Input expected is a single char, e.g. \"%c\".";
    else if(type == FLOAT)         return "Input expected is a floating value, e.g. \"12.02\".";
    else if(type == FOLDER_EXSISTING)        return "Input expected is the name of a directory on the current system, e.g. \"char\".";
    else if(type == FOLDER_NEW)     return "Input expected is a name of a new directory on the current system to be created, e.g. \"char\".";
  }   
  return "";//<..not specified..>";
}

char *get_name_of_type(ARGUMENT_TYPES type) {
  if(type == INTEGER)               return "integer";
  else if(type == BOOLEAN)          return "the boolean (0/1)";
  else if(type == UINT_NOT_NULL)    return "the unsigned integer (number greater than 0)";
  else if(type == CHAR_LIST)        return "string";
  else if(type == CHAR_SINGLE)      return "the letter (char)";
  else if(type == FLOAT)            return "float";
  else if(type == FOLDER_EXSISTING) return "string";
  else if(type == FOLDER_NEW)       return "string";
  else return "<undefined>";

}
//! Prints error message if the parsing did not correspond to the specification, found in the function 'get_correct_usage()'
void cmd_argument::print_wrong_value(char *identifier, char *value) {
  if(false && identifier) ;
  if(type == FOLDER_NEW || type == FOLDER_EXSISTING) {
    fprintf(stderr, "!!\tDid not find path %s argument given as either \"-%s %s\" or \"--%s %s\". The software will therefore abort: If \"ls -l %s\" does not return a file size, then your path is wrong. Alternatively contact the author at oekseth@gmail.com, including the following data: In function %s(..) at %s:%d\n", value, identifier_short, value, identifier_long, value, value, __FUNCTION__, __FILE__, __LINE__);
  } else {
    
    fprintf(stderr, "!!\tDid not find %s %s argument given as either \"-%s %s\" or \"--%s %s\". The input is therefore discarded. If this treatment of %s's is not understood, please contact the author at oekseth@gmail.com, including the following data: In function %s(..) at %s:%d\n", get_name_of_type(type), value, identifier_short, value, identifier_long, value, get_name_of_type(type), __FUNCTION__, __FILE__, __LINE__); 
  }
}

//! Returns true if argument has lphanumeric vlaues, assuming that the uinput has a '\0' at its end
bool cmd_argument::is_alpha_num(char *value) {
  bool is_alpha_num = true;
  char *start = value;
  char *end = strrchr(value, '\0');
  if(end != NULL)
    while(start != end) {
      if(!(isalnum(*start) || ispunct(*start))) return false;
      start++;
    }    
  else is_alpha_num=false; // do not have an alpha numeric value at its end
  if(is_alpha_num) return true; // at this point, the tests are succesful
  else {
    fprintf(stderr, "!!\tError: String '%s' were not in a pure alpha numeric mode: Email %s, sending a copy of this warning.\n", value, DEVELOPER_EMAIL);
    return false;
  }
}

/**! Returns true if argument has numeric values, assuming that the uinput has a '\0' at its end
   @Arg 'decimals' set to true if the number is a floating number
*/
bool cmd_argument::is_numeric(char *value, bool decimals) {
  bool is_alpha_num = true;
  char *start = value;
  char *end = strrchr(value, '\0');
  if(end != NULL) {
    while(start != end) {        // if the number is a floating number
      if(!(isdigit(*start) || (decimals && ispunct(*start)))) return false; 
      start++;
    }
  } else is_alpha_num=false; // do not have an alpha numeric value at its end
  if(is_alpha_num) return true; // at this point, the tests are succesful
  else {
    fprintf(stderr, "!!\tError: string '%s' were not in a pure  numeric mode: Email %s, sending a copy of this warning.\n", value, DEVELOPER_EMAIL);
    return false;
  }
}

/** Inititates the arguments for the constructor.
    @Changed: 17.02.2011 by oekseth
*/
void cmd_argument::init_arguments(ARGUMENT_TYPES type, void *_arg) {
  argument_b = NULL;    argument_f = NULL;
  argument_cc_s = NULL; argument_cc = NULL;  argument_i = NULL;

  // sets the correct parameter
  //    printf("type=%d for '%s'\n", type, identifier_short);
  if(type == INTEGER)                               {argument_i =    (int*)_arg;}
  else if(type == UINT_NOT_NULL)                    {argument_i =    (int*)_arg;}
  else if(type == BOOLEAN)                          {argument_b =    (bool*)_arg;}
  else if(type == CHAR_LIST)                        {argument_cc =   (char**)_arg; }
  else if(type == FOLDER_EXSISTING)                 {argument_cc =   (char**)_arg; }
  else if(type == FOLDER_NEW)                       {argument_cc =   (char**)_arg; }
  else if(type == CHAR_LIST)                        {argument_cc =   (char**)_arg; }
  else if(type == CHAR_SINGLE)                      {argument_cc_s = (char*)_arg;}
  else if(type == FLOAT)                            {argument_f =    (float*)_arg;  }
}
/** Prints the description of the arguments if
    ('sub_section_name!= NULL && 'sub_section_name == 'input_name'):
    - If this hold to true, then the 'sub_section_name' is set to NULL in order
    to avoid printing of this argument at a later stage
    @Argument 'char *input_name' is the name of the header for data to be printed.
    @Argument 'FILE *output_help' is the name of the file where data is to be printed to. Set to 'NULL' if no help to be printed
    @Argument 'FILE *file' is the name of the man page where data is to be printed to. Set to 'NULL'
    if no man page to be written.
    @Return: The 'sub_section_name' 
    @Changed: 17.02.2011 by oekseth
*/
char *cmd_argument::print_cmd_argument(char *input_name, FILE *output_help, FILE *file) {
  if(sub_section_name != NULL) {
    if(input_name == NULL) { // a new 'thing'
      if(output_help != NULL)       fprintf(output_help, "\n%s:\n", sub_section_name);
      if(file!= NULL) {
	fprintf(file, ".SH %s PARAMETERS\n\n", sub_section_name);
	fprintf(file, ".TP 1\n\n");
      }
      if(output_help != NULL){
	fprintf(output_help, "-%-4s --%-20s ", identifier_short, identifier_long);
	fprintf(output_help, "%s. %s\n", description, get_correct_usage());
      }
      if(file!= NULL) { 
	fprintf(file, ".B -%s, --%s\n", identifier_short, identifier_long);
	fprintf(file, " %s. %s\n\n", description, get_correct_usage());
      }
      char *temp = sub_section_name; 
      sub_section_name = NULL; // to avoid extra printing
      return temp; // returns the name
    } else if(0 == strcmp(sub_section_name, input_name)) { // they are equal
      if(output_help != NULL) {
	fprintf(output_help, "-%-4s --%-20s ", identifier_short, identifier_long);
	fprintf(output_help, "%s. %s\n", description, get_correct_usage());
      }
      if(file!= NULL) {
	fprintf(file, ".B -%s, --%s\n ", identifier_short, identifier_long);
	fprintf(file, "%s. %s\n\n", description, get_correct_usage());
      }
      sub_section_name = NULL; // to avoid extra printing
      return input_name; // returns the name
    } else return NULL;//sub_section_name; // if printed, returns NULL;
  } else return NULL;
}

/**! Inserts the value if it conforms to the input type set, else prints an error.
   @arg <char *identifier> is the name of identifier used: in order to produce a correct error message, if neccassary
   @arg <char *value> must not have a null value;
*/
void cmd_argument::insert_value(char *identifier, char *value) {
  data_is_set = true;
  bool input_in_wrong_format = false;
  if(type == INTEGER) 
    if(is_numeric(value, false)) {
      *argument_i = atoi(value);
    } else input_in_wrong_format = true;
  else if(type == BOOLEAN) {
    if(*argument_b == true) *argument_b = false; // Changes om the 
    else *argument_b = true;                // default value
  }
  else if(type == UINT_NOT_NULL) {
    if(is_numeric(value, false)) {
      const uint temp = atoi(value);
      if(temp > 0) {
	*argument_i = temp;
      } else input_in_wrong_format = true;
    } else input_in_wrong_format = true;
  }
  else if(type == CHAR_LIST)
    if(is_alpha_num(value)) *argument_cc = value;
    else input_in_wrong_format = true;
  else if(type == CHAR_SINGLE) {
    if(strlen(value) == 2)
      *argument_cc_s = value[1];
    else       *argument_cc_s = value[0];
  }
  else if(type == FLOAT)
    if(is_numeric(value, true))
      *argument_f = atof(value);
    else input_in_wrong_format = true;
  else if(type == FOLDER_NEW || type == FOLDER_EXSISTING) {       // Make a new directory in the current working directory 
    // Comment: The terminal transforms '~' to the home directory', therby a call to 'getenv' is not needed
    if(is_alpha_num(value)) {
      struct stat st;
      if(stat(value,&st) == 0) { // folder is present
	*argument_cc = value; //        printf(" /tmp is present\n");
      } else {
	if(type == FOLDER_NEW) {
	  if (-1 != mkdir(value, S_IRWXU)) { // 	 grant the owner read, write and execute authority 
	    *argument_cc = value;
	  } else   {
	    print_wrong_value(identifier, value);
	    data_is_set = false;
	  }
	} else { // Should be an exisitng folder, but does not exsist
	  print_wrong_value(identifier, value);
	  data_is_set = false;
	}
	if(!data_is_set) {
	  //! Sets an infomrative end-value, making eror messages later to be better undetood.
	  //! Note: May a small memory leake, but as the problem should abort, we do not regard this as a problem.
	  char *tmp = new char[400]; memset(tmp, '\0', 400);
	  sprintf(tmp, "given as either \"-%s %s\" or \"--%s %s\". The software will therefore abort: If \"ls -l %s\" does not return a file size, then your path is wrong. Alternatively contact the author at oekseth@gmail.com, including the following data: In function %s(..) at %s:%d\n", identifier_short, value, identifier_long, value, value, __FUNCTION__, __FILE__, __LINE__);
	  *argument_cc = tmp;
	}
      } 
    } else input_in_wrong_format = true;
  } else {
    fprintf(stderr, "!!\tError: input-type %d did not exist: Email %s, sending a copy of this warning.\n", type, DEVELOPER_EMAIL);
    data_is_set = false;
  }
  if(input_in_wrong_format) {
    data_is_set = false;
    print_wrong_value(identifier, value);
  }
}
/**! Inserts an argument if the text is found euqal to the input
   @Error: Prints error if the input does not conform to the spcification set. PRogram continues
   @Return: Returns true if the string found, and therby the variable updated.
*/
bool cmd_argument::insert_argument(char *_identifier_short, char *_identifier_long, char *value, bool &arg_a_boolean) {
  //! Updates the argument:
  if(type == BOOLEAN) arg_a_boolean = true; else arg_a_boolean = false;
  //! Sets the data:
  if(!data_is_set && (value != NULL || type == BOOLEAN)) {
    if(_identifier_short != NULL) {
      if(0==strcmp(identifier_short, _identifier_short)) { // the are equal
	insert_value(identifier_short, value);
	return true;
      }
    } else if (_identifier_long != NULL) {
      if(0==strcmp(identifier_long, _identifier_long)) { // the are equal
	insert_value(identifier_long, value);
	return true;
      }
    }
  }
  return false;  
}

/**! Returns true if both the short- and long identifier is not equal
   @Changed: 28.12.2010 by oekseth
*/
bool cmd_argument::is_not_equal(class cmd_argument arg) {
  if(0 == strcmp(identifier_short, arg.identifier_short))       return false;
  if(0 == strcmp(identifier_long, arg.identifier_long))         return false;
  return true;
}
//! Prints the value set
void cmd_argument::print_value_set() {
  if(true) printf("%s\t%s: ", identifier_short, description);
  else printf("-\t. Value set: '%s': ", description);
  if(type == INTEGER)             printf("%d", *argument_i);
  if(type == BOOLEAN)             printf("%d", *argument_b);
  else if(type == UINT_NOT_NULL)  printf("%d", *argument_i);
  else if(type == CHAR_LIST || type == FOLDER_EXSISTING || type == FOLDER_NEW)      printf("%s", *argument_cc);
  else if(type == CHAR_SINGLE)    printf("%c", *argument_cc_s);
  else if(type == FLOAT)          printf("%f", *argument_f);
  else printf("<not specified>");
  printf("\n");
}


cmd_argument::cmd_argument(char *_description, char *_identifier_short, char *_identifier_long,
			   ARGUMENT_TYPES _type, void *_arg, char *DEFAULT_OPTION_NAME, uint &DEFAULT_OPTION_NAME_COUNT) :
  description(_description),
  identifier_short(_identifier_short), identifier_long(_identifier_long),
  type(_type), data_is_set(false),
  sub_section_name(DEFAULT_OPTION_NAME) // the default option
{
  DEFAULT_OPTION_NAME_COUNT++; // increments, saying its set
  init_arguments(type, _arg);
  if(false)   print_cmd_constants();
}

cmd_argument::cmd_argument(char *_description, char *_identifier_short, char *_identifier_long,
			   ARGUMENT_TYPES _type, void *_arg, char *_sub_section_name) :
  description(_description),
  identifier_short(_identifier_short), identifier_long(_identifier_long),
  type(_type), data_is_set(false),
  sub_section_name(_sub_section_name) // the default option
{  init_arguments(type, _arg);
  if(false)   print_cmd_constants();
}

cmd_argument::cmd_argument() : description(NULL),
			       identifier_short(NULL), identifier_long(NULL),
			       type(NA), data_is_set(false)
{argument_f = NULL;           argument_cc_s = NULL; argument_cc = NULL;  argument_i = NULL;}

//! Tests the code with regard to the isenrtion of arguments, and theri validation
#ifdef assert_code
void cmd_argument::assert_insert_argument() {
  /*
    char *c;// = "test";
    class cmd_argument cl = cmd_argument("a test code stump", "l", "link", FOLDER_EXSISTING, &c);
    //  assert(false == cl.insert_argument("l", NULL, "debug_to"));
    assert(true  == cl.insert_argument("l", NULL, "debug"));
    assert(false == cl.insert_argument("p", NULL, NULL));
    assert(0 == strcmp(c, "debug"));

  
    // a list of chars:
    char *c;// = "test";
    class cmd_argument cl = cmd_argument("a test code stump", "l", "link", CHAR_LIST, &c);
    assert(true  == cl.insert_argument("l", NULL, "en_link"));
    assert(false == cl.insert_argument("p", NULL, "en_link"));
    assert(false == cl.insert_argument("p", NULL, NULL));
    assert(0 == strcmp(c, "en_link"));


    // a single char:
    char c_single;// = "test";
    class cmd_argument cl_s = cmd_argument("a test code stump", "l", "link", CHAR_SINGLE, &c_single);
    assert(true  == cl_s.insert_argument("l", NULL, "$_"));
    assert(false == cl_s.insert_argument("p", NULL, "$_"));
    assert(false == cl_s.insert_argument("p", NULL, NULL));
    assert(c_single == '_');


    // a float:
    float data_float;
    class cmd_argument cl2_f = cmd_argument("a test code stump", "m", "max_value", FLOAT, &data_float);
    assert(true  == cl2_f.insert_argument("m", NULL, "123.021"));
    assert(false == cl2_f.insert_argument("p", NULL, "123.021"));
    assert(false == cl2_f.insert_argument("p", NULL, NULL));
    assert(floor(data_float) == floor(123.021)); // uses 'floor(..)' to avoid changes in conversion, raising an error
  

    // an integer:
    int data;
    class cmd_argument cl2 = cmd_argument("a test code stump", "m", "max_value", INTEGER, &data);
    assert(true  == cl2.insert_argument("m", NULL, "123"));
    assert(false == cl2.insert_argument("p", NULL, "123"));
    assert(false == cl2.insert_argument("p", NULL, NULL));
    assert(data == 123);

    // a boolean:
    bool bol = false;
    class cmd_argument clZ = cmd_argument("a test code stump", "Z", "max_value", BOOLEAN, &bol);
    assert(true  == clZ.insert_argument("Z", NULL, NULL));
    assert(false == clZ.insert_argument("p", NULL, "123"));
    assert(false == clZ.insert_argument("p", NULL, NULL));
    assert(bol == true);  
  */
  if(false)   print_cmd_constants();
}
#endif

void cmd_argument::assert_class(const bool print_info) {
  const static char *class_name = "cmd_argument";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  cmd_argument::assert_insert_argument();
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

//! Frees the memory rservered
void cmd_argument::free_mem() {
  printf("function 'free_mem' at line %d in file %s not implemented\n", __LINE__, __FILE__);
}
