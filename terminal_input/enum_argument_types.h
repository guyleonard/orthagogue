#ifndef enum_argument_types_h
#define enum_argument_types_h
/**
   @file
   @enum ARGUMENT_TYPES
   @brief Used for validation of inputs
   @ingroup terminal
   @author Ole Kristian Ekseth (oekseth)
**/
enum ARGUMENT_TYPES {
  UINT_NOT_NULL, // number highter than 0, without decimals. E.g <umber of cpu's>
  INTEGER, // number without decimals
  BOOLEAN, // either true of false: uses  'INTEGER' to store its data on
  CHAR_LIST, // List of chars.
  CHAR_SINGLE, // a single char
  FLOAT, // number with decimals
  FOLDER_NEW, // the name of a new folder
  FOLDER_EXSISTING, // the name of a new folder
  NA // not set
  //HELP, // a single item: lsits the input arguments
};
/**
   @brief Used for validation of inputs
   @ingroup terminal
   @author Ole Kristian Ekseth (oekseth)
**/
typedef enum ARGUMENT_TYPES ARGUMENT_TYPES_T;
#endif
