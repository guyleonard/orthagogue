#ifndef enum_mcl_h
#define enum_mcl_h
/**
   @file
   @brief Defines the output types to be written as files.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
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
 * along with orthAgogue
. If not, see <http://www.gnu.org/licenses/>.
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
 * along with orthAgogue
 . If not, see <http://www.gnu.org/licenses/>.
**/

/**
   @enum mcl
   @brief Defines the output types to be written as files.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
*/
enum mcl {
  inpa,
  inpa_number,
  orth_inpa,
  orth_inpa_number,
  pair_orth,
  pair_orth_number,
  all,
  all_number,
  names_index, // holds the names of the proteins, and their corresponding index
  none // to be used as value when data is not to be specified
};


/**
   @brief Defines the output types to be written as files.
   @ingroup output_format
**/
typedef mcl mcl_t;

/**
   Size of enum
   @ingroup output_format
**/
static const uint mcl_t_size = 9;


/**
   "(9-1)/2" in order to get the name name file as well
   @ingroup output_format
*/
static const uint mcl_list_size = 4; 

/**
   List of those representing integers
   @ingroup output_format
**/
static const mcl_t list_mcl_number[] = {inpa_number, orth_inpa_number, pair_orth_number, all_number};

/**
   List of those representing strings
   @ingroup output_format
*/
static const mcl_t list_mcl[] = {inpa, orth_inpa, pair_orth, all};



/**
   @struct result_format
   @brief Defines the output types to be written as files.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
struct result_format {
  //! @return the formal mcl header format, to be used when building the list:
  static char *get_mcl_header_format() {
    char *mcl_header_format = "(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n(mcldoms\n";
    return mcl_header_format;
  }

  //! @return the length of the used mcl header.
  static uint get_length_of_mcl_header(uint size_proteins) {
    char string[1000]; memset(string, '\0', 1000);
    char *mcl_header_format = get_mcl_header_format();
    sprintf(string, mcl_header_format, size_proteins, size_proteins);
    const uint header_length = strlen(string);
    return header_length;
  }

  //! @return the formal stuff before content startmcl header format, to be used when building the list:
  static char *get_mcl_header_format_before_main() {
    char *mcl_header_format = "$\n)\n(mclmatrix\nbegin\n";
    return mcl_header_format;
  }

  /**
     @return the file names
     @remarks The naming convention follows:
     - Names with integers (numbers in this context) corresponds to mcl type, thereby .mcl'
     - Those files with strings (those other in this contect), the naming covention '.abc' is used
  **/
  static char *get_file_name_result(uint id) {
    const mcl_t t = (mcl_t)id;
    if(t == inpa)                   return "inparalogs.abc";
    else if(t == inpa_number)       return "inparalogs.mci";
    else if(t == orth_inpa)         return "co_orthologs.abc";
    else if(t == orth_inpa_number)  return "co_orthologs.mci";
    else if(t == pair_orth)         return "orthologs.abc";
    else if(t == pair_orth_number)  return "orthologs.mci";
    else if(t == all)               return "all.abc";
    else if(t == all_number)        return "all.mci";
    else if(t == names_index)       return "proteins.map";
    return NULL;
  }
private:
  /**
     @return the path to write/retrieve output result file from:
     @param <identifier> E.g. "all.abc" or "orthologs.mci"
     @param <FILE_BINARY_LOCATION> The result path given by the user.
     @param <SORT_ABC_DATA> If set, adds a the '.s' suffix to the file name
  **/
  static char *__get_storage_path_for_given_file_id(char *identifier, char *FILE_BINARY_LOCATION, bool SORT_ABC_DATA) {
    char *output_file = new char[1000];
    assert(output_file);
    memset(output_file, '\0', 1000);
    if(FILE_BINARY_LOCATION == NULL) FILE_BINARY_LOCATION = "./";
    if(FILE_BINARY_LOCATION != NULL && strlen(FILE_BINARY_LOCATION)>2) {
      if(FILE_BINARY_LOCATION[strlen(FILE_BINARY_LOCATION)-1] != '/') {
	if(SORT_ABC_DATA) {
	  sprintf(output_file, "%s/%s.s", FILE_BINARY_LOCATION, identifier);
	} else {sprintf(output_file, "%s/%s", FILE_BINARY_LOCATION, identifier);}
      } else {
	if(SORT_ABC_DATA) {
	  sprintf(output_file, "%s%s.s", FILE_BINARY_LOCATION, identifier);
	} else {
	  sprintf(output_file, "%s%s", FILE_BINARY_LOCATION, identifier);
	}
      }
    } else {
      if(SORT_ABC_DATA) {
	sprintf(output_file, "%s.s", identifier);
      } else {
	sprintf(output_file, "%s", identifier);
      }
    }
    
    return output_file;    
  }
public:
  
  /**
     @return the path to write/retrieve output result file from:
     @param <identifier> E.g. "all.abc" or "orthologs.mci"
     @param <FILE_BINARY_LOCATION> The result path given by the user.
     @param <SORT_ABC_DATA> If set, adds a the '.s' suffix to the file name
  **/
  static char *get_storage_path_for_given_file_id(char *identifier, char *FILE_BINARY_LOCATION, bool SORT_ABC_DATA=false) {
    if(SORT_ABC_DATA == false) {
      return __get_storage_path_for_given_file_id(identifier, FILE_BINARY_LOCATION, SORT_ABC_DATA);    
    } else {
      char *ret =  __get_storage_path_for_given_file_id(identifier, FILE_BINARY_LOCATION, SORT_ABC_DATA);
    struct stat sb;
    // Then try with the "*.s" sorted extension
    if(stat(ret, &sb) == -1) {
      delete [] ret;
      return __get_storage_path_for_given_file_id(identifier, FILE_BINARY_LOCATION, false);
    } else return ret;
    } 
  }
};
#endif

