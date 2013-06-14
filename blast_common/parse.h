/**
   @file
   @brief Stores- and enumerates filtered blast data for each row of it.
   @ingroup parsing_ops
   @remark This object is later converted to a "file_parse" class structure
   using routines outside the scope of this class.
   @author Ole Kristian Ekseth (oekseth)
**/
#ifndef parse_h 
#define parse_h
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
#include "libs.h"
#include "types.h"
#include "log_builder.h"
#include "../configure.h"
/**
   @brief The column types in the blast file
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
**/
enum c_id {n_in, n_out, t_in, t_out}; 

/**
   @class Parse
   @brief Stores filtered blast data for each row of it.
   @ingroup parsing_ops
   @remark This object is later converted to a "file_parse" class structure
   using routines outside the scope of this class.
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 26.08.2011 by oekseth (asserts)
**/
class Parse{ // Used during parsing
 private:
 public:
  // Variables below public as accessing them direcly do not cause any possiblity of bugs (27.10.2011 by oekseth)
  //! The overlap value
  short int overlap_in; 
  //! The overlap value
  short int overlap_out; 
  //! The distance between the prots
  float distance; 
  private:
  //! name of the inner prot
  char *name_in;
  //! Name of the outer prot
  char *name_out;
  //! The taxa type of the inner
  char *taxon_in;
  //! The taxa type of the outer
  char *taxon_out;
  /**@brief Sets the name given the input       **/
  bool setName(char *start_pos, char *&name, uint index, uint *array, uint array_size);
  //! @return an indexed string based on a line from the blast file.
  uint *buildArray(char *start_pos, char *&pos_column_end, char seperator, uint &size_array, char *logical_end);
 public:
  /**
     @brief Returns the end of the column if the data contains all teh filds: retursn NULL if itdoes not contains all the data
     @param <first_column> True if it's the leftmost protein in the given blast-file line
     @param <start_pos> The memory position to start the work from.
     @param <index_name> The index in the blast file where the protein name is found
     @param <index_taxon> The index in the blast file where the taxon name is found
     @param <seperator>   Field separator in FASTA (blast) headers dividing the taxon from the protein label.
     @param <logical_end> The last address where data may be find; used to avoid errors.
     @param <array_size> An initial guess on how many columns there are in the blast file with regard to the name-label.
  **/
  char *set_names(const bool first_column, char *start_pos, uint index_name, uint index_taxon, char seperator, char *logical_end, uint array_size);

  /**
     @return the inner overlap.
  **/
  short int get_overlap_in(){return overlap_in;}
  /**
     @return the outer overlap.
   **/
  short int get_overlap_out(){return overlap_out;}
  //! Return true if they are equal.
  bool is_equal(char *str, c_id id);  
  //! Return the distance.
  float get_distance(){return distance;}
  //! @return the inner name
  char *getInnerName() {return name_in;}
    //! @return the inner name
  char *get_name_in();

  //! @return the outer name
  char *get_name_out() {return name_out;}
  //! @return the taxon inner name
  char *get_taxon_in(){return taxon_in;}
  //! @return the taxon outer name
  char *get_taxon_out(){return taxon_out;}

  //! @return the inner name's length
  loint get_name_in_length() {
    if(name_in) return strlen(name_in);
    else return 0;
  }
  //! @return the outer name's length
  loint get_name_out_length() {
    if(name_out) return strlen(name_out);
    else return 0;
  }
  //! @return the inner taxon's length
  loint get_taxon_in_length() {
    if(taxon_in) return strlen(taxon_in);
    else return 0;
  }
  //! @return the outer taxon's length
  loint get_taxon_out_length() {
    if(taxon_out) return strlen(taxon_out);
    else return 0;
  }

  //! @return the outer name
  char *getOuterName() {return name_out;}
  //! @return the taxon inner name
  char *getTaxonInName() {return taxon_in;}

  //! @return the taxon outer name
  char *getTaxonOutName() {return taxon_out;}
  //! @return the distance.
  float getDistance() {return distance;}

 private:
//! @return true if it's a char to be used
  bool is_legal_char(char c);
//! @return the correct length of the protein by adjusting the input given.
  loint get_correct_length(char *pos_start, loint length_name);

 //! Sets the name:
  void set_string(char *&string, char *pos_start, loint length_name);

 public:
//! Sets the inner name.
 void set_name_in(char *pos_start, loint length_name) {    
   //   char *name_in = &name_in[0];
   set_string(name_in, pos_start, length_name);
 }
  //! Sets the outer name.
  void set_name_out(char *pos_start, long long int length) {
    //    char *out = &name_out[0];
    set_string(name_out, pos_start, length);
  }
 
  //! Sets the taxon name in.
  void set_taxon_in(char *name, long long int length) {
    //    char *taxo = &taxon_in[0];
    set_string(taxon_in, name, length);
  }
  //! Sets the taxon name out.
  void set_taxon_out(char *name, long long int length) {    
    set_string(taxon_out, name, length);
  }
  //! Sets the distance.
  void setDistance(float dist) {distance = dist;}
  //! Sets the distance.
  void set_distance(float dist) {distance = dist;}
  //  void setOverlap(short int overl){ overlap = overl;}
  //! Sets the overlap in
  void setOverlap_in(short int overl){ overlap_in = overl;}
  //! Sets the overlap out
  void setOverlap_out(short int overl){ overlap_out = overl;}
  //! Sets the overlap in
  void set_overlap_in(short int overl){ overlap_in = overl;}
  //! Sets the overlap out
  void set_overlap_out(short int overl){ overlap_out = overl;}
  //! Prints to param the blast-format 'one'
  void debug_printBlastFormatOne(FILE *file_input, char SEPERATOR);
  //! Prints to param the blast-format 'two'
  void debug_printBlastFormatTwo(FILE *file_input, char SEPERATOR);
  //! Prints info about the class:
  void print_line(char seperator);  
  //! Prints info about the class
  void print(char SEPERATOR);
  //! Deallocates the memory reserved
  void free_memory() {
    delete [] name_in; delete [] name_out;
    name_in = NULL, name_out = NULL;
    delete [] taxon_in; delete [] taxon_out;
    taxon_in = NULL, taxon_out = NULL;
  }
  //! Constructor.
  Parse();
  //! Constructor.
  Parse(char *ni, char *no, char *ti, char *to, short int overl_in, short int overl_out, float di);

  //! The assert function for this class.
  static void assert_class(bool print_info);
};

/**
   @brief Stores filtered blast data for each row of it.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
 */
typedef class Parse Parse_t;
#endif
