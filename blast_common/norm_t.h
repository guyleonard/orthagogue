#ifndef norm_t_h
#define norm_t_h
/**
   @file
   @brief Holds data about a taxon pair.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth).
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
 */
#include "taxa.h"
/**
   @class norm
   @brief Holds data about a taxon pair.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 26.08.2011 by oekseth (asserts)
   @date 24.12.2011 by oekseth (removed calls to 'extern' variables to ease the inclusion of thisclass as a libary)
   @date 10.07.2012 by oekseth (removed 'static' varaibles, and instead included them into fucntion-calls)
**/
struct norm {
  private:
  uint cnt;
  float sum;
  uint cnt_zero;
public:
  //! @return the total score accumulated for this object.
  float get_sum(const float max_input_value);
  //! @return the score for this object.
  float get_sum() {return sum;}
  //! @return the total number of pairs added in this object.
  uint get_cnt();  
  //! @return the total number of pairs added in this object.
  uint get_cnt_zeros() {return cnt_zero;}
  //! @return true if they are equal:
  bool is_equal(norm *obj);
  //! Inserts an given relation in this object.
  void insert(uint taxon_index_in, uint taxon_index_out, float arg, uint protein_in, uint protein_out, struct taxa *listTaxa, bool prot_pair_seen_before, const bool DEBUG_NORM);
  //! Inserts the relations into the structure using local coordinates
  void insert(uint taxon_in, uint taxon_out, float arg, uint protein_in, uint protein_out, struct taxa *listTaxa, const bool DEBUG_NORM);
  //! Takes the max input value found in the blast file pluss one,
  void changeZeroToValue(float max_input_pluss_one);
  //! Merges two structures given the data.
  void merge(uint _cnt, float _sum, uint _zero);
  //!  Merges two structures
  void merge(norm n);
  //! Print data given the params.
  void print(uint in, uint out);
  //! Prints the data.
  void print(struct norm **arrNorm, struct taxa *listTaxa);
  //! Print data given the params.
  static void print_basis(struct norm **arrNorm, int taxon_length);
  /**
     @brief Merging two structures.
     @param <taxon_length> The number of taxa in the collection.
     @param <arrNorm> The object who gets the data.
     @param <arrNorm_add> The object who gives away data
     @remark Do neither release any memory, or allocates any, due to assumption that the sizes of both equals.
   **/
  static void merge_basis(uint taxon_length, struct norm **arrNorm, struct norm **&arrNorm_add);
  /**
     @brief Merging several structures into one.
     @param <n_threads> The total number of squared **norm structures.
     @param <taxon_length> The number of taxa in the collection.
     @param <input> The object who gives away data
     @param <max_input_value> The maximal sim score found in the blast file.
     @param <PRINT_NORMALIXATION_BASIS> Used for initalizing new objects of class norm.
     @param <DEBUG_NORM> Used for initalizing new objects of class norm.
     @return The result object.
     @remark Do not deallocates memory, but transform those with value of '0' according to the rules given.
   **/
  static norm** merge_basis_zero(uint n_threads, uint taxon_length, struct norm ***input, float max_input_value, bool PRINT_NORMALIXATION_BASIS, bool DEBUG_NORM);
  //! Prints the data in form of the labels given in the input file:
  void print_label(struct taxa *listTaxa, uint i, uint out);
  //! Prints the data using the labels found in the input file
  static void print_labels(uint taxon_length, struct norm **arrNorm, struct taxa *listTaxa);
  //! @return the averaged value.
  float get_basis(const float max_input_value);
  
  //! @return the normative value.
  float get_normative_value() { if(cnt > 0 && sum > 0) return (sum/cnt); return 1.0;}

  //!  Builds the basis for the normative proecdure
  static float** build_basis(uint taxon_length, struct norm **arrNorm, bool PRINT_NORMALIXATION_BASIS, float max_input_value); 
  
  //! Allocates- and intiates memory for the 'norm_t' structure
  static struct norm **init_arrNorm(uint taxon_length, bool PRINT_NORMALIXATION_BASIS, bool DEBUG_NORM, float max_input_value);
  //! Allocates- and intiates memory for the 'norm_t' structure
  static struct norm ***init(uint taxon_length);

  //! Allocates- and intiates memory for the 'norm_t' structure
  static norm ***init_list(uint list_size, uint element_size, bool PRINT_NORMALIXATION_BASIS, bool DEBUG_NORM, float max_input_value);
  //! Deallocates memory for the object given.
  static void free_arrNorm(uint element_size, norm **&local_arrNorm);
  //! Deallocates memory for the object given.
  static void free_list(uint list_size, uint element_size, norm ***&local_arrNorm);
  //! Deallocates the object
  static void close(norm ***&obj, uint &taxon_length);
  //! The constructor.
  norm();
  //! The constructor.
  norm(uint _cnt, float _sum, uint _zero);

  //! The main test function for this class  
  static void assert_class(const bool print_info);
};

/**
  @brief Holds a pair of proteins.
  @ingroup blastfile_container
  @author Ole Kristian Ekseth (oekseth)
**/
typedef struct norm norm_t;
#endif
