#ifndef list_norm_h
#define list_norm_h
/**
   @file
   @brief Holds a matrix of norm objects.
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

#include "log_builder.h"
#include "norm_t.h"
/**
   @class list_norm
   @brief Holds a matrix of norm objects.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
class list_norm {
 private:
  norm_t ***list;
  int taxon_length;
  float max_input_value;
  bool PRINT_NORMALIZATION_BASIS;
  bool DEBUG_NORM;
 public:
  //! Should be used with care, as unwise deallocating of memory would be hard to track.
  void get_norm_2d(norm_t ***&obj) {obj = list;}
  //! @return the max input value.
  float get_max_input_value() {return max_input_value;}
  /**
     @brief Sets the data directly
     @remarks To be used with carefullness, as dat is not reserved, ie, function "free_mem()" must not be called.
  **/
  void set_data(norm ***_list, int _taxon_length, float _max_input_value);

  //! Sets the DEBUG_NORM value
  void set_DEBUG_NORM(const bool debug);
  /**
     @brief Prints the memory consumption for the object.
     @remarks For anaylsing the memory fingerprint.
  **/
  void print_getTotal_memoryConsumption_for_object(FILE *f);
  //! Test if condidtion is passed.
  static void test_condition_and_if_not_abort(bool condition, const int line, const char *file, const char *function);
  //! Inserts the relations into the structure using local coordinates
  void insert(uint taxon_in, uint taxon_out, float arg, uint protein_in, uint protein_out, struct taxa *listTaxa);
  //! Inserts the relations into the structure using local coordinates
  void insert(uint taxon_in, uint taxon_out, float arg, uint protein_in, uint protein_out, struct taxa *listTaxa, const bool prot_pair_seen_before);
  //! @return true if data is set.
  bool has_data() {return (list!= NULL);}
  //! @return true if data is set.
  bool has_data(int in);
  //! @return the object at the given position: If not set returns NULL
  norm_t *get_list_element(int in, int out);
  //! @return the object at the given position, and removes the local reference: If not set returns NULL
  norm_t *get_list_element_and_remove_reference(int in, int out);
  //! Prints the element given the coordinates
  void print_element(uint in, uint out);
  //! Print data given the params.
  void print_basis();
  //! Sums the values for the list of objects, updating the input parameters.
  void get_object_sum(loint &cnt, float &sum, uint &cnt_zero);
 private:
  //! Merges two objects
  void merge_basis(list_norm &obj, uint i, uint out);

 public:

  /**
     @brief Merging two structures.
     @param <obj> The object who gives away data.
     @remark 
     # Do neither release any memory, or allocates any, due to assumption that the sizes of both equals.
     # Validateion executed in the sub-call.
  **/
  void merge_basis(list_norm &obj);


 private:
  /*
  bool compare_visually_object(list_norm &obj) {
    //    cnt = 0, sum = 0.0, cnt_zero = 0;
    if(list) {
      for(uint i=0; i<(uint)taxon_length; i++) {
	if(list[i]) {
	  for(uint out=0; out<(uint)taxon_length; out++) {
	    norm_t *element = obj.get_list_element(i, out);
	    if(element) {
	      loint cnt = 0;
	      float sum = 0.0;
	      uint cnt_zero = 0;
	      loint cnt_arg = 0;
	      float sum_arg = 0.0;
	      uint cnt_zero_arg = 0;
	      if(list[i][out]) {
		cnt      += list[i][out]->get_cnt();
		sum      += list[i][out]->get_sum();
		cnt_zero += list[i][out]->get_cnt_zeros();
	      }
	    }
	  }
	}
      }
    }
    }*/
  //! Method below used to avoid the difficulties with round floating-pointer-numbers.
  bool is_not_different(float s1, float s2);
  //! Initiates the list using the already set data
  void init_list();
 public:
  /**
     @brief Merging several structures into one.
     @param <index_root> The part of the index to 'return, ie, to give a way and set to NULL, thereby not removing it when doing the 'cleaning' of the list_norm object.
     @param <n_threads> The total number of squared **norm structures.
     @param <input> The object who gives away data
     @param <max_input_value> The maximal sim score found in the blast file.
     @return The result object.
     @remark Do not deallocates memory, but transform those with value of '0' according to the rules given.
  **/
  void merge_basis_zero(uint index_root, uint n_threads, list_norm **&input, float max_input_value);
  //! Converts the '0's to a mutliplication of the inputs:
  void changeZeroToValue(float max_input_pluss_one);
  //! Prints the data in form of the labels given in the input file:
  void print_labels(struct taxa *listTaxa);    
  //! Prints the details of the normalisation procedure, given the argument list.
  void print_normalization_basis(float **arrAvgNorm);
  private:
  //! Converts the zero-value for the given index, setting the paramters to the updated values.
  void changeZeroToValue(uint taxon_in, uint taxon_out, float max_input_pluss_one, mem_loc &cnt, float &sum);  
  public:
  //!  Builds the basis for the normative proecdure
  float** build_basis();
    //! Changes this new type of matrix to the old format for backwards compability.
  norm_t **change_to_old_matrix();

#ifdef USE_MPI
  /**
     @brief Sends- and receives the list of norm_t* objects accross nodes.
  **/
  void mpi_make_data_consistent_accross_nodes(MPI_Comm comm, int myrank, int number_of_nodes);
#endif
  //! Allocates- and intiates memory for the 'norm_t' structure
  static list_norm **init_list(uint list_size, uint taxon_length, float max_input_value, bool PRINT_NORMALIXATION_BASIS, bool DEBUG_NORM);
  //! Allocates- and intiates memory for the 'norm_t' structure
  static list_norm **init_list(uint list_size);
  //! Deallocates memory for the list of list_norm objects.
  static void close(list_norm **&obj, uint list_size);
  //! Deallocates memory for the list of list_norm objects.
  static void close(list_norm *&obj);
  //! Deallocates memory for the object given.
  void free_mem();  
  //! The constructor.
  list_norm();
  //! The constructor.
  list_norm(uint t_length, float m_value, bool norm_basis, bool d_norm);  
};

/**
   @brief Holds a matrix of norm objects.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
  typedef class list_norm list_norm_t;
#endif
