#ifndef mpi_transfer_p_rel_h
#define mpi_transfer_p_rel_h
#ifdef USE_MPI
/**
   @file
   @brief Holds data about a protein pair, using the struct named p_rel.
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
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


/**
   @struct mpi_transfer_p_rel
   @brief Holds data about a protein pair, using the struct named p_rel.
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
struct mpi_transfer_p_rel {
  //! The overlap of the leftmost protein.
  uint overlap_in;
  //! The overlap of the rightmost protein.
  uint overlap_out;
  //! Counting the number of maxiumum-distance-hits (written as 0 in the blast input file)
  uint cnt_max; 
  //! The index of the outer protein.
  uint ind_out;
  //! The similarity score (high implies greater similarity).
  float distance; 

  /**
     @brief Build a datatype representing the p_rel class.
     @param <mpi_type> The datatype to fill.
  **/
  static void mpi_build_datatype(MPI_Datatype &mpi_type_info) {
    static const int cnt_types = 5;// The number of types to represent.
    int block_lengths[cnt_types] = {1,1,1,1,1};
    MPI_Datatype old_types[cnt_types] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_FLOAT};
    //! The location of each element;
    // TODO: Is the below "mpi initialisation" correct?
    MPI_Aint indices[cnt_types] = {0, sizeof(uint), sizeof(uint), sizeof(uint), sizeof(uint)};
    for(int i = 1; i < cnt_types;i++) {indices[i] +=indices[i-1];}
    //! Builds the mpi representation of the structure:
    MPI_Type_struct(cnt_types, block_lengths, indices, old_types, &mpi_type_info);
    MPI_Type_commit(&mpi_type_info);
  }

  /**
     @brief Reserves- and initiates a list of type mpi_transfer_p_rel.
     @return An initiated list of objects of type mpi_transfer_p_rel.
  **/
  static mpi_transfer_p_rel *init(uint length) {
    mpi_transfer_p_rel *list = new mpi_transfer_p_rel[length];
    log_builder::test_memory_condition_and_if_not_abort(list!=NULL, __LINE__, __FILE__, __FUNCTION__);
    for(uint i = 0; i < length; i++) {list[i]=mpi_transfer_p_rel();}
    return list;
  }

  //! @return true if they are equal:
  static bool is_equal(struct mpi_transfer_p_rel &obj1, struct mpi_transfer_p_rel &obj2) {
    bool all_was_equal = true;    
    if(obj1.overlap_in  != obj2.overlap_in)        all_was_equal = false;    
    else if(obj1.overlap_out != obj2.overlap_out)  all_was_equal = false;
    else if(obj1.cnt_max != obj2.cnt_max)          all_was_equal = false;
    else if(obj1.ind_out != obj2.ind_out)          all_was_equal = false;
    else if(obj1.distance != obj2.distance)        all_was_equal = false;
    return all_was_equal;
  }
  //! Enlarges the given list to a double size:
  static void enlarge(mpi_transfer_p_rel *&list, uint &length) {
    assert(length); // According to our expectations.
    assert(list);
    //! Reserves length
    const uint new_length = 2*length;
    mpi_transfer_p_rel *new_list = new mpi_transfer_p_rel[new_length];
    log_builder::test_memory_condition_and_if_not_abort(new_list!=NULL, __LINE__, __FILE__, __FUNCTION__);
    //! Ensures consistency:
    for(uint i = length; i < new_length; i++) {list[i]=mpi_transfer_p_rel();}
    memcpy(new_list, list, sizeof(mpi_transfer_p_rel)*length);
#ifndef NDEBUG
    //! Ensures that the list is correctly copied; when included handy for inspections:
    for(uint i = 0; i < length; i++) {
      assert(is_equal(list[i], new_list[i]));
    }
    mpi_transfer_p_rel empty = mpi_transfer_p_rel(); // Builds an empty object
    for(uint i = length; i < new_length; i++) {
      assert(is_equal(list[i], empty));
    }
#endif
    //! Updates the parameters given:
    delete [] list; list = new_list; length = new_length;
  }
  //! Deallocates the list given:
  static void close(mpi_transfer_p_rel *&list, uint &length) {
    if(list) {
      delete [] list; list = NULL;
      length = 0;
    }
  }
  //! The constructor:
  mpi_transfer_p_rel() : overlap_in(0), overlap_out(0),
			 cnt_max(), ind_out(0), distance(0.0)
  {};
  //! The constructor:
  mpi_transfer_p_rel(uint _out, float _dist, short int _overlap_in, short int _overlap_out, unsigned char _cnt_max) :
    overlap_in(_overlap_in), overlap_out(_overlap_out),
    cnt_max(_cnt_max), ind_out(_out), distance(_dist) {}
};

/**
   @brief Holds data about a protein pair, using the struct named p_rel.
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
typedef struct mpi_transfer_p_rel mpi_transfer_p_rel_t;

#endif
#endif
