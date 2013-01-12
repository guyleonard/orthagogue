#ifndef mpi_transfer_norm_h
#define mpi_transfer_norm_h
#ifdef USE_MPI
/**
   @file
   @brief Holds data about a taxon pair, using the struct named norm.
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
//#include "norm_t.h"

/**
   @struct mpi_transfer_norm
   @brief Holds data about a taxon pair, using the struct named norm.
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
struct mpi_transfer_norm {
  //! 1d-value of the 2d-original-value.
  uint id;
  //! The number of elements: Corresponds the the variable of the same name in struct norm
  uint cnt;
  //! The sum of elements scores: Corresponds the the variable of the same name in struct norm
  float sum;
  //! The number of elements having '0' as value: Corresponds the the variable of the same name in struct norm
  uint cnt_zero;    
  /**
     @brief Reserves- and initiates a list of type mpi_transfer_norm.
     @return An initiated list of objects of type mpi_transfer_norm.
  **/
  static mpi_transfer_norm *init(uint length) {
    mpi_transfer_norm *list = new mpi_transfer_norm[length];
    log_builder::test_memory_condition_and_if_not_abort(list!=NULL, __LINE__, __FILE__, __FUNCTION__);
    for(uint i = 0; i < length; i++) {list[i]=mpi_transfer_norm();}
    return list;
  }

  //! Enlarges the given list to a double size:
  static void enlarge(mpi_transfer_norm *&list, uint &length) {
    assert(length); // According to our expectations.
    assert(list);
    //! Reserves length
    const uint new_length = 2*length;
    mpi_transfer_norm *new_list = new mpi_transfer_norm[new_length];
    log_builder::test_memory_condition_and_if_not_abort(new_list!=NULL, __LINE__, __FILE__, __FUNCTION__);
    //! Ensures consistency:
    for(uint i = length; i < new_length; i++) {list[i]=mpi_transfer_norm();}
    memcpy(new_list, list, sizeof(mpi_transfer_norm)*length);
#ifndef NDEBUG
    //! Ensures that the list is correctly copied; when included handy for inspections:
    for(uint i = 0; i < length; i++) {
      assert(list[i].id  == new_list[i].id);
      assert(list[i].cnt == new_list[i].cnt);
      assert(list[i].sum == new_list[i].sum);
      assert(list[i].cnt_zero == new_list[i].cnt_zero);
    }
    for(uint i = length; i < new_length; i++) {
      assert(list[i].id  == UINT_MAX);
      assert(list[i].cnt == 0);
      assert(list[i].sum == 0.0);
      assert(list[i].cnt_zero == 0);
    }
#endif
    //! Updates the parameters given:
    delete [] list; list = new_list; length = new_length;
  }
  //! Deallocates the list given:
  static void close(mpi_transfer_norm *&list, uint &length) {
    if(list) {
      delete [] list; list = NULL;
      length = 0;
    }
  }
  //! The constructor:
  mpi_transfer_norm() : id(UINT_MAX), cnt(0), sum(0.0), cnt_zero(0) {};
  //! The constructor:
  mpi_transfer_norm(uint _id, uint _cnt, float _sum, uint zero) :
    id(_id), cnt(_cnt), sum(_sum), cnt_zero(zero) {};
};

/**
   @brief Holds data about a taxon pair, using the struct named norm.
   @remarks Used for MPI transfer, and included with the macro variable USE_MPI
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
typedef struct mpi_transfer_norm mpi_transfer_norm_t;

#endif
#endif
