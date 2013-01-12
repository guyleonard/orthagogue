#ifndef pre_read_blast_h
#define pre_read_blast_h
/**
   @file
   @brief Describes the blast preprocessing.
   @ingroup mpi_blast
 **/
#include "../configure.h"
#include "libs.h"
#include "types.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "read_file.h"
#include "tsettings_input.h"
#include "string_section.h"
#include "taxon_list_t.h"
#include "log_builder.h"
/**
   @class pre_read_blast
   @brief Describes the blast preprocessing.
   @ingroup mpi_blast
   @author Ole Kristian Ekseth (oekseth)
   @date 02.01.2012 by oekseth (initial)
**/
class pre_read_blast {
 private:
  uint disk_buffer_size;
  const uint MACHINE_CNT;
  tsettings_input_t blast_settings;
  loint *file_start_pos;
  loint *file_end_pos;
  //! @return the file pointer.
  FILE *get_file();
  //! @return true if exact positions will need to be calculated.
  bool calculate_exact_positions(const loint file_size, const loint approx_size_each);
 public:
  //! @return the file size.
  static loint get_file_size(char *file) ;
  // -----------------------------------------------------------
  // Standard proecedures for getting the results:

  //! @return the starting position in the file for the machine id given.
  loint get_start_pos(uint i);
  //! @return the end position in the file for the machine id given.
  loint get_end_pos(uint i);
  //! @return number fo chars to read in the file for the machine id given.
  loint get_read_length(uint i);
  //! Prints the elmement given the index as param
  void print_element(uint index); 
  // -----------------------------------------------------------
  // The interesting part of this class (public):

  /**
     @brief Produces a safe division of the blast file (into blocks).
     @remarks Procedure is:
     -# Tests some general properties set by the constructor, e.g. file exsistence.
     -# Reserves memory and checks if exact positions must be caluclated.
     -# Sets the approximate size for each machine id.
     -# Blocks are tests, ensuring there's a shift from one type of left-most protein to another. This avoids the need for the merging procedure in the mpi-implemntation having to consider overlapping regions.
     @author Ole Kristian Ekseth (oekseth)
     @date 02.01.2012 by oekseth (initial)
     @date 03.01.2012 by oekseth (asserts)
  **/
  void start_preprocessing(int myrank, char seperator);

  // -----------------------------------------------------------
  // Initiators and de-allocations procedures:
  //! Initiates an object of type pre_read_blast.
  static pre_read_blast *init(uint disk_buffer_size, uint MACHINE_CNT, tsettings_input_t blast_settings);
  //! De-allocates the memory reserved the object of type pre_read_blast.
  static void close(pre_read_blast *&obj);
  //! De-allocates the memory reserved for this object.
  void free_memory();
  /**
     @brief The constructor:
     @param <_disk_buffer_size> The number of chars to approximate read for each 'run'
     @param <_MACHINE_CNT> The number of machines (nodes) to use for the MPI operation.
     @param <_blast_settings> Defines the properties for the blast file.
     @author Ole Kristian Ekseth (oekseth)
  */
  pre_read_blast(uint _disk_buffer_size, uint _MACHINE_CNT,  tsettings_input_t _blast_settings);


  // -----------------------------------------------------------
  // Methods for testing:
  
  //! The test function for the private parts of this class.
  void assert_private_parts();

  /**
     @brief The main test function for this class.
     @remarks This method includes:
     - Formalized tests for verification that a block starts at a line.
     - Valgrind used verifying memory usage.
     - Examples of how this class may be used.
     @author Ole Kristian Ekseth (oekseth)
     @date 03.01.2012 by oekseth.
  **/
  static void assert_class(const bool print_info); 
};

/**
   @brief Describes the blast preprocessing.
   @ingroup mpi_blast
   @author Ole Kristian Ekseth (oekseth)
   @date 02.01.2012 by oekseth (initial)
**/
typedef class pre_read_blast pre_read_blast_t;
#endif
