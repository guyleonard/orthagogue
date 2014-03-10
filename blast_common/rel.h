/**
   @file
   @brief Defines a specific relation of proteins for a given 'mother protein'.
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
 */
#ifndef rel_h
#define rel_h
#include "types.h"
#include "tbb_libs.h"
#include "../configure.h"
#include "log_builder.h"
/**
   @class rel
   @brief Defines a specific relation of proteins for a given 'mother protein'.
   @ingroup blastfile_container 
   @author Ole Kristian Ekseth (oekseth)
   @date 19.08.2011 by oekseth (initial)
**/
class rel {
 public:
  //! The id of the outer protein.
  uint ind_out;
  //! The similarity score.
  float distance; 

  //! Prints the values set for this object.
  void  print() {
    printf("ind_out(%d), distance(%f)\n",
	   ind_out, distance);
  }
  //! Prints to stdout the values for the elements in this list:
  static void print_buffer(rel *buffer, uint size) {
    if(!buffer) {
      printf("(empty)\tbuffer at line %d in file %s is empty.\n", __LINE__, __FILE__);
    } else {
      printf("# Prints the buffer at line %d in file %s:\n", __LINE__, __FILE__);
      for(uint i = 0; i < size; i++) {
	buffer[i].print();
      }
    }
  }
  //! Prints the given protein pair
  void print(uint protein_in, char **name_in, char **name_out);  
  //! Used for verifying existence of data:
  bool is_equal(uint index) {return (ind_out == index);}
  //! Used for verifying existence of data:
  bool equals(rel arg) {return is_equal(arg);}
  //! Used for verifying existence of data:
  bool is_equal(rel arg) {
    if(arg.get_ind_out()==ind_out && (arg.get_distance() == distance))
      return true;
    else {
      return false;      
    }
  }

  //! Sets the data for the pair.
  void set_data(rel obj) {
    ind_out=obj.get_ind_out(), distance=obj.get_distance();
  }
  //! Sets the data for the pair.
  void set_data(uint ind_o, float sim_score, short int o_in) {
    if(false && !o_in) ; // In order to have a standarised interface.
    ind_out = ind_o;
    distance = sim_score;	
  }
  //! Sets the data given the input
  void set_data(uint ind_o, float sim_score, short int o_in, short int o_out) {
    if(false && !o_out) ; // In order to have a standarised interface.
    set_data(ind_o, sim_score, o_in);
  }

  //! @returns the index out.
  uint get_ind_out() {return ind_out;}
  //! @returns the similarity score.
  float get_distance(){return distance;}
  //! @returns the overlap for the inner.
  uint get_overlap_in(){return 0;}
  //! @returns the overlap for the outer.
  uint get_overlap_out(){return 0;}
  //! @returns the number of '0' (max) hits.
  uint get_cnt_max(){return 0;}
 
  /**
     @brief In order to know the object type when writing temporary files.
     @returns the unique class name     
  **/
  static char *get_class_name(){return "internal";}


  //! Reserves a list of these.
  static rel *reserve_list(uint buffer_size) {
    assert(buffer_size != UINT_MAX);
    if(buffer_size) {
      rel *ret = NULL;
      try {ret = new rel[buffer_size]();} 
      catch (std::exception& ba) {
	fprintf(stderr, "!!\t Tried reserving a \"rel\" class set of length=%u, [%s]:%s:%d\n",  buffer_size,  __FUNCTION__, __FILE__, __LINE__);
	if(!log_builder::catch_memory_exeception(buffer_size, __FUNCTION__, __FILE__, __LINE__, true)) {
	  fprintf(stderr, "!!\t An interesting error was discovered: %s."
		  "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		  "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
	}
      }
      log_builder::test_memory_condition_and_if_not_abort(ret !=NULL, __LINE__, __FILE__, __FUNCTION__);
      return ret;
    }
    else return NULL;
  }

  /**
     @brief Reserves- and initiates a list of rel objects.
     @param <buffer_size> The number of objects to produce, i.e. the length of the list returned.
     @param <buffer_in_mem_pos> A strict way of ensuring that the psotion starts at the start of it.
     @returns a list of objects of this class using the input as length of it
  **/
  static rel *init_list(uint buffer_size, loint &buffer_in_mem_pos) {
    if(buffer_size > 0) {
      rel *buffer = NULL;
      try {buffer = new rel[buffer_size]();} 
      catch (std::exception& ba) {
	fprintf(stderr, "!!\t Tried reserving a \"rel\" class set of length=%u, [%s]:%s:%d\n",  buffer_size,  __FUNCTION__, __FILE__, __LINE__);
	if(!log_builder::catch_memory_exeception(buffer_size, __FUNCTION__, __FILE__, __LINE__, true)) {
	  fprintf(stderr, "!!\t An interesting error was discovered: %s."
		  "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		  "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
	}
      }
      //      rel *buffer = new rel[buffer_size]();
      log_builder::test_memory_condition_and_if_not_abort(buffer!=NULL, __LINE__, __FILE__, __FUNCTION__);
      memset(buffer, 0, sizeof(rel)*buffer_size); // Prevents valgrind from complaining, as valgrind chekcs if its set (ie, initalized)
      for(loint i = 0; i<buffer_size; i++) buffer[i] = rel();
      buffer_in_mem_pos = 0; // A strict way of ensuring that the psotion starts at the start of it
      return buffer;
    } else return NULL;
  }

  /**
     @brief Init a list of this class.  
     @param <buffer_size> The number of objects to produce, i.e. the length of the list returned.
     @param <buffer> The list contaning objects of this type (rel) to initialise.
  **/
  static void init(loint buffer_size, rel *&buffer) {
    assert(buffer_size < (loint)UINT_MAX);
    init(0, (uint)buffer_size, buffer);
    //    for(uint i = 0; i<buffer_size; i++) buffer[i] = rel(_USE_BEST_BLAST_PAIR_SCORE);
  }
  /**
     @brief Init a list of this class.
     @param <start> The start position in the buffer to initialize.
     @param <end> The (last-1) position in the buffer to initialize.
     @param <buffer> The list contaning objects of this type (rel) to initialise.
     @remarks
     - Consitutes the standard method in this class for reserving memory for several objects of this type. 
  **/
  static void init(uint start, uint end, rel *&buffer) {
    assert(buffer);
    for(;start<end; start++) buffer[start] = rel();
  }
  
  //! Deletes a list containing elements of this class.
  static void delete_list(rel *&buffer, loint &pos_in_buffer, loint &size_of_buffer ) {
    close(buffer, pos_in_buffer, size_of_buffer);
  }
  //! Deletes a list containing elements of this class.
  static void close(rel *&buffer, loint &pos_in_buffer, loint &size_of_buffer) {
    if(buffer) {delete [] buffer, buffer = NULL;} pos_in_buffer = 0, size_of_buffer = 0;
  } 
  //! Copies a list of these from second argument to first.
  static void copy(loint size_destination, rel *&destination, rel *source, loint size_source) {
    assert(size_destination >= size_source);
    assert(source != NULL);
    if(true)	memcpy(destination, source, sizeof(rel)*size_source);
    else {
      for(loint i = 0; i< size_source; i++) destination[i] = source[i];
    }
  }
  //! Copies a list of these from second argument to first.
  static void copy(rel *&destination, rel *source, uint size) {
    if(true)	memcpy(destination, source, sizeof(rel)*size);
    else {
      for(mem_loc i = 0; i< size; i++) destination[i] = source[i];
    }
  }

  //! Copies a list of these from second argument to first at specific position.
  static void copy(loint size_destination, rel *&destination, uint start_pos_destination, rel *source, uint start_pos_source, uint size) {
    assert(size_destination >= (start_pos_source + size));
    assert(source != NULL);
    memcpy(destination+start_pos_destination, source+start_pos_source, sizeof(rel)*size); // Copies the data
  }

  /**
     @brief Enlarges the given list.
     @remarks If not set, allocates a new list.
     @param <buffer> The buffer to be extended at its lower end.
     @param <size_old> Used during DEBUG to assert the size is correct.
     @param <size_new> The length of the empty room to be made in front of the (new) extended list.
     @param <length_of_data_set> The number ob rel objects to be copied from the old buffer into the new one.
     @param <USE_BEST_BLAST_PAIR_SCORE> The default value to use when intiallizing the new elements
  **/
  static void enlarge(rel *&buffer, loint size_old, loint size_new, loint length_of_data_set, bool USE_BEST_BLAST_PAIR_SCORE) {
    if(false && USE_BEST_BLAST_PAIR_SCORE) ; // In order to use a standarised interface for enlarging the list.
    enlarge(buffer, size_old, size_new, length_of_data_set);
  }
  /**
     @brief Enlarges the given list.
     @remarks If not set, allocates a new list.
     @param <buffer> The buffer to be extended at its lower end.
     @param <size_old> Used during DEBUG to assert the size is correct.
     @param <size_new> The length of the empty room to be made in front of the (new) extended list.
     @param <length_of_data_set> The number ob rel objects to be copied from the old buffer into the new one.
  **/
  static void enlarge(rel *&buffer, loint size_old, loint size_new, loint length_of_data_set) {
    if(size_new > size_old) { // If it's not, there is not point in enlarging it.
      assert(size_old <= size_new);
      assert(length_of_data_set <= size_new);
      if(buffer) {
	rel *buff_temp = reserve_list(size_new); 
	log_builder::test_memory_condition_and_if_not_abort(buff_temp!=NULL, __LINE__, __FILE__, __FUNCTION__);
	memset(buff_temp, 0, sizeof(rel)*size_new); // Prevents valgrind from complaining, as valgrind chekcs if its set (ie, initalized)
	memcpy(buff_temp, buffer, sizeof(rel)*length_of_data_set);    
	delete [] buffer; buffer = buff_temp; // Updates the pointer with the new data.
      } else {
	buffer = reserve_list(size_new);
	memset(buffer, 0, sizeof(rel)*size_new); // Prevents valgrind from complaining, as valgrind chekcs if its set (ie, initalized)
      }
      if(length_of_data_set) {
	for(loint i = length_of_data_set; i < size_new; i++) {buffer[i]=rel();}
      }
    }

  }
  /**
     @brief Enlarges the given list an puts the 'current' at the back of it (reciprocal to what's normal).
     @remarks This procedure makes empty room in front of it, in order to enable copying of data from a file into it.
     @param <buffer> The buffer to be extended at its lower end.
     @param <new_size> Used during DEBUG to assert the size is correct.
     @param <start_pos> The length of the empty room to be made in front of the (new) extended list.
     @param <copy_length> The number ob rel objects to be copied from the old buffer into the new one.
     @param <_USE_BEST_BLAST_PAIR_SCORE> If set will not sum the scores when multiple scores are given for a pair (object of the p_rel type).
  **/
  static void enlarge_and_put_current_at_back(rel *&buffer, loint new_size, loint start_pos, loint copy_length, bool _USE_BEST_BLAST_PAIR_SCORE) {
    if(false && !(_USE_BEST_BLAST_PAIR_SCORE)) ; // Added to standarise the function.
    rel *buff_temp = reserve_list(new_size); //(rel*)malloc(sizeof(rel)*new_size);
    if(buff_temp) {
      assert(new_size <= (start_pos+copy_length)); // If procedure is used correct, this should hold. 
      for(loint i = 0; i < start_pos; i++) {buff_temp[i]=rel();}    
    // Inserts the buffer in the front of the file (i.e. in
      // order for the indexes to maintain their correct position pointers
      memcpy(buff_temp +start_pos, buffer, sizeof(rel)*copy_length);    
      delete [] buffer; 
      buffer = buff_temp;
    } else fprintf(stderr, "Not enough memory reserving %lld elements of type rel (having size %d). This message was generated at line %d in file%s, found in method %s. Please contact the developer at oekseth@gmail.com if this message is seen.\n", new_size, (int)sizeof(rel), __LINE__, __FILE__, __FUNCTION__);
  }
  

  /**
     @brief Writes the buffer given to a file:
     @param <buffer> The object to write. Do not change any of its data, i.e. the data is still kept in the buffer after the file writing.
     @param <file_name> The identifier later used reading the file now generated (or extended).
     @param <buffer_size> The number of relevant objects in the buffer (i.e. the 'size'-those allocated, but not used.)
     @param <file_length> The length of the current file.
  **/
  static void write_buffer(rel *&buffer, char *file_name, loint buffer_size, loint file_length) {
    assert(file_name);
    if(buffer!= NULL && (buffer_size > 0)) {
      FILE *file = NULL;
      if(file_length == 0) file = fopen(file_name, "wb");
      else file = fopen(file_name, "ab");
      if(file != NULL) {
	printf("writes %lld rel elements to file %s at line %d in soruce-file %s\n", buffer_size, file_name, __LINE__, __FILE__);
	fwrite(buffer, sizeof(rel), buffer_size, file);
	fclose(file);
      } else {
	fprintf(stderr, "!!\tUnable to open '%s' in %s at line %d at function %s: Aborts\n", file_name, __FILE__, __LINE__, __FUNCTION__);
	exit(2);
      }
    }
  }

  /**
     @brief Writes the buffer given to a file:
     @param <file_name> The identifier later used reading the file now generated (or extended).
     @param <buffer_one> The 1. object to write. Do not change any of its data, i.e. the data is still kept in the buffer after the file writing.
     @param <buffer_one_size> The number of relevant objects in the 1. buffer (i.e. the 'size'-those allocated, but not used.)
     @param <buffer_two> The 2. object to write. Do not change any of its data, i.e. the data is still kept in the buffer after the file writing.
     @param <buffer_two_size> The number of relevant objects in the 2. buffer (i.e. the 'size'-those allocated, but not used.)
     @param <file_length> The length of the current file.
     @attention The data is written in the same order as parameters given, ie, buffer_one is written BEFORE buffer_two.
   @author Ole Kristian Ekseth (oekseth).
  **/
  static void write_buffers(char *file_name, rel *buffer_one, uint buffer_one_size, rel *buffer_two, uint buffer_two_size, uint file_length) {    
    FILE *file = NULL;
    assert(file_name);
    if(file_length == 0) file = fopen(file_name, "wb");
    else file = fopen(file_name, "ab");
    if(file != NULL) {
      if(buffer_one && buffer_one_size > 0) {
	printf("writes %u rel elements to file %s at line %d in soruce-file %s\n", buffer_one_size, file_name, __LINE__, __FILE__);
	fwrite(buffer_one, sizeof(rel), buffer_one_size, file);
      }
      if(buffer_two && buffer_two_size > 0) {
	printf("writes %u rel elements to file %s at line %d in soruce-file %s\n", buffer_two_size, file_name, __LINE__, __FILE__);
	fwrite(buffer_two, sizeof(rel), buffer_two_size, file);
      }
      fclose(file);         
    } else {
      fprintf(stderr, "!!\tUnable to open '%s': Aborts\n", file_name);
      exit(2);
    }
  }

  /**
     @brief Reads data from a file into the buffer
     @param <file_name> The name of the file containing the data to get.
     @param <buffer> The buffer to put data into.
     @param <size> The number of objects to put into the given buffer. Requires that the size is at least equal to the number of elements allocated for the buffer.
     @return The buffer (same as given in the input);
  **/
  static rel *read_file(char *file_name, rel *buffer, loint size) {
    if(buffer != NULL) {
      FILE *file = fopen(file_name, "rb");
      if (file != NULL) {
	const loint size_read = fread(buffer, sizeof(rel), size, file); // to avoid the intial char
	if(size_read < size) init(size_read, size, buffer);
	assert(size_read <= size);
	fclose(file);
	return buffer;
      } else {
	fprintf(stderr, "!!\tCould not open '%s': input given has buffer_set=%d, rel_count(%lld) and total_size %d B; Aborts at line %d in %s\n", file_name, (int)(buffer!=NULL), size, (int)(sizeof(rel)*size), __LINE__, __FILE__);
	exit(2);
	return NULL;
      } 
    } else {
      fprintf(stderr, "!!\tCould not reserve %u amounts of class 'rel'; Aborts at line %d in %s!\n", (uint)size, __LINE__, __FILE__);
      exit(2);
    }
    return buffer;
  }

  /**
     @brief Reads a bunch of rel objects from the file given.
     @param <file_name> The name of the file to get the data from.
     @param <start_pos> The index in the file to read the first object from.
     @param <length> The number of objects to read.
     @returns A list containing the objects read.
  **/
  static rel* get_file_block(char *file_name, uint start_pos, uint length) {
    FILE *file = fopen(file_name, "rb");
    if (file != NULL) {
      rel *buffer_ret =  NULL; //       tbb::tbb_allocator<rel>().allocate(length);
      try {buffer_ret = tbb::tbb_allocator<rel>().allocate(length);} 
      catch (std::exception& ba) {
	if(!log_builder::catch_memory_exeception(length, __FUNCTION__, __FILE__, __LINE__)) {
	  fprintf(stderr, "!!\t An interesting error was discovered: %s."
		  "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		  "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
	}
      }
      //
      if(buffer_ret != NULL) {
	fseek(file, sizeof(rel)*start_pos, SEEK_SET);
	memset(buffer_ret, 0, sizeof(rel)*length);
	const uint read_length = fread(buffer_ret, sizeof(rel), length, file); // to avoid the intial char
	assert(read_length <= length);
	fclose(file);
	//	for(uint i = 0; i <length; i++) buffer_ret[i].print();
	return buffer_ret;
      } else { fprintf(stderr, "!!\tCould not allocate %u types of class rel; Aborts at line %d in %s!\n", length, __LINE__, __FILE__);
	exit(2);
      }
    } else {
      fprintf(stderr, "!!\tCould not open %s; Aborts at line %d in %s!\n", file_name, __LINE__, __FILE__);
      exit(2);
    }
    return NULL;
  }
  /**
     @brief Deallocates the object given.
     @param <block> The block to de-allocate.
     @param <size>  The number of elements to deallocate.
  **/
  static void free_file_block(rel *&block, uint size) {
    {//						if(block =! NULL) {
      tbb::tbb_allocator<rel>().deallocate((rel*)block, size);
      block = NULL;
    }
  }

  //! Reallocates the list, but make the room available in the start, instead of the latter as normal.
  static void realloc_making_free_space_in_front_of_list(rel *&buffer, loint old_size, loint new_size, loint &buffer_in_mem_pos) {
    if(!(old_size > new_size)) {
      assert(!(old_size > new_size));
    }    
    if(buffer) {
      assert(old_size != 0); // If it's correct that the buffer has data, this must also have value
      {
	const loint adjustment_length = new_size-old_size;
	rel *buff_temp = reserve_list(new_size); //(rel*)malloc(sizeof(rel)*new_size);
	if(buff_temp) {
	  //	  for(loint i = 0; i < adjustment_length; i++) {buff_temp[i]=rel();}    
	  for(loint i = 0; i < new_size; i++) {buff_temp[i]=rel();}    
	  assert((adjustment_length+old_size) <= new_size);
	  memcpy(buff_temp+adjustment_length, buffer, (sizeof(rel)*old_size));
	  delete [] buffer;     buffer = buff_temp;
	  buffer_in_mem_pos = adjustment_length; // TODO: verify this!
	} else fprintf(stderr, "Not enough memory reserving %lld elements of type rel (having size %d). This message was generated at line %d in file%s, found in method %s. Please contact the developer at oekseth@gmail.com if this message is seen.\n", new_size, (int)sizeof(rel), __LINE__, __FILE__, __FUNCTION__);
      }
    } else { // First set here:
      assert(old_size == 0); // If it's correct that the buffer does not have data, this should not hold data
      buffer = reserve_list(new_size);
      buffer_in_mem_pos = 0;
      init(0, new_size, buffer);
    }
  }
  /**
     @brief Enlarges the given block to its new size.
     @param <buffer> The block to extend.
     @param <old_size> The size of the buffer given as input.
     @param <new_size>  The number of elements in the new list
     @param <buffer_in_mem_pos> Included as a safeguard, setting the correct value.
  **/
  static void realloc_list(rel *&buffer, loint old_size, loint new_size, loint &buffer_in_mem_pos) {
    if(!(old_size > new_size)) {
      assert(!(old_size > new_size));
    }    
    if(buffer) {
      assert(old_size != 0); // If it's correct that the buffer has data, this must also have value
      rel::enlarge(buffer, old_size, new_size, /*length_of_data_set=*/buffer_in_mem_pos);
    } else { // First set here:
      assert(old_size == 0); // If it's correct that the buffer does not have data, this should not hold data
      buffer = reserve_list(new_size);
      buffer_in_mem_pos = 0;
      init(new_size, buffer);
    }
  }

  /**
     @brief Tests if the object given has the protein (search index) stated.
     @param <search_index> The index of the rightmost protein.
     @param <buffer> The list to search through.
     @param <start> The first index to do the searching at.
     @param <end> The (last-1) index to do the searching at.
     @return true if key-index is found in the data-set.
  **/
  static bool has_data(uint search_index, rel *buffer, uint start, uint end) {
    assert(start < end);
    for(uint k = start; k< end; k++) {
      if(buffer[k].is_equal(search_index)) {
	return true;
      }
    }
    return false;
  }

  //! The constructor.
  rel(uint _out, float _dist)   :
    ind_out(_out), distance(_dist) {}

  //! The constructor.
 rel() : ind_out(0), distance(0.0) {}

  //! The constructor.
   rel(bool temp) : ind_out(0), distance(0.0) {if(temp) ;}

  //! The assert method for this class
  static void assert_class(bool print_info);

};

/**
   @brief Defines a specific relation of proteins for a given 'mother protein'.
   @ingroup blastfile_container 
   @author Ole Kristian Ekseth (oekseth)
   @date 19.08.2011 by oekseth (initial)
**/
typedef class rel rel_t;

#endif
