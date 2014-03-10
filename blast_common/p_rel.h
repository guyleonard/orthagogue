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
#ifndef p_rel_h
#define p_rel_h
#include "types.h"
#include "tbb_libs.h"
#include "../configure.h"
#include "log_builder.h"
/**
   @file
   @class p_rel
   @brief Like class 'rel', but with some extension for storing temporary blast-data.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
   @date 19.08.2011 by oekseth (initial)
**/
class p_rel {
 private:
  uint overlap_in;
  uint overlap_out;
  uint cnt_max; // Counting the number of maxiumum-distance-hits (written as 0 in the blast input file)
  uint ind_out; 
  float distance; 
  //  bool USE_BEST_BLAST_PAIR_SCORE; // If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
 public:
  //  //! Return true if the OrthoMcl default setting is used.
  //  bool uses_best_blast_pair_score() {return USE_BEST_BLAST_PAIR_SCORE;}
  /**
     @brief In order to know the object type when writing temporary files.
     @returns the unique class name     
  **/
  static char *get_class_name(){return "parse";}
  //! Used for verifying existence of data:
  bool is_equal(uint index) {return (ind_out == index);}
  //! Used for verifying existence of data:
  bool equals(p_rel arg) {return is_equal(arg);}
  //! Used for verifying existence of data:
  bool is_equal(p_rel arg) {
    if(arg.get_ind_out()==ind_out && (arg.get_distance() == distance) && (arg.get_cnt_max() == cnt_max) /*&& (arg.get_overlap_out() == overlap_out) */&& (arg.get_overlap_in() == overlap_in)
       //&& (arg.uses_best_blast_pair_score() == USE_BEST_BLAST_PAIR_SCORE)
       )
      return true;
    else {
      printf("err\tthis.overlap_in(%d), arg.overlpa_in(%d)  ", overlap_in, arg.get_overlap_in());
      printf("err\tthis.overlap_out(%d), arg.overlpa_out(%d)\n", overlap_out, arg.get_overlap_out());
      return false;
      
    }
  }
  //! Inits the data to empty values.
  void init() {ind_out=0, distance=0.0, overlap_in = 0, overlap_out = 0, cnt_max = 0; }
  //! Sets the outer index.
  void set_ind_out(uint _ind) {ind_out = _ind;}
  //! Sets the similiarity score for this object.
  void set_distance(float d) {distance = d;}
  //! Sets the overlap for the inner protein.
  void set_overlap_in(short int d) {overlap_in = d;}
  //! Sets the outer overlap value.
  void set_overlap_out(short int d) {overlap_out = d;}
  //! Sets the number of maximal scores (0-hits).
  void set_cnt_max(uint d) {cnt_max = (char)d;}
  //! Sets the data for the pair.
  void set_data(p_rel obj) {
    ind_out=obj.get_ind_out(), distance=obj.get_distance();
    overlap_in = obj.get_overlap_in(), overlap_out = obj.get_overlap_out(), cnt_max = obj.get_cnt_max();
  }
  //! Sets the data for the pair.
  void set_data(uint ind_o, float sim_score, short int o_in) {
    ind_out = ind_o;
    if(sim_score == 0.0) distance = 0.0, cnt_max =1; // The first protein met
    else cnt_max = 0, distance = sim_score;	
    overlap_in = o_in;
  }
  //! Sets the data given the input
  void set_data(uint ind_o, float sim_score, short int o_in, short int o_out) {
    set_data(ind_o, sim_score, o_in);
    overlap_out = o_out;
  }

  
  /**
     @brief Sets the scores for this object (p_rel pair).
     @remarks if variable USE_BEST_BLAST_PAIR_SCORE is set, do not sum together the different scores.
     @param <sim_score> The similarity value to use
     @param <o_in> The overlap for the leftmost protein.
     @param <USE_BEST_BLAST_PAIR_SCORE> If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
  **/
  bool increment_data(float sim_score, uint o_in, const bool USE_BEST_BLAST_PAIR_SCORE) {
    if(USE_BEST_BLAST_PAIR_SCORE) {
      if(cnt_max == 0) {
	if(sim_score == 0) {cnt_max = 1, overlap_in = o_in; return true;}
	else if(sim_score > distance) {distance = sim_score, overlap_in = o_in; return true;}
      }
      return false;
    } else {
      if(sim_score == 0.0) cnt_max++;
      else distance += sim_score;
      overlap_in += o_in;
      return true;
    }
  }
  /**
     @brief Sets the scores for this object (p_rel pair).
     @remarks if variable USE_BEST_BLAST_PAIR_SCORE is set, do not sum together the different scores.
     @param <sim_score> The similarity value to use
     @param <o_in> The overlap for the leftmost protein.
     @param <o_out> The overlap for the rightmost protein.
     @param <USE_BEST_BLAST_PAIR_SCORE> If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
  **/
  void increment_data(float sim_score, uint o_in, uint o_out, const bool USE_BEST_BLAST_PAIR_SCORE) {
    const bool sim_score_updated = increment_data(sim_score, o_in, USE_BEST_BLAST_PAIR_SCORE);
    if(USE_BEST_BLAST_PAIR_SCORE) {
      if(sim_score_updated) overlap_out = o_out;
    } else overlap_out += o_out;
  }
  //! @returns the index out.
  uint get_ind_out() {return ind_out;}
  //! @returns the similarity score.
  float get_distance(){return distance;}
  //! @returns the overlap for the inner.
  uint get_overlap_in(){return overlap_in;}
  //! @returns the overlap for the outer.
  uint get_overlap_out(){return overlap_out;}
  //! @returns the number of '0' (max) hits.
  uint get_cnt_max(){return (uint)cnt_max;}

  /**
     @return The similarity score, based upon the data: used in pipe_binary.cpp.
     @date 01.01.2011 by oekseth
  */
  float getSimScore(float max_input_value) {
    return (distance + (float)((int)cnt_max*(max_input_value+1)));
  }
  //! Prints the values set for this object.
  void  print() {
    printf("ind_out(%d), distance(%f), overlap_in(%d), overlap_out(%d), cnt_max(%d)\n",
	   ind_out, distance, overlap_in, overlap_out, (int)cnt_max);
  }
  //! Prints the values set for this object.
  void print(uint index) {
    printf("[%u]\tind_out(%d), distance(%f), overlap_in(%d), overlap_out(%d), cnt_max(%d)\n",
	   index, ind_out, distance, overlap_in, overlap_out, (int)cnt_max);
  }
  //! Prints to stdout the values for the elements in this list:
  static void print_buffer(p_rel *buffer, uint size) {
    if(!buffer) {
      printf("(empty)\tbuffer at line %d in file %s is empty.\n", __LINE__, __FILE__);
    } else {
      printf("# Prints the buffer at line %d in file %s:\n", __LINE__, __FILE__);
      for(uint i = 0; i < size; i++) {
	buffer[i].print(i);
      }
    }
  }
  //! Prints the values set for this object using the names given as argument.
  void print(char **name_out) {
    printf("\t%s %u %u \t%f (%d)\n",
	   name_out[ind_out], overlap_in, overlap_out, distance, (int)cnt_max);
  }

  /**
     @brief Reserves a list for 'this'
     @remarks 
     - Do not initiate the data for the objects. 
     - Consitutes the standard method in this class for reserving memory for several objects of this type.
     - Do not trhow an error if data was not allocted. (A bunch of methods using this one, ie, har to know where it was cast from.)
     @param <buffer_size> The number of objects to produce, i.e. the length of the list returned.
     @returns a list of objects of this class using the input as length of it
  **/
  static p_rel *reserve_list(uint buffer_size) {
    assert(buffer_size != UINT_MAX);
    if(buffer_size) {
      p_rel *ret = NULL; //new p_rel[buffer_size]();
      try {ret = new p_rel[buffer_size]();} 
      catch (std::exception& ba) {
	fprintf(stderr, "!!\t Tried reserving a \"p_rel\" class set of length=%u, [%s]:%s:%d\n",  buffer_size,  __FUNCTION__, __FILE__, __LINE__);
	if(!log_builder::catch_memory_exeception(buffer_size, __FUNCTION__, __FILE__, __LINE__, true)) {
	  fprintf(stderr, "!!\t An interesting error was discovered: %s."
		  "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		  "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
	}
      }
      return ret;
    }
    else return NULL;
  }

  /**
     @brief Reserves- and initiates a list of p_rel objects.
     @param <buffer_size> The number of objects to produce, i.e. the length of the list returned.
     @param <buffer_in_mem_pos> A strict way of ensuring that the psotion starts at the start of it.
     @returns a list of objects of this class using the input as length of it
  **/
  static p_rel *init_list( loint buffer_size, loint &buffer_in_mem_pos) {
    if(buffer_size > 0) {
      p_rel *buffer = NULL; //new p_rel[buffer_size]();
      try {buffer = new p_rel[buffer_size]();} 
      catch (std::exception& ba) {
	fprintf(stderr, "!!\t Tried reserving a \"p_rel\" class set of length=%llu, [%s]:%s:%d\n",  buffer_size,  __FUNCTION__, __FILE__, __LINE__);
	if(!log_builder::catch_memory_exeception(buffer_size, __FUNCTION__, __FILE__, __LINE__, true)) {
	  fprintf(stderr, "!!\t An interesting error was discovered: %s."
		  "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		  "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
	}
      }
      memset(buffer, 0, sizeof(p_rel)*buffer_size); // Prevents valgrind from complaining, as valgrind chekcs if its set (ie, initalized)
      for(loint i = 0; i<buffer_size; i++) buffer[i] = p_rel(); // This line states if the default- or OrthoMcl procedure is to be followed.
      buffer_in_mem_pos = 0; // A strict way of ensuring that the psotion starts at the start of it
      return buffer;
    } else return NULL;
  }  
  /**
     @brief Init a list of this class.  
     @param <buffer_size> The number of objects to produce, i.e. the length of the list returned.
     @param <buffer> The list contaning objects of this type (p_rel) to initialise.
  **/
  static void init(uint buffer_size, p_rel *&buffer) {
    init(0, buffer_size, buffer);
    //    for(uint i = 0; i<buffer_size; i++) buffer[i] = p_rel(_USE_BEST_BLAST_PAIR_SCORE);
  }
  /**
     @brief Init a list of this class.
     @param <start> The start position in the buffer to initialize.
     @param <end> The (last-1) position in the buffer to initialize.
     @param <buffer> The list contaning objects of this type (p_rel) to initialise.
     @remarks
     - Consitutes the standard method in this class for reserving memory for several objects of this type. 
  **/
  static void init(uint start, uint end, p_rel *&buffer) {
    for(;start<end; start++) buffer[start] = p_rel();
  }
  
  //! Deletes a list containing elements of this class.
  static void delete_list(p_rel *&buffer, loint &pos_in_buffer, loint &size_of_buffer ) {
    close(buffer, pos_in_buffer, size_of_buffer);
  }
  //! Deletes a list containing elements of this class.
  static void close(p_rel *&buffer, loint &pos_in_buffer, loint &size_of_buffer) {
    if(buffer && size_of_buffer) {delete [] buffer, buffer = NULL;}
    pos_in_buffer = 0, size_of_buffer = 0;
  } 
  //! Copies a list of these from second argument to first.
  static void copy(loint size_destination, p_rel *&destination, p_rel *source, loint size_source) {
    assert(size_destination >= size_source);
    assert(source != NULL);
    if(true)	memcpy(destination, source, sizeof(p_rel)*size_source);
    else {
      for(loint i = 0; i< size_source; i++) destination[i] = source[i];
    }
  }
  //! Copies a list of these from second argument to first at specific position.
  static void copy(loint size_destination, p_rel *&destination, uint start_pos_destination, p_rel *source, uint start_pos_source, uint size) {
    assert(size_destination >= (start_pos_source + size));
    assert(source != NULL);
    memcpy(destination+start_pos_destination, source+start_pos_source, sizeof(p_rel)*size); // Copies the data
  }

  /**
     @brief Enlarges the given list.
     @remarks If not set, allocates a new list.
     @param <buffer> The buffer to be extended at its lower end.
     @param <size_old> Used during DEBUG to assert the size is correct.
     @param <size_new> The length of the empty room to be made in front of the (new) extended list.
     @param <length_of_data_set> The number ob p_rel objects to be copied from the old buffer into the new one.
     @remarks The size_new parameter is not updated, ie, calling function must perform the update. (This message due to an difficult error made earlier, ie, of apriori knowledge.)
  **/
  static void enlarge(p_rel *&buffer, loint size_old, loint size_new, loint length_of_data_set) {
    if(size_new > size_old) { // If it's not, there is not point in enlarging it.
      assert(size_old <= size_new);
      assert(length_of_data_set <= size_new);
      if(buffer) {
	p_rel *buff_temp = reserve_list(size_new); 
	memset(buff_temp, 0, sizeof(p_rel)*size_new); // Prevents valgrind from complaining, as valgrind chekcs if its set (ie, initalized)
	memcpy(buff_temp, buffer, sizeof(p_rel)*length_of_data_set);    
	delete [] buffer; buffer = buff_temp; // Updates the pointer with the new data.
      } else {
	buffer = reserve_list(size_new);
	memset(buffer, 0, sizeof(p_rel)*size_new); // Prevents valgrind from complaining, as valgrind chekcs if its set (ie, initalized)
      }
      if(length_of_data_set) {
	for(loint i = length_of_data_set; i < size_new; i++) {buffer[i]=p_rel();}
      }
    }

  }
  /**
     @brief Enlarges the given list an puts the 'current' at the back of it (reciprocal to what's normal).
     @remarks This procedure makes empty room in front of it, in order to enable copying of data from a file into it.
     @param <buffer> The buffer to be extended at its lower end.
     @param <new_size> Used during DEBUG to assert the size is correct.
     @param <start_pos> The length of the empty room to be made in front of the (new) extended list.
     @param <copy_length> The number ob p_rel objects to be copied from the old buffer into the new one.
     @param <_USE_BEST_BLAST_PAIR_SCORE>  If set, uses the best score found in the blast file, i.e. do not merges multiple scores for the same protein.
  **/
  static void enlarge_and_put_current_at_back(p_rel *&buffer, loint new_size, loint start_pos, loint copy_length, bool _USE_BEST_BLAST_PAIR_SCORE) {
    p_rel *buff_temp = reserve_list(new_size); //(p_rel*)malloc(sizeof(p_rel)*new_size);
    if(buff_temp) {
      assert(new_size <= (start_pos+copy_length)); // If procedure is used correct, this should hold. 
      for(loint i = 0; i < start_pos; i++) {buff_temp[i]=p_rel();}    
    // Inserts the buffer in the front of the file (i.e. in
      // order for the indexes to maintain their correct position pointers
      memcpy(buff_temp +start_pos, buffer, sizeof(p_rel)*copy_length);    
      delete [] buffer; 
      buffer = buff_temp;
    } else fprintf(stderr, "Not enough memory reserving %lld elements of type p_rel (having size %d). This message was generated at line %d in file%s, found in method %s. Please contact the developer at oekseth@gmail.com if this message is seen.\n", new_size, (int)sizeof(p_rel), __LINE__, __FILE__, __FUNCTION__);
  }
  

  /**
     @brief Writes the buffer given to a file:
     @param <buffer> The object to write. Do not change any of its data, i.e. the data is still kept in the buffer after the file writing.
     @param <file_name> The identifier later used reading the file now generated (or extended).
     @param <buffer_size> The number of relevant objects in the buffer (i.e. the 'size'-those allocated, but not used.)
     @param <file_length> The length of the current file.
  **/
  static void write_buffer(p_rel *&buffer, char *file_name, loint buffer_size, loint file_length) {
    assert(file_name);
    if(buffer!= NULL && (buffer_size > 0)) {
      FILE *file = NULL;
      if(file_length == 0) file = fopen(file_name, "wb");
      else file = fopen(file_name, "ab");
      if(file != NULL) {
	printf("writes %lld p_rel elements to file %s at line %d in source-file %s\n", buffer_size, file_name, __LINE__, __FILE__);
	fwrite(buffer, sizeof(p_rel), buffer_size, file);
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
  static void write_buffers(char *file_name, p_rel *buffer_one, uint buffer_one_size, p_rel *buffer_two, uint buffer_two_size, uint file_length) {    
    FILE *file = NULL;
    assert(file_name);
    if(file_length == 0) file = fopen(file_name, "wb");
    else file = fopen(file_name, "ab");
    if(file != NULL) {
      if(buffer_one && buffer_one_size > 0) {
	fwrite(buffer_one, sizeof(p_rel), buffer_one_size, file);
      }
      if(buffer_two && buffer_two_size > 0) {
	fwrite(buffer_two, sizeof(p_rel), buffer_two_size, file);
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
  static p_rel *read_file(char *file_name, p_rel *buffer, loint size) {
    if(buffer != NULL) {
      FILE *file = fopen(file_name, "rb");
      if (file != NULL) {
	const loint size_read = fread(buffer, sizeof(p_rel), size, file); // to avoid the intial char
	if(size_read < size) init(size_read, size, buffer);
	assert(size_read <= size);
	fclose(file);
	return buffer;
      } else {
	fprintf(stderr, "!!\tCould not open '%s': input given has buffer_set=%d, p_rel_count(%lld) and total_size %d B; Aborts at line %d in %s\n", file_name, (int)(buffer!=NULL), size, (int)(sizeof(p_rel)*size), __LINE__, __FILE__);
	exit(2);
	return NULL;
      } 
    } else {
      fprintf(stderr, "!!\tCould not reserve %u amounts of class 'p_rel'; Aborts at line %d in %s!\n", (uint)size, __LINE__, __FILE__);
      exit(2);
    }
    return buffer;
  }

  /**
     @brief Reads a bunch of p_rel objects from the file given.
     @param <file_name> The name of the file to get the data from.
     @param <start_pos> The index in the file to read the first object from.
     @param <length> The number of objects to read.
     @returns A list containing the objects read.
  **/
  static p_rel* get_file_block(char *file_name, uint start_pos, uint length) {
    FILE *file = fopen(file_name, "rb");
    if (file != NULL) {
      //      p_rel *buffer_ret =  tbb::tbb_allocator<p_rel>().allocate(length);
      p_rel *buffer_ret =  NULL; //       tbb::tbb_allocator<rel>().allocate(length);
      try {buffer_ret = tbb::tbb_allocator<p_rel>().allocate(length);} 
      catch (std::exception& ba) {
	if(!log_builder::catch_memory_exeception(length, __FUNCTION__, __FILE__, __LINE__)) {
	  fprintf(stderr, "!!\t An interesting error was discovered: %s."
		  "The tool will therefore crash, though if you update the developer at [oekseth@gmail.com]."
		  "Error generated at [%s]:%s:%d\n",  ba.what(), __FUNCTION__, __FILE__, __LINE__);
	}
      }

      if(buffer_ret != NULL) {
	fseek(file, sizeof(p_rel)*start_pos, SEEK_SET);
	memset(buffer_ret, 0, sizeof(p_rel)*length);
	const uint read_length = fread(buffer_ret, sizeof(p_rel), length, file); // to avoid the intial char
	assert(read_length <= length);
	fclose(file);
	//	for(uint i = 0; i <length; i++) buffer_ret[i].print();
	return buffer_ret;
      } else { fprintf(stderr, "!!\tCould not allocate %u types of class p_rel; Aborts at line %d in %s!\n", length, __LINE__, __FILE__);
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
  static void free_file_block(p_rel *&block, uint size) {
    {//						if(block =! NULL) {
      tbb::tbb_allocator<p_rel>().deallocate((p_rel*)block, size);
      block = NULL;
    }
  }

  //! Reallocates the list, but make the room available in the start, instead of the latter as normal.
  static void realloc_making_free_space_in_front_of_list(p_rel *&buffer, loint old_size, loint new_size,  loint &buffer_in_mem_pos) {
    if(!(old_size > new_size)) {
      assert(!(old_size > new_size));
    }    
    if(buffer) {
      assert(old_size != 0); // If it's correct that the buffer has data, this must also have value
      {
	const loint adjustment_length = new_size-old_size;
	p_rel *buff_temp = reserve_list(new_size); //(p_rel*)malloc(sizeof(p_rel)*new_size);
	if(buff_temp) {
	  //	  for(loint i = 0; i < adjustment_length; i++) {buff_temp[i]=p_rel();}    
	  for(loint i = 0; i < new_size; i++) {buff_temp[i]=p_rel();}    
	  assert((adjustment_length+old_size) <= new_size);
	  memcpy(buff_temp+adjustment_length, buffer, (sizeof(p_rel)*old_size));
	  delete [] buffer;     buffer = buff_temp;
	  buffer_in_mem_pos = adjustment_length; // TODO: verify this!
	} else fprintf(stderr, "Not enough memory reserving %lld elements of type p_rel (having size %d). This message was generated at line %d in file%s, found in method %s. Please contact the developer at oekseth@gmail.com if this message is seen.\n", new_size, (int)sizeof(p_rel), __LINE__, __FILE__, __FUNCTION__);
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
  static void realloc_list(p_rel *&buffer, loint old_size, loint new_size, loint &buffer_in_mem_pos) {
    if(!(old_size > new_size)) {
      assert(!(old_size > new_size));
    }    
    if(buffer) {
      assert(old_size != 0); // If it's correct that the buffer has data, this must also have value
      p_rel::enlarge(buffer, old_size, new_size, /*length_of_data_set=*/buffer_in_mem_pos);
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
  static bool has_data(uint search_index, p_rel *buffer, uint start, uint end) {
    assert(start < end);
    for(uint k = start; k< end; k++) {
      if(buffer[k].is_equal(search_index)) {
	return true;
      }
    }
    return false;
  }

  //! The constructor.
  p_rel(uint _out, float _dist, short int _overlap_in, short int _overlap_out, unsigned char _cnt_max)   :
   /*ind_out(_out), distance(_dist),*/ overlap_in(_overlap_in), overlap_out(_overlap_out),
   cnt_max(_cnt_max), ind_out(_out), distance(_dist) {}

  //! The constructor.
 p_rel() : overlap_in(0), overlap_out(0),
   cnt_max(), ind_out(0), distance(0.0) {}


  //! The assert method for this class
  static void assert_class(bool print_info) {
    const static char *class_name = "p_rel_t";
    if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
    p_rel test = p_rel();
    const bool use_best_blast_score = false;
    // Asserts when proteins equals:
    test.set_data(0, 0.0, 2);
    assert(test.get_overlap_in() == 2);
    assert(test.get_distance() == 0.0);
    assert(test.get_cnt_max() == 1);
    test.increment_data(0.0, 2, use_best_blast_score);
    assert(test.get_overlap_in() == 4);
    assert(test.get_distance() == 0.0);
    assert(test.get_cnt_max() == 2);
    assert(test.getSimScore(10.0) == (2*(10+1)));

    // Asserts when proteins differs:
    test.set_data(0, 2.0, 2);
    assert(test.get_overlap_in() == 2);
    assert(test.get_distance() == 2.0);
    assert(test.get_cnt_max() == 0);
    test.increment_data(2.0, 2, use_best_blast_score);
    assert(test.get_overlap_in() == 4);
    assert(test.get_distance() == 4.0);
    assert(test.get_cnt_max() == 0);
    assert(test.getSimScore(1.0) == 4.0); // Sim-score should be the same.
#ifndef NDEBUG
    p_rel test_2 = p_rel((uint)2, (float)3, 4, 5, 6);
    p_rel test_3 = p_rel((uint)2, (float)3, 4, 5, 6);
    assert(test_2.is_equal(test_3));
#endif
#endif
    if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
  }
};

/**
   @brief Like class 'rel', but with some extension for storing temporary blast-data.
   @ingroup blastfile_container
   @author Ole Kristian Ekseth (oekseth)
**/
typedef class p_rel p_rel_t;
#endif
