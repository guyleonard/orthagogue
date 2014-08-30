/**
   @file
   @brief A class accessing- and bulding the perfect hash.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
**/
#ifndef prot_list_h
#define prot_list_h
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

#include "parse.h"
#include "types.h"
#include "cmph.h"
#include "taxa.h"

#include "log_builder.h"
/**
   @class prot_list
   @brief A class accessing- and bulding the perfect hash.
   @ingroup parsing_ops
   @remark Used for name retrieval during second pass of the blast file:
   The function 'getIndex(...)' is the most important function of this class, 
   returning the index for each name requested. Tests focused on the hashing.
   @date 21.12.2010 by oekseth (initial)
   @date 26.08.2011 by oekseth (asserts)
   @author Ole Kristian Ekseth (oekseth)
**/
class prot_list {
  // TODO: Seperate out in own cxx- file when located in own folder due to those cmkae-things...
 private:
  //! All the variables (and functions) will be thread safe (and thereby is accessible for multiple threadsat once
  bool variables_transfered_to_taxa; // If set, limits the varaibles being able to be deleted.
  mem_loc prots_length;
  char *taxon_name; // The name of this taxon
  cmph_t *hash;
  bool USE_RECIPROCAL_TEST;
  bool use_hash; // if the number of keys are to low, a hash can not be generated, and thereby a will instead be retrived from a list
  uint *back_track_index;  /*back_track_indexes set if a bakctracking index should be returned.*/
  cmph_io_adapter_t *source;
 public:
  //! Holds the index->name-relations for proteins of this taxon
  char **arrKey; 
  //! If set, limits the varaibles being able to be deleted.
  void lock_variables_transfered_to_taxa_list() {variables_transfered_to_taxa = true;}
  //! @return the back tracking indexes
  uint get_backtracking_index(int index) {
    assert(index < prots_length);
    if(back_track_index) return back_track_index[index];
    else return 0;
  }
  /**
     Initiates a list of prot_list objects of length given by param.
     @remarks Only allocates, but does not initialize, so use with care.
  **/
  static prot_list* init(uint size) {return (prot_list*)malloc(sizeof(prot_list)*size);}
  //! Deallocates the object given as param.
  static void free_list(prot_list *&obj) {free(obj); obj= NULL;}
 
  //! Returns the name of the taxon:
  char *getName() {return taxon_name;}

  //! Returns the length of the proteins:
  mem_loc getLength() {return prots_length;}

  //! Executes the building of the hash'es: Called from the constructor
  void buildHash(char **vector) {
    log_builder::test_memory_condition_and_if_not_abort(vector!=NULL, __LINE__, __FILE__, __FUNCTION__);
    log_builder::test_memory_condition_and_if_not_abort(prots_length>0, __LINE__, __FILE__, __FUNCTION__);
    const mem_loc nkeys = prots_length; // The number of hash-keys to produce
    source = cmph_io_vector_adapter((char **)vector, nkeys);   // Source of keys                                                                                                    
    //Create minimal perfect hash function using the brz algorithm:
    cmph_config_t *config = cmph_config_new(source);
    cmph_config_set_algo(config, CMPH_CHD);
    assert(config);
    log_builder::test_memory_condition_and_if_not_abort(config!=NULL, __LINE__, __FILE__, __FUNCTION__);
    //cmph_config_set_algo(config, CMPH_BRZ); // An alternative hashing algo
    hash = cmph_new(config);
    if(hash==NULL) {
      fprintf(stderr, "!!\t Hash was not constructed, probably due to rececing a non-sorted Blast input file.\n"
	      "-\t The execution-time will therefore be severly hampered;"
	      "-\t We would be utmost thankful if you would forward this message to the developer,\n"
	      " \t either through orthAgogue's issue-page (at our home-page), or directly to the developer at [oekseth@gmail.som].\n"
	      "This message was printed at [%s]:%s:%d\n",
	      __FUNCTION__, __FILE__, __LINE__); 
      use_hash = false;
    }  else {
      //printf("(ok: constructed hash, at prot_list:%d.)\n", __LINE__);
      cmph_config_destroy(config), use_hash = true;
    }
  }

  //! Builds the 'overlap' and 'list of keys' vectors: Called from the constructor
  void fillDataStructure(char **buffer) {
    // The other arrays:
    arrKey = (char**)malloc(sizeof(char*)*prots_length);
    int i = 0;
    if(use_hash) {
      while(i < prots_length) { // Making a reciprocal array:
	const char *key = buffer[i];
	const mem_loc id = (mem_loc)cmph_search(hash, key, (cmph_uint32)strlen(key));
	if(id > -1 && key) {
	  const loint size = strlen(key); arrKey[id] = new char[size+1]; arrKey[id][size] = '\0';
	  strncpy(arrKey[id], key, size);
	  if(back_track_index) back_track_index[id] = (uint)i;
	  //	  printf("back_track_index[%d] = %d for %s in prot_list.h\n", id, i, key);
	  //	strncpy(arrKey[id], key, SIZE_PROTEIN-1);
	} else arrKey[i] = NULL;
	i++;
      }
    } else while(i < prots_length) {
      if(buffer[i]) {
	const loint size = strlen(buffer[i]); arrKey[i] = new char[size+1]; arrKey[i][size] = '\0';
	strncpy(arrKey[i], buffer[i], size);
	//	strncpy(arrKey[i], buffer[i], SIZE_PROTEIN-1);
      } else arrKey[i] = NULL;
      i++;
    }
  }
  
  //! Intiaizes the data structures
 prot_list(mem_loc _prots_length, char **buffer, char *_taxon_name, const bool _USE_RECIPROCAL_TEST) :
   variables_transfered_to_taxa(false),
  prots_length(_prots_length),
    USE_RECIPROCAL_TEST(_USE_RECIPROCAL_TEST) // if set to true, tests if the protein is a new one
   , back_track_index(NULL)
    {
      if(_taxon_name) {
	const loint size = strlen(_taxon_name);
	taxon_name = new char[size+1];
	taxon_name[size] = '\0'; // A trailing char.
	strncpy(taxon_name, _taxon_name, size); 
      } else taxon_name = NULL;
      assert(buffer);
      assert(prots_length);
      buildHash(buffer);   //Builds the hash
      fillDataStructure(buffer); // Fills the datastruture
    }

  //! Intiaizes the data structures
   prot_list(mem_loc _prots_length, char **buffer, char *_taxon_name, const bool _USE_RECIPROCAL_TEST, bool use_back_track_indexes) :
   variables_transfered_to_taxa(false),
  prots_length(_prots_length),
    USE_RECIPROCAL_TEST(_USE_RECIPROCAL_TEST) // if set to true, tests if the protein is a new one
   , back_track_index(NULL)
    {
      if(use_back_track_indexes) {
	back_track_index = new uint[prots_length];
	for(uint i = 0; i < prots_length; i++) back_track_index[i] = 0;
      }
      if(_taxon_name) {
	const loint size = strlen(_taxon_name);
	taxon_name = new char[size+1];
	taxon_name[size] = '\0'; // A trailing char.
	strncpy(taxon_name, _taxon_name, size); 
      } else taxon_name = NULL;

      buildHash(buffer);   //Builds the hash
      fillDataStructure(buffer); // Fills the datastruture
    }

  //! Returns true if taxa are equal
  bool is_equal(char *taxon) {
    if(taxon_name) {
      if(taxon) return (0==strncmp(taxon_name, taxon, SIZE_TAXO-1));
      else return false;
    } else if(!taxon) return true;
    return false;
  }
  
  //! @return the string corresponding to the index.
  const char *get_string(const int index) const {
    if(arrKey && (index < (int)prots_length) && arrKey[index]) {
      return arrKey[index];
    } else {return NULL;}
  }

  /**
     @brief Calculates the index of the given string.
     @param <key>  The name of the protein (witout the taxon identifier)
     @param <index> The index found in the hash table (and therefore the unique identifier inside this taxon)
     @param <debug> if set to true, then print debugging-information at special poins in the exeuction.
     @return True if reciprocal; the index given as a parameter is set
     @date 26.08.2011 by oekseth (asserts).
  */
  bool getIndex(char *key, int &index, const bool debug=false) {
    if(key != NULL) {
      if(*key == '\n') key++; // Sometimes a 'newline' is set at its start.
      if(use_hash) {
	index = (mem_loc)cmph_search(hash, key, (cmph_uint32)strlen(key));
	if(index < prots_length) {
	  if(USE_RECIPROCAL_TEST) { // Tests if reciprocal: practical if blast format has errors in it (found during the intial testing of this code)
	    if(arrKey[index]) { // is set:
	      if(strlen(arrKey[index]) != strlen(key)) return false;
	      else if(0==(strncmp(arrKey[index], key, strlen(arrKey[index])))) {return true;}
	      else return false; // the index found is not a correct one
	    } else { // The index was not found for the leftmost protein:
	      if(false) { // TODO: Consider including this error with modifications, remopving noise when filer is used.
		char temp[100]; memset(temp, '\0', 100);
		sprintf(temp, "Label(%s) not defined, therefore discards the line(s) containing it", key);
		log_builder::throw_warning(blastp_syntax, __LINE__, __FILE__, __FUNCTION__, temp);
	      }
	      return false;
	    }
	  } else return true;
	} else { // Not found: // TODO: below some possible ways of stating these. Consider the best way of expressing this.
	  //	fprintf(stderr, "!!\tIndex(%ld) for char (%s) !< prots_length(%ld). Error located in class 'prot_list'. Implies that blast file do not have a (%s--%s) relation. Probably caused by abnormalities in blast-file. Contact %s at %s for questions. Updates the taxon_id.\n", index, key, prots_length, key, key, DEVELOPER_NAME, DEVELOPER_EMAIL);
	  // TODO: If problem is reaccureing, consider a user-specific option ignoring these warnings.
	  /*	fprintf(stderr, "!!\tIndex(%ld) for char (%s) !< prots_length(%ld). Error located in class 'prot_list'. Implies that blast file do not have a (%s--%s) relation. Probably caused by abnormalities in blast-file. Contact %s at %s for questions.\n", index, key, prots_length, key, key, DEVELOPER_NAME, DEVELOPER_EMAIL);
		}
		printf("# '%s' not found (total of %d of type '%s'), As seen below, the above key is not in the lsit.\n", key, (int)prots_length, getName()); print_data(); 
	  */
	}
      } else {
	if(debug) {printf("\t get the key \"%s\" through the slow approach, as prot_list:%d\n", key, __LINE__);} // FIXME: remove.
	for(mem_loc i = 0;i<prots_length;i++) { // Iterates thorugh the list
	  if(0==(strncmp(arrKey[i], key, strlen(arrKey[i])))) {index=i; return true;}
	  //	  if(0==(strncmp(arrKey[i], key, SIZE_PROTEIN-1))) {index=i; return true;}
	}
      }
    }
    if(false) { // TODO: Consider including this error with modifications, remopving noise when filer is used.
      char temp[100]; memset(temp, '\0', 100);
      sprintf(temp, "Label(%s) not defined, therefore discards the line containing it", key);
      log_builder::throw_warning(blastp_syntax, __LINE__, __FILE__, __FUNCTION__, temp);
    }
    return false; // key (protein label) not found
  }

  /**
     @brief Calculates the index of the given string.
     @param <key>  The name of the protein (witout the taxon identifier)
     @param <index> The index found in the hash table (and therefore the unique identifier inside this taxon)
     @return True if reciprocal; the index given as a parameter is set
     @date 26.08.2011 by oekseth (asserts).
  */
  bool getIndex(char *key, mem_loc &index) {
    int ind = (int)index;
    const bool ret = getIndex(key, ind);
    index = (mem_loc)ind;
    return ret;
  }

  //! Allocates- and intitates the data for a taxon input:
  static prot_list *init(taxa_t *listTaxa, int taxon_length) {    
    char **buffer = taxa::get_taxa_names(listTaxa, taxon_length);
    prot_list *obj = new prot_list(taxon_length, buffer, NULL, true, /*back_track_indexes=*/true);
    for(int i = 0; i < taxon_length; i++) {delete [] buffer[i];} delete [] buffer; buffer = NULL;    
    return obj;
  }
  //! Deallocates the memory.
  void delete_buffer() {
    if(hash) {cmph_destroy(hash); hash = NULL;}
    if(source) {cmph_io_vector_adapter_destroy(source), source=NULL;}; 
    if(back_track_index) {delete [] back_track_index, back_track_index = NULL;}
    if(taxon_name) {delete [] taxon_name; taxon_name = NULL;}
    if(!variables_transfered_to_taxa) {
      for(loint i = 0; i < (loint)prots_length; i++) {
	delete [] arrKey[i]; arrKey[i] = NULL;
      }
      //    free(arrKey[0]);
      free(arrKey); arrKey = NULL;
    } else {
      arrKey = NULL;
      taxon_name = NULL;
    }
  }

  //! Deallocates the reserved memory for the class:
  static void close(prot_list *&obj) {
    if(obj) {
      obj->delete_buffer();
      delete obj, obj = NULL;
    }
  }

  //! Deallocates the reserved memory for the list of objects
  static void close(prot_list *&obj, int taxon_length) {
    if(obj) {
      for(int i = 0; i <taxon_length; i++) {
	obj[i].delete_buffer();
      }
      free(obj); obj = NULL;
      //      delete [] obj, obj = NULL;
    }
  }
};
/**
   @brief The prot list.
   @ingroup parsing_ops
   @date 26.12.2011 by oekseth (asserts)
   @author Ole Kristian Ekseth (oekseth)
**/
typedef prot_list prot_list_t;

#endif
