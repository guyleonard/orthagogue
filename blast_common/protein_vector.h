#ifndef protein_vector_h
#define protein_vector_h
#include "../configure.h"
 
#include "protein_relation.h"
#include "prot_list.h"

/**
   @file
   @class protein_vector
   @brief Holds the data for the whole bunch of a taxon during parsing.
   @ingroup parsing_ops
   @remark Object takes a block of text formatted in the blast-standard: Goal
   of optimizing the retreival of the indexes in the hash build for the
   scoping (of the variables) required duiring paralelisation on a
   single core (multi cpu environment).
   @author Ole Kristian Ekseth (oekseth)
   @date 18.03.2011 by oekseth (initial)
   @date 16.09.2011 by oekseth (asserts)  
**/
class protein_vector {
 private:
/*   char string_taxon_in_prev[SIZE_TAXO]; */
/*   char string_taxon_out_prev[SIZE_TAXO]; */
/*   char string_name_in_prev[SIZE_PROTEIN]; // Used as reference in order to know if index is locally known */
/*   char string_name_out_prev[SIZE_PROTEIN]; // Used as reference in order to know if index is locally known */
  char *string_taxon_in_prev;
  char *string_taxon_out_prev;
  char *string_name_in_prev; // Used as reference in order to know if index is locally known
  char *string_name_out_prev; // Used as reference in order to know if index is locally known
  loint string_taxon_in_prev_size;
  loint string_taxon_out_prev_size;
  loint string_name_in_prev_size; // Used as reference in order to know if index is locally known
  loint string_name_out_prev_size; // Used as reference in order to know if index is locally known

  mem_loc this_index_in_prev; // the left protein
  mem_loc this_index_out_prev; // the right protein
  int this_taxon_in_prev;
  int this_taxon_out_prev;
  prot_list_t *hashProtein;
  prot_list_t *hashTaxa;
  int taxon_length;
  bool is_initiated; // set to true when the first protein of the input is met
//! Gets the protein index for the taxon specified
bool getProteinIndex(int taxon_id, char *protein, mem_loc &protein_ind) {
  if(taxon_id < taxon_length) {
    if(hashProtein) {
      if(hashProtein[taxon_id].getIndex(protein, protein_ind)) return true;
      else {
	protein_ind = INT_MAX;
	return false;
      }
    } else {
      fprintf(stderr, "!!\tError: 'hashProtein' at line %d in file %s found in method %s is not set. If this error is seen, please contact oekseth@gmail.com\n",
	      __LINE__, __FILE__, __FUNCTION__); return false;
    }
  } else {
    protein_ind = INT_MAX;
    return false;
  }
}

 public:
  //! The lsit holding the overlap data.
  overlap_t **arrOverlap;
  //! Inserted: If in debug-mode, this variable corresponds to the total sum of inserted ovelap values.
  uint debug_sum_overlap;  
  //! Ensures that only one "global overlap" is inserte:
  uint overlaps_inserted_by_curret_call;  
 private:
  //! 'has_data' set to true after first insert of values, used ensuring correct behaviour when updating overlap values.
  bool has_data;
  //! Uses the 'prot_list_t hashProtein' to retrive the index for the taxon
  bool getTaxonIndex(char *taxon, int &taxon_ind/*, const bool debug*/);
  //! Gets the protein index for the taxon specified
  bool getProteinIndex(int taxon_id, char *protein, int &protein_ind);

  /**
     Sets the taxon and returns the id of it
     @Return: false if the relation is not consisten (in accordance with the rules defined in the assumptions of the blast file)
  */
  bool set_taxon(Parse p, const bool is_inner, const bool is_equal);
  //! Sets the overlap for the given protein in the list of overlaps
  void setOverlap(int this_taxon_in_prev, mem_loc this_index_in_prev, overlap_t overlap);

  //! Copies protein name from 'src' to 'dest', updating the 'dest_size' variable.
  void copy_name(char *&dest, loint &size_dest, char *src, const loint size_src);
  //!  @return True if the strings given are equal.
  bool compare_strings(char *name_1, char *name_2);
  //!  @return True if the strings given are equal.
  bool compare_strings(char *name_1, char *name_2, uint length_string);

 public:
  //! Sets the previous pair value on the argument reference given.
  void get_prev_pair(Parse &p);
  //! Frees the memory for 'arrOverlap'
  void free_overlap();
  /**
     @brief Sets the protein and returns the id of it
     @return false if the relation is not consisten (in accordance with the rules defined in the assumptions of the blast file)
  **/
  bool set_protein(Parse p, const bool is_inner, const bool is_outer, const bool is_equal, bool &overlap_is_inserted);
  //! Returns a list of objects using 'this' type.
  static protein_vector* init(uint size);
  /**
     @brief Deallocates memory for this object.
     @remarks Does not included the hash list in the deallocation procedure (provided as a 'global external' argument for this object).
  **/
  void free_memory();  
  //! Deallocates memory of the object given.
  static void free_list(protein_vector *&obj, uint length);
  //! Updates the setting for the inner (left)/ and outer (right)  protein:
  struct protein_relation get_protein_indexes(Parse p, bool &overlap_is_inserted, taxa *listTaxa);
  //! The constructor; does not reserve any 'extra' memory.
  protein_vector();
  //! The constructor.
  protein_vector(prot_list_t *_hashProtein, prot_list_t *hashTaxa, int _taxon_length);

 public:
#ifdef assert_code
  void assert_private_parts();
#endif
  //! The main test function for this class  
  static void assert_class(const bool print_info) {
    const static char *class_name = "protein_vector";
    if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
    //uint taxon_length = 2, n_threads = 2;
    //protein_vector temp = protein_vector();
    //    temp.assert_private_parts();
    // temp.free_memory();
#endif
    if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
  }    


};

/**
   @brief Holds the data for the whole bunch of a taxon during parsing.
   @ingroup parsing_ops
   @author Ole Kristian Ekseth (oekseth)
*/
typedef class protein_vector protein_vector_t;

#endif
