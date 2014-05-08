#include "protein_vector.h"
//! Copies protein name from 'src' to 'dest', updating the 'dest_size' variable.
void protein_vector::copy_name(char *&dest, loint &size_dest, char *src, const loint size_src) {
  if(src) {
    //! Sometimes a newline is found in the front of the string:
    char *tmp = strchr(src, '\n');
    char *end = src + size_src;
    loint cnt = 0;
    while(tmp && (src != end) && (cnt < size_src)) {
      if(tmp) src = tmp +1;
      cnt++;
    }
#ifndef NDEBUG
  //! Asserts that no newlines are found in the names:
  for(uint i = 0; i < (uint)strlen(src); i++) {
    assert(src[i] != '\n');
  }
#endif
    if(dest) {
      if(size_dest < size_src) {
	delete [] dest; dest = new char[size_src+1]; 
	dest[size_src] = '\0';
	memset(dest, '\0', size_src+1);
      }
    } else {dest = new char[size_src+1];
      //dest[size_src] = '\0';
      memset(dest, '\0', size_src+1);
    }
    strncpy(dest, src, size_src);
#ifndef NDEBUG
  //! Asserts that no newlines are found in the names:
  for(uint i = 0; i < (uint)strlen(dest); i++) {
    assert(dest[i] != '\n');
  }
#endif
  }
}


//!  @return True if the strings given are equal.
bool protein_vector::compare_strings(char *name_1, char *name_2, uint length_string) {
  if(name_1 && name_2) {
    return (0 == strncmp(name_1, name_2, length_string));
  } else { // If both are not set, they are equal.
    if((name_1 == NULL) && (name_2 == NULL)) return true;
    else return false; // only one of them is set ie, they are not equal.
  }
}

//!  @return True if the strings given are equal.
bool protein_vector::compare_strings(char *name_1, char *name_2) {
  if(name_1 && name_2) {
    return (0 == strcmp(name_1, name_2));
  } else { // If both are not set, they are equal.
    if((name_1 == NULL) && (name_2 == NULL)) return true;
    else return false; // only one of them is set ie, they are not equal.
  }
}

//! Frees the memory for 'arrOverlap'
void protein_vector::free_overlap() {
  for(int i = 0; i< taxon_length;i++) 
    free(arrOverlap[i]);
  free(arrOverlap); arrOverlap = NULL;
}

//! Uses the 'prot_list_t hashProtein' to retrive the index for the taxon
bool protein_vector::getTaxonIndex(char *taxon, int &taxon_ind/*, const bool debug*/) {
  assert(taxon); // If not set, the usage of this method is wrong.
  if(hashTaxa) { 
    int temp = -1;
    if(hashTaxa->getIndex(taxon, temp)) {
      taxon_ind = hashTaxa->get_backtracking_index(temp);
      return true;
    }
  } else {
    for(int taxon_id = 0;taxon_id<taxon_length;taxon_id++) {
      if(hashProtein[taxon_id].is_equal(taxon)) {
	taxon_ind = taxon_id;
// 	if(true) {
// 	  if(hashTaxa) {
// 	    int temp = -1;
// 	    hashTaxa->getIndex(taxon, temp);
// 	    const uint back_track = hashTaxa->get_backtracking_index(temp);
// 	    printf("%s\treal_taxon(%d) and hash(%d --> %u) in protein_vector\n", taxon, taxon_ind, temp, back_track);
// 	  }
// 	}
	return true;
      }
    }
  }
  if(false) { // TODO: Consider including this error with modifications, remopving noise when filer is used.
    char temp[100]; memset(temp, '\0', 100);
    sprintf(temp, "Taxonlabel(%s) not defined, therefore discards the line containing it", taxon);
    log_builder::throw_warning(blastp_syntax, __LINE__, __FILE__, __FUNCTION__, temp);
  }
  taxon_ind = INT_MAX; 
  return false;
}

//! Gets the protein index for the taxon specified
bool protein_vector::getProteinIndex(int taxon_id, char *protein, int &protein_ind) {
  if(taxon_id < taxon_length) {
    if(hashProtein) {
      if(hashProtein[taxon_id].getIndex(protein, protein_ind)) return true;
      else {
	protein_ind = INT_MAX;
	return false;
      }
    } else {
      if(false) { // TODO: Consider including this error with modifications, remopving noise when filer is used.
	char temp[100]; memset(temp, '\0', 100);
	sprintf(temp, "Proteinlabel(%s) not defined, therefore discards the line containing it", protein);
	log_builder::throw_warning(blastp_syntax, __LINE__, __FILE__, __FUNCTION__, temp);
      }
      protein_ind = INT_MAX;
      return false;
    }
  } else {
    protein_ind = INT_MAX;
    return false;
  }
}


/**
   Sets the taxon and returns the id of it
   @Return: false if the relation is not consisten (in accordance with the rules defined in the assumptions of the blast file)
*/
bool protein_vector::set_taxon(Parse p, const bool is_inner, const bool is_equal) {
  bool format_is_correct = true; 
  //  const bool debug = false;
  if(is_equal) { // both has changed
    if(!getTaxonIndex(p.get_taxon_in(), this_taxon_in_prev)) format_is_correct=false;
    this_taxon_out_prev  = this_taxon_in_prev;
    if(string_taxon_out_prev) memset(string_taxon_out_prev, '\0', string_taxon_out_prev_size); // Do not set name due to the fact that an equal hit will only occure once
    // Updates the name:
    copy_name(string_taxon_in_prev, string_taxon_in_prev_size, p.get_taxon_in(), p.get_taxon_in_length());
  } else if(is_inner) { // left has changed; Updates the name.
    copy_name(string_taxon_in_prev, string_taxon_in_prev_size, p.get_taxon_in(), p.get_taxon_in_length());
    if(!getTaxonIndex(p.get_taxon_in(), this_taxon_in_prev)) format_is_correct=false;
  } else { // right has changed: string_taxon_out_prev --> p.taxo_out;  Updates the name.
    copy_name(string_taxon_out_prev, string_taxon_out_prev_size, p.get_taxon_out(), p.get_taxon_out_length());
    if(!getTaxonIndex(p.get_taxon_out(), this_taxon_out_prev)) {format_is_correct=false;}
  }
  return format_is_correct;
}

//! Sets the overlap for the given protein in the list of overlaps
void protein_vector::setOverlap(int this_taxon_in_prev, mem_loc this_index_in_prev, overlap_t overlap) {
  arrOverlap[this_taxon_in_prev][this_index_in_prev] += overlap;
#ifndef NDEBUG
  debug_sum_overlap += (uint)overlap;
#endif
  //assert(overlaps_inserted_by_curret_call == 0);
  overlaps_inserted_by_curret_call++;
  

}
/**
   Sets the protein and returns the id of it
   @REQURES: Executes after the update of the taxon
   @Return: false if the relation is not consisten (in accordance with the rules defined in the assumptions of the blast file)
*/
bool protein_vector::set_protein(Parse p, const bool is_inner, const bool is_outer, const bool is_equal, bool &overlap_is_inserted) {
  bool format_is_correct = true;
  #ifndef NDEBUG
  if(is_equal) {
    //! Asserts the proteins are equal:
    char *name_in = p.get_name_in();
    char *name_out = p.get_name_out();
    assert(0 == strcmp(name_in, name_out));
  } 
  #endif

  if(is_equal) { // both has changed; updates:
    if(!(getProteinIndex(this_taxon_in_prev, p.get_name_in(), this_index_in_prev))) format_is_correct=false;
    else {
      setOverlap(this_taxon_in_prev, this_index_in_prev, p.get_overlap_in());
      overlap_is_inserted = true;
      has_data = true;
    }
    this_index_out_prev = this_index_in_prev; // equal
    // Updates the name:
    copy_name(string_name_in_prev, string_name_in_prev_size, p.get_name_in(), p.get_name_in_length());
    copy_name(string_name_out_prev, string_name_out_prev_size, p.get_name_out(), p.get_name_out_length());
  } else {
    if(is_inner) { // left has changed
      if(!(getProteinIndex(this_taxon_in_prev, p.get_name_in(), this_index_in_prev))) format_is_correct=false;
      // Updates the name:
      copy_name(string_name_in_prev, string_name_in_prev_size, p.get_name_in(), p.get_name_in_length());
    } 
    if(is_outer) { // right has changed
      if(!(getProteinIndex(this_taxon_out_prev, p.get_name_out(), this_index_out_prev))) format_is_correct=false;
      // Updates the name:
      copy_name(string_name_out_prev, string_name_out_prev_size, p.get_name_out(), p.get_name_out_length());
    }
  }
  return format_is_correct;
}

//! Note: Sometimes a newline is found in the front of the string:
char *get_string_with_newlines_cleared(char *string) {
  char *src = string;
  char *tmp = strchr(src, '\n');
  const loint size_src = strlen(src);
  char *end = src + strlen(src);
  loint cnt = 0;
  while(tmp && (src != end) && (cnt < size_src)) {
    if(tmp) src = tmp +1;
    cnt++;
  }
  return src;
}
//! Updates the setting for the inner (left)/ and outer (right)  protein:
struct protein_relation protein_vector::get_protein_indexes(Parse p, bool &overlap_is_inserted, taxa *listTaxa) {
  struct protein_relation data = protein_relation();
  // Updates the setting for the inner (left) protein:
  overlap_is_inserted = false;
  //! Ensures that only one "global overlap" is inserte:
  overlaps_inserted_by_curret_call = 0;

  //! Our summary of the parsing:
  enum {id_taxon_in, id_taxon_out, id_name_in, id_name_out} field_id;
  if(false && field_id) ; // To hide the indirect usage of the above enum.
  bool field_has_changed[4] = {false, false, false, false};
  bool input_is_valid[4] = {true, true, true, true};

#ifndef NDEBUG
      char *name_in = p.get_name_in(); char *name_out = p.get_name_out();
      //! Asserts that no newlines are found in the names:
      for(uint i = 0; i < (uint)strlen(name_in); i++) {
	assert(name_in[i] != '\n');
      }
      //! Asserts that no newlines are found in the names:
      for(uint i = 0; i < (uint)strlen(name_out); i++) {
	assert(name_out[i] != '\n');
      }
#endif
  //! Note: If not initiated, finds the index for both the inner- and the outer protein.
  //! Proteins: Tests if the names has changed:
  if(!is_initiated || !(compare_strings(p.get_name_in(), string_name_in_prev)))    field_has_changed[id_name_in]   = true;
  if(!is_initiated || !(compare_strings(p.get_name_out(), string_name_out_prev)))  field_has_changed[id_name_out]  = true;
  //! taxa: Tests if the names has changed: &string_taxon_in_prev[0]--> p.taxon_in || &string_taxon_out_prev[0]--> p.taxon_out
  if(!is_initiated ||  !(compare_strings(p.get_taxon_in(), string_taxon_in_prev)))   field_has_changed[id_taxon_in]  = true; // Taxon inner changed: 
  if(!is_initiated ||  !(compare_strings(p.get_taxon_out(), string_taxon_out_prev))) {
    field_has_changed[id_taxon_out] = true; // right has changed
  } 
  if(false) {
    // TODO: Below code is included if further analyse of the code is required.
    char *str[2] = {"no", "yes"};
    printf("changed:\t taxon_in_changed(%s), taxon_out_changed(%s)[%s->%s], label_in_changed(%s)[%s->%s], label_out_changed(%s)[%s->%s], is_initiated(%s), at protein_vector:%d\n",
	   str[field_has_changed[id_taxon_in]],
	   str[field_has_changed[id_taxon_out]],
	   string_taxon_out_prev, p.get_taxon_out(), 
	   str[field_has_changed[id_name_in]],
	   string_name_in_prev, p.get_name_in(), 
	   str[field_has_changed[id_name_out]],
	   string_name_out_prev, p.get_name_out(), 
	   str[is_initiated],
	   __LINE__
	   ); // TODO: remove this printf!
  }
  //! If the inner names has changed, 
  if(field_has_changed[id_name_in]) { // Not equal; Assumes that the names are unique
    // The inner protein has changed: tests if the taxon has changed:
    if(field_has_changed[id_taxon_in]) {
      if(!set_taxon(p, true, true)) {
	data.set_as_no_exsists(); 
	input_is_valid[id_taxon_in] = false;
      }
    } else {
      if(this_taxon_out_prev == INT_MAX) {data.set_as_no_exsists();} 
      else {
	if(0 == strcmp(p.get_name_in(), p.get_taxon_out())) {
	  this_taxon_out_prev = this_taxon_in_prev; // Protein in changed from string_name_out_prev --> p.get_name_out());
	  field_has_changed[id_taxon_out] = true;
	} else {
	  //! Then it is a special case where the protein does not have a self-comparison as its first case.
	  //! This could be due to parallell exeuction when constructing the input Blast file.
	  getTaxonIndex(p.get_taxon_out(), this_taxon_out_prev);
	  field_has_changed[id_taxon_out] = true;
	}
	int temp_index = 0;  getTaxonIndex(p.get_taxon_out(), temp_index);
	assert(temp_index == this_taxon_out_prev);
      }
    }
  } else if(this_index_in_prev == INT_MAX) data.set_as_no_exsists();  
  // Updates the setting for the outer (right) protein:
  // NOTE: At frist look at the blastp-file it seems like it follows the label-id-format (AA->AA, AA->BB, AA->CC, DD->BB), but sometimes (AA->BB, AA->AA, AA->CC, DD->BB), therefore the last label must be checked independently of the first label.    
  if(input_is_valid[id_taxon_in] && field_has_changed[id_name_out]) { // Assumes that the names are unique
    // Protein out changed from string_name_out_prev --> p.get_name_out()
    if(field_has_changed[id_name_out]) { // Not equal; Assumes that the names are unique
      if(!set_taxon(p, false, false)) { // only the outer taxon has changed
	data.set_as_no_exsists();
	input_is_valid[id_taxon_out] = false;
      }
    } else if(this_taxon_out_prev == INT_MAX) data.set_as_no_exsists();
  } else {if(this_index_out_prev == INT_MAX) data.set_as_no_exsists();}

  //! Updates the data for the inner relation:
  //! Note: A bunch of taxa may not be valid, as the user may choose to ignroe txa having less than a threshold-size of proteins.
  if(data.exsists() && input_is_valid[id_taxon_out] && input_is_valid[id_taxon_in]) {
    //! Tests if the labels are equal
    bool is_equal = false;
    char *name_in = p.get_name_in(); char *name_out = p.get_name_out();
    if(0 == strcmp(name_in, name_out)) is_equal = true;

    //    bool format_is_correct = true;
    if(!set_protein(p, field_has_changed[id_name_in], field_has_changed[id_name_out], is_equal, overlap_is_inserted)){ // sets the protein
      data.set_as_no_exsists();
      //      format_is_correct = false;
      // input_is_valid[id_name_in] = false;
#ifndef NDEBUG
      //! Note: This error is hidden, as it tells about the files property, and not about any errros in the source:
      if(false) { 
	char temp[100]; memset(temp, '\0', 100);
	sprintf(temp, "Label(%s) for taxon(%s) not defined for leftmost relation, therefore discards the line containing it", p.get_name_in(), p.get_taxon_in());
	log_builder::throw_warning(blastp_syntax, __LINE__, __FILE__, __FUNCTION__, temp);
      }
#endif
    }
    is_initiated = true; // Update the variable for the next run:

    if(this_taxon_out_prev == INT_MAX) data.set_as_no_exsists();
    data.protein_in = this_index_in_prev;
    data.taxon_in = this_taxon_in_prev;
    data.protein_out = this_index_out_prev;
    data.taxon_out = this_taxon_out_prev;
#ifndef NDEBUG
    if(this_taxon_out_prev != INT_MAX) {
      assert(hashProtein[this_taxon_out_prev].is_equal(p.get_taxon_out()));
    }

   if(data.exsists()) {
      //! Asserts that the names has been updated for the next run:
      if(string_taxon_in_prev) {
	char *string = get_string_with_newlines_cleared(p.get_taxon_in());
	if(!((0 == strcmp(string, string_taxon_in_prev)))) {
	  printf("!!\tstring(%s)!= previous_string(%s), at %s:%d\n", p.get_taxon_in(), string_taxon_in_prev, __FILE__, __LINE__);
	  assert(0 == strcmp(p.get_taxon_in(), string_taxon_in_prev));
	}
      }
      if(string_name_out_prev) assert(0 == strcmp(p.get_taxon_out(), string_taxon_out_prev));
      if(string_taxon_in_prev) assert(0 == strcmp(p.get_name_in(), string_name_in_prev));
      if(string_name_out_prev) assert(0 == strcmp(p.get_name_out(), string_name_out_prev));
      
      //! Perform an assertion of the mapping
      assert(listTaxa);
      assert(0 == strcmp(p.get_name_in(), listTaxa[data.taxon_in].getArrKey(data.protein_in)));
      assert(data.protein_out  != INT_MAX);
      const bool out_is_equal = (0 == strcmp(p.get_name_out(), listTaxa[data.taxon_out].getArrKey(data.protein_out))); 
      if(out_is_equal == false) {
	char *str[2] = {"no", "yes"};
	fprintf(stderr, 
		"!!\t An error in identification of the rightmost keys, where we expected the vertex \"%s\"==\"%s\":\n"
		"-\t Evaluate a blast-row covering pair (%s, %s)\n"
		"-\t The error may be due to the format of the Blast data, wrong estimation of the protein-hash-key in question,\n"
		" \t or errors in our usage of the CMPH hash map library.\n" 
		"\t To help this investigation, we inspect the type of changed keys:\n"
		"\t- \t taxon_in_changed(%s),         taxon_out_changed(%s)[%s->%s],\n"
		"\t- \t label_in_changed(%s)[%s->%s]," " label_out_changed(%s)[%s->%s],\n" 
		"\t- \t w.r.t. the asserted validness of input ot the taxa, taxa-left=%s while taxa-right=%s\n"
		//"is_initiated(%s)\n",
		 "-\t We would be utmost thankful if you would forward this message to the developer,\n"
		 " \t either through orthAgogue's issue-page (at our home-page), or directly to the developer at [oekseth@gmail.som].\n"
		"This message was printed at [%s]:%s:%d\n",
		p.get_name_out(), listTaxa[data.taxon_out].getArrKey(data.protein_out),
		p.get_name_in(), p.get_name_out(), 
		str[field_has_changed[id_taxon_in]], str[field_has_changed[id_taxon_out]], string_taxon_out_prev, p.get_taxon_out(),
		str[field_has_changed[id_name_in]],   string_name_in_prev, p.get_name_in(), 
		str[field_has_changed[id_name_out]] ,  string_name_out_prev, p.get_name_out(),
		str[input_is_valid[id_taxon_in]], str[input_is_valid[id_taxon_out]], 
		__FUNCTION__, __FILE__, __LINE__); 
	assert(0 == strcmp(p.get_name_out(), listTaxa[data.taxon_out].getArrKey(data.protein_out)));       
      }
    }
#endif
  } else if(data.exsists()) {data.set_as_no_exsists();}

  return data;
}


//! Sets the previous pair value on the argument reference given.
void protein_vector::get_prev_pair(Parse &p) {
  if(string_name_in_prev)   p.set_name_in(string_name_in_prev, string_name_in_prev_size);
  if(string_name_out_prev)  p.set_name_out(string_name_out_prev, string_name_out_prev_size);
  if(string_taxon_in_prev)  p.set_taxon_in(string_taxon_in_prev, string_taxon_in_prev_size);
  if(string_taxon_out_prev) p.set_taxon_out(string_taxon_out_prev, string_taxon_out_prev_size);
}

//! Returns a list of objects using 'this' type.
protein_vector* protein_vector::init(uint size) {
  protein_vector *obj = (protein_vector*)malloc(sizeof(protein_vector)*size);
  assert(obj); // Cast an error if memory not allocated
  for(uint i = 0; i < size; i++) obj[i] = protein_vector();
  return obj;
}

/**
   @brief Deallocates memory for this object.
   @remarks Does not included the hash list in the deallocation procedure (provided as a 'global external' argument for this object).
**/
void protein_vector::free_memory() {
  if(string_name_in_prev)  {   delete [] string_name_in_prev,   string_name_in_prev = NULL,   string_name_in_prev_size = 0;}
  if(string_name_out_prev)  {  delete [] string_name_out_prev,  string_name_out_prev = NULL,  string_name_out_prev_size = 0;}
  if(string_taxon_in_prev)  {  delete [] string_taxon_in_prev,  string_taxon_in_prev = NULL,  string_taxon_in_prev_size = 0;}
  if(string_taxon_out_prev)  { delete [] string_taxon_out_prev, string_taxon_out_prev = NULL, string_taxon_out_prev_size = 0;}
}

//! Deallocates memory of the object given.
void protein_vector::free_list(protein_vector *&obj, uint length) {
  if(obj) {
    for(uint i = 0; i < length; i++)  obj->free_memory();
    free(obj); obj= NULL;
  }
}

protein_vector::protein_vector() :
  string_taxon_in_prev(NULL), string_taxon_out_prev(NULL), 
  string_name_in_prev(NULL),  string_name_out_prev(NULL),  
  string_taxon_in_prev_size(0), string_taxon_out_prev_size(0),
  string_name_in_prev_size(0),  string_name_out_prev_size(0), 
  this_index_in_prev(0),
  this_index_out_prev(0),
  this_taxon_in_prev(0),
  this_taxon_out_prev(0),
  hashProtein(NULL),
  hashTaxa(NULL),
  taxon_length(0),
  is_initiated(false),
  arrOverlap(NULL),
  debug_sum_overlap(0),
  overlaps_inserted_by_curret_call(0)
{
  //  for(int i = 0; i<(int)SIZE_TAXO; i++) string_taxon_in_prev[i] = string_taxon_out_prev[i]= '\0';
  //  for(int i = 0; i<(int)SIZE_PROTEIN; i++) string_name_in_prev[i] = string_name_out_prev[i]= '\0';

}

protein_vector::protein_vector(prot_list_t *_hashProtein, prot_list_t *_hashTaxa, int _taxon_length) :
  string_taxon_in_prev(NULL), string_taxon_out_prev(NULL), 
  string_name_in_prev(NULL),  string_name_out_prev(NULL),  
  string_taxon_in_prev_size(0), string_taxon_out_prev_size(0),
  string_name_in_prev_size(0),  string_name_out_prev_size(0), 
  this_index_in_prev(0),
  this_index_out_prev(0),
  this_taxon_in_prev(0),
  this_taxon_out_prev(0),
  hashProtein(_hashProtein),
  hashTaxa(_hashTaxa),
  taxon_length(_taxon_length),
  is_initiated(false),
  arrOverlap(NULL),
  debug_sum_overlap(0),
  overlaps_inserted_by_curret_call(0),
  has_data(false)
{
  arrOverlap = (overlap_t**)malloc(sizeof(overlap_t*)*taxon_length);
  log_builder::test_memory_condition_and_if_not_abort(arrOverlap!=NULL, __LINE__, __FILE__, __FUNCTION__); 
  for(int i = 0; i< taxon_length;i++) {
    arrOverlap[i] = (overlap_t*)malloc(sizeof(overlap_t)*hashProtein[i].getLength());
    log_builder::test_memory_condition_and_if_not_abort(arrOverlap[i]!=NULL, __LINE__, __FILE__, __FUNCTION__); 
    //      memset(arrOverlap[i], 0, hashProtein[i].getLength());
    for(int k = 0; k< hashProtein[i].getLength(); k++) arrOverlap[i][k] = 0;
  }    
}

#ifdef assert_code

void protein_vector::assert_private_parts() {
  Parse p("in_1", "out_1", "taxon_1", "taxon_2", 0, 0, 0);
  // Verifies emptynes
  assert(!p.is_equal(string_name_in_prev, n_in));
  assert(0 != compare_strings(string_name_in_prev, p.getInnerName()));
  assert(0 != compare_strings(string_name_out_prev, p.getOuterName()));
  assert(0 != compare_strings(string_taxon_in_prev, p.getTaxonInName()));
  assert(0 != compare_strings(string_name_out_prev, p.getTaxonOutName()));

  // Tests setting the protein name:
  bool is_inner = true, is_equal = true;
  bool overlap_is_inserted = false;
  if(false) printf("is_inner(%d), is_equal(%d)\n", is_inner, is_equal);
  assert(!set_protein(p, is_inner, !is_inner, is_equal, overlap_is_inserted));
  assert(0 == compare_strings(string_name_in_prev, p.getInnerName()));
  assert(0 != compare_strings(string_name_out_prev, p.getOuterName()));
  assert(0 != compare_strings(string_taxon_in_prev, p.getTaxonInName()));
  assert(0 != compare_strings(string_name_out_prev, p.getTaxonOutName()));

  is_inner = false, is_equal = false;
  assert(!set_protein(p, is_inner, !is_inner, is_equal, overlap_is_inserted));
  assert(0 == compare_strings(string_name_in_prev, p.getInnerName()));
  assert(0 == compare_strings(string_name_out_prev, p.getOuterName()));
  assert(0 != compare_strings(string_taxon_in_prev, p.getTaxonInName()));
  assert(0 != compare_strings(string_taxon_out_prev, p.getTaxonOutName()));

  // Tests setting the taxon name:
  is_inner = true, is_equal = true;
  assert(!set_taxon(p, is_inner, is_equal));
  assert(0 == compare_strings(string_name_in_prev, p.getInnerName()));
  assert(0 == compare_strings(string_name_out_prev, p.getOuterName()));
  assert(0 == compare_strings(string_taxon_in_prev, p.getTaxonInName()));
  assert(0 != compare_strings(string_name_out_prev, p.getTaxonOutName()));

  is_inner = false, is_equal = false;
  assert(!set_taxon(p, is_inner, is_equal));
  assert(0 == compare_strings(string_name_in_prev, p.getInnerName()));
  assert(0 == compare_strings(string_name_out_prev, p.getOuterName()));
  assert(0 == compare_strings(string_taxon_in_prev, p.getTaxonInName()));
  assert(0 == compare_strings(string_taxon_out_prev, p.getTaxonOutName()));

}

/*
  void protein_vector::assert_class(const bool print_info) {

  //  if(false)   print_constants();
  }
*/
#endif
