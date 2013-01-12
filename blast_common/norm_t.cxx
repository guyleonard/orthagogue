#include "norm_t.h"

void norm_t::print(uint i, uint out) {
  //printf("norm\t\t[%d][%d]\t%ld\t%f\t%f", i, out, cnt, sum, sum/cnt);
  printf("norm\t\t[%d][%d]\t%d\t%f", i, out, cnt, sum);
  if(cnt_zero > 0) printf("\thas %c zeros", cnt_zero);
  printf("\n");
}

//! Prints the data in form of the labels given in the input file
void norm_t::print_label(struct taxa *listTaxa, uint i, uint out) {
  printf("%s %s has %d elements, with a totalscore of %.3f. ",listTaxa[i].getName(), listTaxa[out].getName(),  cnt, sum);
  if(cnt_zero > 0) printf("\thas %c zeros", cnt_zero);
  printf("\n");
}

//! Prints the data using 'print'
void norm_t::print_basis(struct norm **arrNorm, int taxon_length) { 
  for(uint i=0; i<(uint)taxon_length; i++) {
    for(uint out=0; out<(uint)taxon_length; out++) {
      arrNorm[i][out].print(i, out);
    }
  }
}

/**! Prints the data using the labels found in the input file
   @Changed: 18.03.2011 by oekseth*/
void norm_t::print_labels(uint taxon_length, struct norm **arrNorm, struct taxa *listTaxa) {
  for(uint i=0; i<taxon_length; i++) {
    for(uint out=0; out<(uint)taxon_length; out++) {
      arrNorm[i][out].print_label(listTaxa, i, out);
      arrNorm[i][out].print_label(listTaxa, out, i);
    }
  }
}

/**! Inserts the relations into the structure using local coordinates
   - May print out information if debug information set
   @Changed: 17.11.2011 by oekseth.
*/
void norm_t::insert(uint taxon_index_in, uint taxon_index_out, float arg, uint protein_in, uint protein_out, struct taxa *listTaxa, bool prot_pair_seen_before, const bool DEBUG_NORM) {
  if(arg!= 0) sum+= arg;
  else {
    cnt_zero++;
  }
  if(!prot_pair_seen_before) cnt++; // only updates the count if its a new protein.
  if(DEBUG_NORM) {
    if(listTaxa) printf("%-20s %s(%s)\t%s(%s) \t %.2f.\n","Normbasis:",   listTaxa[taxon_index_in].getArrKey(protein_in), listTaxa[taxon_index_in].getName(), listTaxa[taxon_index_out].getArrKey(protein_out), listTaxa[taxon_index_out].getName(),  arg);
    else printf("%-s %d(%d)\t%d(%d) \t %.2f.\n","Normbasis:",  protein_in, taxon_index_in, protein_out, taxon_index_out,  arg);
  }
}
/**! Inserts the relations into the structure using local coordinates
   - May print out information if debug information set
   @Changed: 18.03.2011 by oekseth.
*/
void norm_t::insert(uint taxon_in, uint taxon_out, float arg, uint protein_in, uint protein_out, struct taxa *listTaxa, const bool DEBUG_NORM) {
  insert(taxon_in, taxon_out, arg, protein_in, protein_out, listTaxa, false, DEBUG_NORM);
}


/*! Takes the max input value found in the blast file pluss one,
  @Changed: 18.03.2011 by oekseth 
*/
void norm_t::changeZeroToValue(float max_input_pluss_one) {
  sum += (float)cnt_zero * max_input_pluss_one; 
  cnt_zero = 0;
}

// Merges two structures given the data.
void norm::merge(uint _cnt, float _sum, uint _zero) {
  cnt += _cnt; sum += sum; cnt_zero += _zero;
}

/*!
  @Name: merge(..) -- Merges two structures
  @Tested: 27.08.2011 by oekseth.
  @Changed: 18.03.2011 by oekseth  
*/
void norm_t::merge(norm n) {
  cnt+= n.cnt, sum+= n.sum, cnt_zero+=n.cnt_zero;
}

void norm_t::merge_basis(uint taxon_length, struct norm **arrNorm, struct norm **&arrNorm_add) {
  if(arrNorm && arrNorm_add) {
    for(uint i=0; i<taxon_length; i++) {
      if(arrNorm[i] && arrNorm_add[i]) {
	for(uint out=0; out<taxon_length; out++) {
	  arrNorm[i][out].merge(arrNorm_add[i][out]);
	}
      }
    }
  }
  //  norm_t::free_arrNorm(taxon_length, arrNorm_add);
}

//! Merges nrom-arrays from several threads into on returning structure.
norm** norm_t::merge_basis_zero(uint n_threads, uint taxon_length, struct norm ***input, float max_input_value, bool PRINT_NORMALIXATION_BASIS, bool DEBUG_NORM) {
  const float zero_add = max_input_value + 1;
  norm_t **norm = norm_t::init_arrNorm(taxon_length, PRINT_NORMALIXATION_BASIS, DEBUG_NORM, max_input_value);
  for(uint in = 0; in<(uint)taxon_length; in++) {
    for(uint out = 0; out<(uint)taxon_length; out++) {
      for(uint i = 0; i< (uint)n_threads; i++) {
	//  const float zero_add = max_input_value + 1;
	input[i][in][out].changeZeroToValue(zero_add), norm[in][out].changeZeroToValue(zero_add); // updates.
	norm[in][out].merge(input[i][in][out]);
      }
    }
  }
  return norm;
}

//! Returns the averaged value
float norm_t::get_basis(const float max_input_value) {
  changeZeroToValue(max_input_value+1); // Change possible zeros to 'normal' values;
  if(cnt > 0 && sum > 0) return (sum/cnt); return 1.0;}

//! Getters sum, transform zeros to 'normal' values.
float norm_t::get_sum(const float max_input_value) {
  changeZeroToValue(max_input_value+1); // Change possible zeros to 'normal' values;
  return sum;
}
uint norm_t::get_cnt()  {return cnt;} // Get the number of proteins occuring in the sample

//! @return true if they are equal:
bool norm::is_equal(norm *obj) {
  if(!obj) assert(obj != NULL);
  if(cnt != obj->get_cnt()) return false;
  if(sum != obj->get_sum()) return false;
  if(cnt_zero != obj->get_cnt_zeros()) return false;
  return true;
}
/**
   @Name: build_basis(..) -- Builds the basis for the normative proecdure.
   -- The code is documenteted more than usual, due to its importancy.
   @TODO: Need some cleanup: probably move most of the outprints to an external function
   @Changed: 18.03.2011 by oekseth (init).
   @Changed: 27.08.2011 by oekseth (asserts).
*/
extern taxa_t *listTaxa;
float** norm_t::build_basis(uint taxon_length, struct norm **arrNorm, bool PRINT_NORMALIXATION_BASIS, float max_input_value) { 
  bool print_debug = false; 
  if(DEBUG || PRINT_NORMALIXATION_BASIS) print_debug = true;
  float **arrAvgNorm = (float**)malloc(sizeof(float*)*taxon_length);
  for(uint i = 0; i< taxon_length; i++) {    arrAvgNorm[i] = (float*)malloc(sizeof(float)*taxon_length);  }
  if(print_debug) {printf("\n# The basis for the normalization procedure (defined in file 'norm_t.cxx' at line %d) is:\n", __LINE__);}
  for(uint i = 0; i< (uint)taxon_length; i++) {
    arrAvgNorm[i][i] = arrNorm[i][i].get_basis(max_input_value); // (sum of sim-scores)/(total number of protein-pairs)
    //    if(print_debug) printf("   arrAvgNorm[%i][%i] = %f\n", i, i, arrNorm[i][i].get_basis());
  }
  for(uint i = 0; i< (uint)taxon_length; i++) { 
    // Iterates over the leftmost (inner) taxon in the pair.
    // (As specified in the documentation, a unique number represents a taxon).
    for(uint out = i+1; out< (uint)taxon_length; out++) {
      // Iterates over the rightmost (outer) taxon in the pair.  (As specified in the documentation, a unique number represents a taxon).
      // The 'arrNorm' contains information about the similarity scores for a given collection of protein pairs.
      if((arrNorm[i][out].get_cnt() + arrNorm[out][i].get_cnt()) > 0) {
	const float zero_add = max_input_value + 1;
	arrNorm[i][out].changeZeroToValue(zero_add); // Changes possible values from zero to real value
	arrNorm[out][i].changeZeroToValue(zero_add); // Changes possible values from zero to real value
	arrAvgNorm[out][i] // Calculates the averaged sum:
	  =  (arrNorm[i][out].get_sum(max_input_value) + arrNorm[out][i].get_sum(max_input_value)) // The 'sum' corresponds to all of the similarity scores summed together
	  /  (arrNorm[i][out].get_cnt() + arrNorm[out][i].get_cnt()) // The 'cnt' corresponds to number of proteins used as a basis for the summing.
	  ;
	arrAvgNorm[i][out] = arrAvgNorm[out][i]; // Makes them reciprocal.
      } else {
	arrAvgNorm[i][out] = 1.0; // No data set
	arrAvgNorm[out][i] = 1.0; // No data set
	// TODO: Document how reciprocal values are set.
      }
      //      if(print_debug) printf("   arrAvgNorm[%i][%i] = %f\n", i, out, arrAvgNorm[i][out]);
    }
  }
  if(print_debug) {
    printf("# The resulting arrAvgNorm (i.e. what the similarity-score(i,j) is to be divided upon) is:\n");
    if(arrAvgNorm) {
      for(uint i = 0; i< (uint)taxon_length; i++) {
	if(arrAvgNorm[i]) {
	  for(uint out = 0; out< (uint)taxon_length; out++) {
	    printf("\tarrAvgNorm[%i][%i] = %f\n", i, out, arrAvgNorm[i][out]);
	  }	  
	} else printf("(not set for [%d])\n", i);
      }
    } else printf("(not set)\n");
  }
  return arrAvgNorm;
}

/**
   @Name: init_arrNorm -- Allocates- and intiates memory for the 'norm_t'structure
   @ProgramPoint: At the start of each thread in the 'void *operator'
*/
struct norm **norm::init_arrNorm(uint taxon_length, bool PRINT_NORMALIXATION_BASIS, bool DEBUG_NORM, float max_input_value) {
  struct norm **arrNorm = (struct norm**)malloc(sizeof(struct norm*)*taxon_length);
  for(uint i = 0; i< taxon_length; i++) {
    arrNorm[i] = (struct norm*)malloc(sizeof(struct norm)*taxon_length);
    for(uint out = 0; out < taxon_length; out++) {
      arrNorm[i][out] = norm();
      //      arrNorm[i][out] = norm(i, out, taxon_length, PRINT_NORMALIXATION_BASIS, DEBUG_NORM, max_input_value);
    }
  }
  return arrNorm;
}
//! Allocates- and intiates memory for the 'norm_t' structure
struct norm ***norm::init(uint taxon_length) {
  assert(taxon_length); // Verifies its set
  norm ***recip_list = new norm**[taxon_length];
  assert(recip_list != NULL);
  if(recip_list) {
    for(uint i = 0; i< (uint)taxon_length; i++) {
      recip_list[i] = new norm*[taxon_length];
      if(recip_list[i]) {
	for(uint out = 0; out < (uint)taxon_length; out++) recip_list[i][out] = NULL;
      } else {
	fprintf(stderr, "!!\tWas not able allocating data due to memory constraints. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", __LINE__, __FILE__, __FUNCTION__);
	exit(2);
      }
      assert(recip_list[i]);
    }
  } else {
    fprintf(stderr, "!!\tWas not able allocating data due to memory constraints. If this error was not as expected, please contact the developer at oekseth@gmail.com adding the following information: Error generated at line %d in file %s for method %s\n", __LINE__, __FILE__, __FUNCTION__);
    exit(2);
  }
  return recip_list;
}

  //! Deallocates the object
void norm::close(norm ***&obj, uint &taxon_length) {
  if(obj) {
    for(uint i = 0; i< (uint)taxon_length; i++) {
      if(obj[i]) {
	for(uint out = 0; out < (uint)taxon_length; out++) {
	  if(obj[i][out]) delete obj[i][out]; obj[i][out] = NULL;
	}
      }
      delete [] obj[i]; obj[i] = NULL;
    }
    delete [] obj; obj = NULL; taxon_length = 0;
  }
}
norm ***norm::init_list(uint list_size, uint element_size, bool PRINT_NORMALIXATION_BASIS, bool DEBUG_NORM, float max_input_value) {
  norm ***local_arrNorm = (norm***)malloc(sizeof(norm**)*list_size);
  for(uint i = 0; i< list_size; i++) local_arrNorm[i] = norm::init_arrNorm(element_size, PRINT_NORMALIXATION_BASIS, DEBUG_NORM, max_input_value);
  return local_arrNorm;
}
void norm::free_arrNorm(uint element_size, norm **&local_arrNorm) {
  if(local_arrNorm != NULL) {
    for(uint taxon_id = 0; taxon_id< element_size; taxon_id++)
      free(local_arrNorm[taxon_id]);
    free(local_arrNorm); local_arrNorm  = NULL;
  }
}
void norm::free_list(uint list_size, uint element_size, norm ***&local_arrNorm) {
  if(local_arrNorm != NULL) {
    for(uint i = 0; i < list_size; i++) {
      norm::free_arrNorm(element_size, local_arrNorm[i]);
    }
    free(local_arrNorm); local_arrNorm = NULL;
  }
}
norm::norm()
  :
  cnt(0), sum(0.0), cnt_zero(0)
{}

//! The constructor.
norm::norm(uint _cnt, float _sum, uint _zero) :
  cnt(_cnt), sum(_sum), cnt_zero(_zero)
{}
void norm::assert_class(const bool print_info) {
  const static char *class_name = "norm";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  //  const float real_max_input_value = max_input_value;
  //  float max_input_value = 900;
  float max_input_value = 0.0;
//   int taxon_length = 1;
//   bool PRINT_NORMALIXATION_BASIS = false;
   bool DEBUG_NORM = false;
   norm temp_1 = norm();
  assert(0 == temp_1.get_sum(max_input_value));  assert(0 == temp_1.get_cnt());
  temp_1.insert(0, 0,40, 0, 0, NULL, DEBUG_NORM);
  assert(40 == temp_1.get_sum(max_input_value));  assert(1 == temp_1.get_cnt());
  temp_1.insert(0, 0,40, 0, 0, NULL, DEBUG_NORM);
  assert(80 == temp_1.get_sum(max_input_value));  assert(2 == temp_1.get_cnt());
  temp_1.insert(0, 0,20, 0, 0, NULL, DEBUG_NORM);
  assert(100 == temp_1.get_sum(max_input_value));  assert(3 == temp_1.get_cnt());
  temp_1.insert(0, 0, 0, 0, 0, NULL, DEBUG_NORM);
  assert(1001 == temp_1.get_sum(max_input_value));  assert(4 == temp_1.get_cnt());

  norm temp_2 = norm();
  temp_2.merge(temp_1);
  assert(temp_2.get_sum(max_input_value) == temp_1.get_sum(max_input_value)); 
  assert(temp_2.get_cnt() == temp_1.get_cnt());
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}


