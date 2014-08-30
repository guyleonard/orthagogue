#include "algo_overlap.h"
/**
   aboveOverlapLimit -- returns true if data is to be used.
   -- Default input is taken from column 4 in the standard blast output format.
*/
bool aboveOverlapLimit(overlap_t overlap_left_left, overlap_t overlap_left_right,
		       overlap_t overlap_right_right, overlap_t overlap_right_left, short int AMINO_LIMIT, bool PRINT_OVERLAP_VALUES_ABOVE) {
  overlap_t largest_overlap = max(overlap_left_left, overlap_right_right);
  // TODO: Verify the below sentence!
  if(largest_overlap < 1) largest_overlap = 1;
  const overlap_t overlap_proteins = min(overlap_left_right, overlap_right_left);
  const overlap_t final_bastard = (100*overlap_proteins)/largest_overlap;
  // The input gives that the 'overlap_proteins' is greater than the 'largest_overlap' cusing numbers grater than 100:
  if(false) {
    if(!(final_bastard < (100+1))) {
      fprintf(stderr, "final_bastard(%u)\n", final_bastard);
      fprintf(stderr, "overlap_proteins(%u), largest_overlap(%u)\n", overlap_proteins, largest_overlap);
    }
    //  assert((float)(final_bastard < (float)(100+1)));
    // if(!((float)final_bastard < (float)(100+1))) printf("----.sadfsa\fn");
    // printf("__final_bastard(%u)\n", final_bastard);
    //assert(((float)final_bastard < (float)(100+1)));
  }
  //  assert(final_bastard < (100+1));
  if(final_bastard < (overlap_t)AMINO_LIMIT) {

    return false;
  }  else if(final_bastard >= 101) {
    //! Then, if the parssing is correct, two different proteins have a bigger overlap the the proteins overlap with itself, which might be a bit strange.
    return true; // FIXME: validate that this si correct. <-- change made by oekseth at Wed. of August 27. 2014
    // if(PRINT_OVERLAP_VALUES_ABOVE || true) {
    //   printf("overlap=%u due to ", final_bastard);
    //   if(overlap_left_left == largest_overlap) printf("left_protein(%u) > right_protein(%u) ", overlap_left_left, overlap_right_right);
    //   else  printf("right_protein(%u) > left_protein(%u) ", overlap_right_right, overlap_left_left);
    //   printf(" and ");
    //   if(overlap_left_right == overlap_proteins) printf("left_right_protein(%u) < right_left_protein(%u) ", overlap_left_right, overlap_right_left);
    //   else printf("right_left_protein(%u) < left_right_rotein(%u) ", overlap_right_left, overlap_left_right);
    //   printf("\t at %s:%d\t ", __FILE__, __LINE__);
    // }
    // //assert(((float)final_bastard < (float)(100+1)));
    // //    printf("\tabove\t(%u)\t", final_bastard);
    // return false;
  }   else return true; // Conclusion: greater or equal implies holding it';
}

/**
   aboveOverlapLimit_improved -- returns true if data is to be used.
   -- Default input is: diff(column(8)-column_(7)) and diff(column(10)-column_(9))
   Comment: The value could at some point be negative: Do not take the absolute difference.
*/
bool aboveOverlapLimit_improved(overlap_t overlap_left_left, overlap_t overlap_left_right,
				overlap_t overlap_right_right, overlap_t overlap_right_left,   overlap_t recip_left_right_left, overlap_t recip_right_left_right, short int AMINO_LIMIT, bool PRINT_OVERLAP_VALUES_ABOVE) {
  const overlap_t avg_left_right =  (overlap_left_right + recip_right_left_right)/2;
  //  assert(avg_left_right >= 0);
  const overlap_t avg_right_left =  (overlap_right_left + recip_left_right_left)/2;
  //  assert(avg_right_left >= 0);
  return aboveOverlapLimit(overlap_left_left,avg_left_right, overlap_right_right, avg_right_left, AMINO_LIMIT, PRINT_OVERLAP_VALUES_ABOVE);
}
void assert_overlap_module(const bool print_info) {
  const static char *class_name = "algo_overlap";
  if(print_info) printf("--\tStarts testing module %s\n", class_name);
#ifdef assert_code
  // Iinitial tests:
  short int AMINO_LIMIT = 10;
  const short int old_AMINO_LIMIT = AMINO_LIMIT;
  AMINO_LIMIT = 10;
  assert(aboveOverlapLimit(100, 10, 20, 30, AMINO_LIMIT, true));
  //  if(false)   print_constants();
  AMINO_LIMIT = old_AMINO_LIMIT;
#endif
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
}

