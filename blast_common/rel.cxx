#include "rel.h"

  //! Prints the given protein pair
  void rel::print(uint protein_in, char **name_in, char **name_out) {
    if(name_in == NULL)
      printf("\tprotein_in=%d, protein_out=%d, similarity=%f \n", protein_in, ind_out, distance);
    else
      printf("\t%10s\t%s \t%f \n", name_in[protein_in], name_out[ind_out], distance);
  }



  //! The assert method for this class
void rel::assert_class(bool print_info) {
    const static char *class_name = "rel_t";
    if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
    ;
#endif
    if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
  }
