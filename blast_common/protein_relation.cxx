#include "protein_relation.h"
protein_relation::protein_relation(int _taxon_in, mem_loc _protein_in, int _taxon_out, mem_loc _protein_out) : 
  protein_relation_exsists(true), 
  taxon_in(_taxon_in), protein_in(_protein_in),
  taxon_out(_taxon_out), protein_out(_protein_out)
{}
protein_relation::protein_relation(): 
  protein_relation_exsists(true), 
  taxon_in(-1), protein_in(-1), // values (-1) used to check if strucure is set
  taxon_out(0), protein_out(0)
{}

/*
  protein_relation(struct protein_relation prot)  {
  taxon_in = prot.taxon_in, taxon_out = prot.taxon_out;
  protein_in = prot.protein_in, protein_out = prot.protein_out;
  };
*/

//! Prints the fields set
void protein_relation::print() {  printf("(%d)[%d] -> (%d)[%d]\t", taxon_in, (int)protein_in, taxon_out, (int)protein_out);}

bool protein_relation::is_set() {return (taxon_in != -1);}
//! Used in order to knwo what proteins to merge
void protein_relation::copy(protein_relation prot) {
  taxon_in = prot.taxon_in, taxon_out = prot.taxon_out;
  protein_in = prot.protein_in, protein_out = prot.protein_out;
}

void protein_relation::assert_class(const bool print_info) {
  const static char *class_name = "protein_relation";
  if(print_info) printf("--\tStarts testing class %s\n", class_name);
#ifdef assert_code
  protein_relation rel = protein_relation();
  assert(!rel.is_set());
  protein_relation rel_1 = protein_relation(0, 10, 2, 20);
  assert(rel_1.is_set());   
  protein_relation rel_2 = protein_relation(1, 100, 2, 200);
  assert(rel_2.is_set());
  rel_1.copy(rel_2);
  assert(rel_1.get_protein_in() == rel_2.get_protein_in());
  assert(rel_1.get_protein_out() == rel_2.get_protein_out());
  if(print_info) printf("ok\tCompleted testing class %s\n", class_name);
#endif
}




