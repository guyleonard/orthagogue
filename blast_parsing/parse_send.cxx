#include "parse_send.h"

//! Prints the content of the structs
void parse_send::print_parse_struct() {
  parseData->print_data(false, false);
}

//! Prints the data stored in this object, i.e. the protein relations.
void parse_send::print() {
  if(parseData!= NULL) parseData->print_data(false, false);
  else printf("(list_file_parse_t is NULL)\n");
  prot.print();
  if(arrOverlap != NULL) {
    bool is_not_set = false;
    for(uint i =0; i< (uint)taxon_length; i++) if (arrOverlap[i] == NULL) is_not_set = true;
    if(is_not_set)      printf("arrOverlap is NOT set\n");
    else  printf("arrOverlap is set\n");
  }
  printf("max_sim_value=%f\n", max_sim_value);
}

//! De-allocates the memory reserved for this object.
void parse_send::finalize() {
  if(parseData != NULL) {
    parseData->free_memory(false);
    delete parseData; parseData = NULL;
  }
}

//! The constructor
parse_send::parse_send() :
  taxon_length(0), arrOverlap(NULL), max_sim_value(0),parseData(NULL), found_at_index_pos(0) // moves the data
{}

//! The constructor
parse_send::parse_send(list_file_parse_t *_parseData, protein_relation _prot, overlap_t **arr, float _max_sim_value, int _taxon_length, loint _found_at_index_pos) :
  taxon_length(_taxon_length), arrOverlap(arr), max_sim_value(_max_sim_value),
  parseData(_parseData) // moves the data
  , found_at_index_pos(_found_at_index_pos)
{
  prot.copy(_prot);
}


