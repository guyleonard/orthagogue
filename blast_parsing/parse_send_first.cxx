#include "parse_send_first.h"

//! Deallocates the memory
void parse_send_first::free_taxon_list_mem() {
  if(t_list != NULL) {
    taxon_list::close(t_list);
  //   // t_list->delete_buffer();
  //   // free(t_list); t_list = NULL;
  }
}

//! Deallocates the pointer
void parse_send_first::finalize() {
  if(t_list) {free(t_list); t_list = NULL;}
}

//! The constructor
parse_send_first::parse_send_first(taxon_list_t *list, long long int _last, uint _block_cnt) :
  t_list(list), last_block_pos(_last), block_number(_block_cnt), data_resides_in_mem(false)
{}
