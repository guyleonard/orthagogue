#ifndef enum_write_list_h
#define enum_write_list_h
/**
   @file
   @enum write_list;
   @brief In order to identify the files.
   @ingroup filter_container
   @author Ole Kristian Ekseth (oekseth
*/
enum write_list {
  blast_ortho_pairs_list_numbers, // the list of ortho pairs
  blast_inparalog_list_numbers, // the list of the inparalogs
  blast_complete_list_numbers, // the complte list of all the data
  blast_ortho_relations_list_numbers // the list of the ortholog relations
};
/**
   @brief In order to identify the files.
   @ingroup filter_container
   @author Ole Kristian Ekseth (oekseth
*/
typedef enum write_list write_list_t;
/**
   @brief The number of elements in the enum write_list.
   @ingroup filter_container
   @author Ole Kristian Ekseth (oekseth
*/
static const uint size_write_list_t = 4; 
#endif
