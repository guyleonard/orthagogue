#ifndef enum_logid_h
#define enum_logid_h
//! @ingroup log_builder
//! \{

/**
   @file
   @brief Identifies the different log types
   @enum logid
   @author Ole Kristian Ekseth (oekseth)
   @date 17.12.2011 by oekseth (initial)
**/
enum logid {
  read_first, collect_read_medio, read_second, 
  collect_read_final, collect_lfp, collect_overlap,
  filter_orthologs_init, filter_orthologs, filter_inparalogs, filter_co_orthologs, write_resultfile,
  complete_runningtime,
  complete_blast_parsing, complete_blast_filtering
};

//! The number of elements in enum 'logid'
static const uint logid_size = 14; // Must be updated iaw above count.

/**
   @brief Identifies the different log types
   @author Ole Kristian Ekseth (oekseth)
   @date 17.12.2011 by oekseth (initial)
**/
typedef enum logid logid_t;

//! \}

#endif
