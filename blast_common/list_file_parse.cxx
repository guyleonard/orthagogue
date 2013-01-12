#include "list_file_parse.h"

/**
   Inserts a new relation into the list
   @return true if protein was inserted as a new protein, and not updating an exsisting one.
   @remarks Returnvalue to be used in conjunction with the method named 'getTotalLengthOfData()' found in class list_file_parse.
**/
template <> bool list_file_parse<p_rel_t>::insert_new_rel(uint taxon_in, uint taxon_out, Parse p, protein_relation prot) {
  assert(list);
  assert(list[taxon_in]);
  assert(taxon_in < (uint)taxon_length);
  assert(taxon_out < (uint)taxon_length);
  assert(hashProtein);
  if(list && list[taxon_in]) {
    if(!list[taxon_in][taxon_out]) {
      list[taxon_in][taxon_out] = new file_parse<p_rel>(taxon_in, taxon_out, hashProtein[taxon_in].getLength(), 0,  MAX_BUFFER_SIZE, FILE_BINARY_LOCATION, USE_BEST_BLAST_PAIR_SCORE); // UINT_MAX to avoid writing of data to the file
      if(!list[taxon_in][taxon_out]) { // If memory was not allocated:
	fprintf(stderr, "!!\tMemory was not allocated for taxon-pair[%u, %u]. Aborts. Please refer to the manual or contact the developer at oekset@gmail.com. (This error is found in %s at line %d in method %s\n", taxon_in, taxon_out, __FILE__, __LINE__, __FUNCTION__);
	exit(2);
      }
    }
    return list[taxon_in][taxon_out]->insert_new_rel(prot.protein_in,
						     prot.protein_out,
						     p.get_distance(),
						     p.get_overlap_in(), p.get_overlap_out());
  } else {
    fprintf(stderr, "!!\tlist in class 'list_file_parse' not set at line %d in function %s found in file %s. Please contact the developer, oekseth@gmail.com.\n", __LINE__, __FUNCTION__, __FILE__);
    return false;
  }
}

#ifdef USE_MPI
#include "mpi_transfer_list_file_parse.h"

/**
   @brief Sends only selective data accross nodes.
   @todo Requires that at least the same number of taxa must be found in the file as the number of nodes used for this processing.
**/
template <>   void list_file_parse<p_rel_t>::mpi_send_only_selective_data(MPI_Comm comm, int myrank, int number_of_nodes, int t_taxon_length) {
  assert(taxon_length == t_taxon_length);
  if(!list) return;
  mpi_transfer_list_file_parse<p_rel> mpi_transfer = mpi_transfer_list_file_parse<p_rel>(this);
  mpi_transfer.mpi_send_only_selective_data();
  //  mpi_transfer.mpi_send_only_selective_data(comm, myrank, number_of_nodes, this);
  mpi_transfer.free_memory();
}


/**
   @brief Sends- and receives the list of p_rel_t* objects accross nodes.
   @remarks Each node only sends- and receives those taxon-pairs regarded as interesting in the building of putative orthologs- and inparalogs.
**/
template <>   void list_file_parse<p_rel_t>::mpi_make_data_consistent_accross_nodes(MPI_Comm comm, int myrank, int number_of_nodes, uint index_start, uint index_end_pluss_one) {
  if(!list) return;
  mpi_transfer_list_file_parse<p_rel> mpi_transfer = mpi_transfer_list_file_parse<p_rel>(this);
  mpi_transfer.mpi_make_data_consistent_accross_nodes(comm, myrank, number_of_nodes, index_start, index_end_pluss_one, this);
  mpi_transfer.free_memory();
}

/**
   @brief Sends- and receives the list of rel_t* objects accross nodes.
   @remarks Each node only sends- and receives those taxon-pairs regarded as interesting in the building of co-orthologs.
**/
template <>   void list_file_parse<rel_t>::mpi_make_data_consistent_accross_nodes(int taxon_length) {
  if(!list) return;
  mpi_transfer_list_file_parse<rel> mpi_transfer = mpi_transfer_list_file_parse<rel>(this);
  mpi_transfer.tx_rx_all_data_to_all_nodes();
  //  mpi_transfer.mpi_send_only_selective_data(comm, myrank, number_of_nodes, this);
  //  mpi_transfer.mpi_make_data_consistent_accross_nodes(comm, myrank, number_of_nodes, index_start, index_end_pluss_one, this);
  mpi_transfer.free_memory();
}
#endif
