set(CXXFILES_BLAST_COMMON
  taxa.cxx parse.cxx 
  list_norm.cxx
  norm_t.cxx
  prot_list.cxx 
  algo_overlap.cxx
  id_simil_list.cxx
  list_file_parse.cxx   file_parse.cxx relations_list.cxx
  mcl_format.cxx mcl_bunch.cxx 
  list_file_chunk.cxx
  protein_relation.cxx # taxon_pair 
  protein_vector.cxx  
  blast_memory.cxx
  index_list.cxx index.cxx rel.cxx p_rel.cxx
  mpi_read.cxx
  )

# Those that not shall be accessable as librries are set here:
set(CXXFILES_BLAST_COMMON_LOCAL )

# set(CXXFILES_BLAST_COMMON
#   taxa.cxx parse.cxx norm_t.cxx
#   prot_list.cxx 
#   algo_overlap.cxx
#   id_simil_list.cxx
#   list_file_struct.cxx  file_struct.cxx list_file_parse.cxx   file_parse.cxx
#   mcl_format.cxx mcl_bunch.cxx protein_relation.cxx # taxon_pair 
#   protein_vector.cxx  
  
#   rel.cxx p_rel.cxx
#   blast_memory.cxx
#   )

# # Those that not shall be accessable as librries are set here:
# set(CXXFILES_BLAST_COMMON_LOCAL 
# #  prot_list.cxx
# #  algo_overlap.cxx  id_simil_list.cxx    list_file_struct.cxx 
# # mcl_format.cxx  protein_relation.cxx  taxa.cxx #taxon_pair.cxx
#  # file_parse.cxx         
# #x norm_t.cxx    
# #  protein_vector.cxx    file_struct.cxx   list_file_parse.cxx  mcl_bunch.cxx   parse.cxx 

# )