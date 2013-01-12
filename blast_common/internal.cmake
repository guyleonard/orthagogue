add_library(
  internal
  file_parse.cxx file_struct.cxx list_file_struct.cxx list_file_parse.cxx norm_t.cxx
  protein_relation.cxx prot_list.cxx
  mcl_bunch.cxx mcl_format.cxx
  taxa.cxx
  # Below files belogning to the parse operation:
#  blast_parsing.cxx  
internal_blast.cxx  
#pipe_parse_file.cxx  pipe_parse_merge.cxx  pipe_parse_parse.cxx 
 protein_vector.cxx  taxon_data.cxx
  taxon_list_t.cxx 
 taxon_pair.cxx pipe_binary.cxx id_simil_list.cxx 
  # Below are files belonging to cluster operation:
  pipe_bucket.cxx  pipe_merge.cxx  pipe_norm.cxx  pipe_struct.cxx  pipe_write.cxx  pipe_write_row.cxx 
  # Below are files beldonging to the cmd_list-folder:
#  cmd_list.cxx cmd_argument.cxx
blast_filtering.cxx
#  read_file.cxx
algo_overlap.cxx

)

