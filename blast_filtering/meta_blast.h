/**
   @file
   @ingroup filtering
   @brief Holds proeprties about the blast file- and filtering
   @remarks Data only added if explicit set to true.
   @author Ole Kristian Ekseth (oekseth)
**/
#ifndef meta_blast_h
#define meta_blast_h
/*
 * Copyright 2012 Ole Kristian Ekseth (oekseth@gmail.com)
 *
 * This file is part of orthAgogue.
 *
 * orthAgogue is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * orthAgogue is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with orthAgogue
. If not, see <http://www.gnu.org/licenses/>.
 */
/**
   @class meta_blast
   @brief Holds proeprties about the blast file- and filtering
   @ingroup filtering
   @remarks Data only added if explicit set to true.
   @author Ole Kristian Ekseth (oekseth)
**/
class meta_blast {
 private:
  //! @return the percentage difference between two values
  float get_percentage_change(uint cnt_before, uint cnt_after) {
    //! To avoid division by zero.
    if(cnt_before == 0) {
      if(cnt_after == 0) return 0;
      else return FLT_MAX; // It has increased to infinite.
    }
    const float ret_val = (float)(100*(float)(cnt_before - cnt_after)/(float)cnt_before);
    return ret_val;
  }
 public:
  /**
     @brief If not set, the container is not updated
     @remarks For performance measurements, it's preferable having this variable set to false.
  **/
  bool build_meta_blast;
  //! The name of the input file.
  char *FILE_INPUT_NAME;
  //! the number of taxa in the collection.
  uint taxon_length;
  //! The total number of unique proteins, stpread over all the taxa.
  uint cnt_unique_proteins;
  //! The number of putative orthologs before reciprocal matchins is performed.
  uint cnt_relations_orthologs_Before_filtering;
  //! The number of putative orthologs after reciprocal matchins is performed.
  uint cnt_relations_orthologs_After_filtering;
  //! The number of possible inparalogs before threshold filtering is performed, i.e. is the limits generated during the ortholog operation.
  uint cnt_relations_inparalogs_Before_filtering;
  //! The number of possible inparalogs after threshold filtering is performed, i.e. is the limits generated during the ortholog operation.
  uint cnt_relations_inparalogs_After_filtering;
  //! The total number of relations before threshold filtering is applied, i.e. is the limits generated during the ortholog operation.
  uint cnt_relations_Total_Before_filtering;
  //! The total number of relations before co-ortholog filtering is applied.
  uint cnt_relations_Total_After_filtering;
  //! Printes the metadata:
  void print_meta_data(FILE *f) {
    assert(f);
    if(build_meta_blast) {
      uint cnt = 1;
      fprintf(f, "# Processed input file '%s':\n", FILE_INPUT_NAME);
      fprintf(f, "(%u) \t The file had %u taxa- and %u unique proteins in total.\n", cnt++, taxon_length, cnt_unique_proteins);
      float change_factor = get_percentage_change(cnt_relations_Total_Before_filtering, cnt_relations_Total_After_filtering);
      fprintf(f, "(%u) \t The number of possible %s had a %.2f%% change, i.e. from '%u' before filtering, into '%u' after the filtering.\n", cnt++,
	      "relations in total", change_factor, cnt_relations_Total_Before_filtering, cnt_relations_Total_After_filtering);
      //! The inparalogs:
       change_factor = get_percentage_change(cnt_relations_inparalogs_Before_filtering, cnt_relations_inparalogs_After_filtering);
      fprintf(f, "(%u) \t The number of possible %s had a %.2f%% change, i.e. from '%u' before filtering, into '%u' after the filtering.\n", cnt++,
	      "inparalogs", change_factor, cnt_relations_inparalogs_Before_filtering, cnt_relations_inparalogs_After_filtering);
      //! The orthologs:
       change_factor = get_percentage_change(cnt_relations_orthologs_Before_filtering, cnt_relations_orthologs_After_filtering);
      fprintf(f, "(%u) \t The number of possible %s had a %.2f%% change, i.e. from '%u' before filtering, into '%u' after the filtering.\n", cnt++,
	      "orthologs", change_factor, cnt_relations_orthologs_Before_filtering, cnt_relations_orthologs_After_filtering);
    }
  }
  /**
     @brief The constructor.
  **/
  meta_blast(bool _build_meta_blast=false) :
      build_meta_blast(_build_meta_blast), FILE_INPUT_NAME(NULL)
    , taxon_length(0), cnt_unique_proteins(0)
    , cnt_relations_orthologs_Before_filtering(0)
    , cnt_relations_orthologs_After_filtering(0)
    , cnt_relations_inparalogs_Before_filtering(0)
    , cnt_relations_inparalogs_After_filtering(0)
    , cnt_relations_Total_Before_filtering(0)
    , cnt_relations_Total_After_filtering(0)
    {}
};
/**
   @brief Holds proeprties about the blast file- and filtering
   @ingroup filtering
   @remarks Data only added if explicit set to true.
   @author Ole Kristian Ekseth (oekseth)
**/
typedef  meta_blast meta_blast_t;

#endif
