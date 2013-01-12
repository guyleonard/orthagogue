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
#ifndef mcl_print_settings_h
#define mcl_print_settings_h
/**
   @file
   @struct mcl_print_settings
   @brief Holds the mcl output settings.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
struct mcl_print_settings {
  //! The number of taxa set for the parsed blast file
  int taxon_length;
  //! Information about the taxa collection from the parsed blast file.
  taxa_t *listTaxa;
  //! If unset, uses the raw similarity score for the abc-files.
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_ABC_FORMAT;
  //! If unset, uses the raw similarity score for the mci-files.
  bool DIVIDE_BY_NORMALIZATION_VALUE_FOR_MCL_FORMAT;
  //! for the abc file: if set, the data out pairwise in stead of as in a row
  bool MODE_PAIRWISE_OUTPUT_ABC; 
  //! for the mcl file: if set, the data out pairwise in stead of as in a row
  bool MODE_PAIRWISE_OUTPUT_MCL; 
  //! If unset, the *.abc (i.e. the file with <ype as protein names) discardes
  bool PRINT_IN_ABC_FORMAT;
  //! If unset, the *.mci (i.e. the file with type as protein names) discardes
  bool PRINT_IN_MCL_FORMAT;
  //! If true, sorts teh abc files before outprint
  bool SORT_ABC_DATA; 
};
/**
   @brief Holds the mcl output settings.
   @ingroup output_format
   @author Ole Kristian Ekseth (oekseth)
**/
typedef mcl_print_settings mcl_print_settings_t;
#endif
