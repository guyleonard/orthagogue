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
   Tests the properites by running them independent; handy when testing isolated cases, avoiding the need for remopilation of several layers of makefiles.
*/
#include "blast_common_main.h"
using namespace std;
#include "prot_list.h"
#include "mcl_bunch.h"
#include "mcl_format.h"
#include "taxa.h"
#include <cmph.h>
#include "prot_list.h"
#include "protein_relation.h"
#include "index.h"
#include "norm_t.h"
#include "list_file_parse.h"
#include "file_parse.h"
#include "protein_vector.h"
#include "relations_list.h"
#include "p_rel.h"
#include "rel.h"
int main() {
#ifdef assert_code
  if(true) {
    bool print_info = true;
    mcl_bunch::assert_class(print_info);
    mcl_format::assert_class(print_info);
    index::assert_class(print_info);
    norm::assert_class(print_info);
    file_parse<p_rel>::assert_class(print_info);
    list_file_parse<p_rel>::assert_class(print_info);
    protein_vector::assert_class(print_info);    
    relations_list<rel>::assert_class(print_info);
    relations_list<p_rel>::assert_class(print_info);
    // TOODO: Below code not tested:
    //  protein_relation::assert_class(print_info);
    //    prot_list::assert_class(print_info);
  } 
#endif
}
