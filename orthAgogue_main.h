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
 * along with orthAgogue. If not, see <http://www.gnu.org/licenses/>.
 */
#include "configure.h" 
#include "blast_parsing.h"  
#include "blast_filtering.h"
#include "taxa.h"
#include "norm_t.h"
#include "list_file_parse.h"
#include "log_builder.h"
#include "tsettings_input.h"
/**
   @file
   @brief The launcher of the orthAgogue
   @ingroup filtering parsing
   @author Ole Kristian Ekseth (oekseth)
   @date 01.01.2012 by oekseth (initial)
**/

// ------------------------Terminal----------------------------

/**
    @defgroup terminal "Terminal Library"
    Handling of terminal inputs.
 */

// ------------------------Log-Builder----------------------------

/**
   @defgroup log_builder "Log Builder"
   Functionality for doing the loging operation with regard to time consumption.
 */


// ------------------------Parsing----------------------------

/**
   @defgroup parsing "Parsing a raw blast"
   Executing code for parsing
 */

/**
   @defgroup pipe_parsing "Thread safe parsing code"
   Parallel blocks used during the blast parsing process.
 */

/**
   @defgroup parsing_container "Parsing Containers"
   Data container used during the blast parsing process.
 */

/**
   @defgroup transfer_container "Transfer Containers"
   Transfer container used in communication between two pipes.
 */

/**
   @defgroup parsing_ops "Parsing Operations"
   Data container used for parsing operations during the blast parsing process.
 */


// ------------------------FILTERING----------------------------

/**
   @defgroup filtering "Filtering a parsed blast"
   Executing code for filtering
 */


/**
   @defgroup pipe_filters "Parallel filtering algorithms"
   Parallel blocks used during the filtering- and output process.
 */


/**
   @defgroup filter_container "Storage of Filtered blast data"
   Data container used during the filter- and output process.
 */

/**
   @defgroup filter_transport_container "Intermediary transport Containers"
   Data container stored intermediary as transported between two pipes.
 */

/**
   @defgroup filtering_ops "Algorithms for serial Filtering"
   Data container used for specific filtering operations.
 */

// ------------------------COMMON----------------------------

/**
   @defgroup blastfile_container "Containers storing input Blast file"
   Data container used storing the blast file.
 */

/**
   @defgroup output_format "Strutured output Formating"
   Data container used buidling the output format required.
 */


// ------------------------MPI---------------------------------

/**
   @defgroup mpi_blast "MPI wrapper code for blast parsing"
   Parallel blocks used during the blast parsing process.
 */


/** 
    @brief The Ececutable of the software
    @ingroup parsing
    @remarks  This file is regarded as the launcher of the software, implying using the two libraries 'blast_parsing' and 'blast_filtering'. In order to understand the aboe libraries, this file might be of use the the developer.
    @author Ole Kristian Ekseth (oekseth)
    @date 27.12.2011 by oekseth (cleanup)
**/
int main_operation_gogue(int array_cnt, char **array);
