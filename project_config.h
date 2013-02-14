//! \{ @ingroup common

/**
   @file
   @brief Configuration settings for OrthaGogue
   @author Ole Kristian Ekseth (oekseth)
**/

//! The name of the directory this projet is build from.
#define PROJECT_BINARY_DIR /home/klatremus/Dokumenter/Work/code/orthAgogue/src
#define FUNC_PROJECT_BINARY_DIR FUNCTION_NAME(PROJECT_BINARY_DIR) // The variable to be called.

//! The major version number of the software:
#define CPACK_PACKAGE_VERSION_MAJOR 1
#ifndef CPACK_PACKAGE_VERSION_MAJOR
#define CPACK_PACKAGE_VERSION_MAJOR 0
#endif

//! The minor version number of the software:
/* #undef CPACK_PACKAGE_VERSION_MINOR */
#ifndef CPACK_PACKAGE_VERSION_MINOR
#define CPACK_PACKAGE_VERSION_MINOR 9
#endif

//! The patch version number of the software:
#define CPACK_PACKAGE_VERSION_PATCH 1
#ifndef CPACK_PACKAGE_VERSION_PATCH
#define CPACK_PACKAGE_VERSION_PATCH 9
#endif

//! Defines if mpi is to be used
/* #undef USE_MPI */

//! The seperator to be used in the output in general
static const char MCL_HEAD_SEPERATOR = ' ';
//! The seperator used after start of the output row.
const static char MCL_SEPERATOR_AFTER_ROW_START = '\t'; 
//! The seperator used after start of rows in the output format
const static char MCL_SEPERATOR_IN_ROWS = '\t'; 
//! The seperator to be used in general.
const static char MCL_SEPERATOR = '\t'; 


// Below macro for handling a konstant as a string:
#define STR_VALUE(arg)      #arg
#define FUNCTION_NAME(name) STR_VALUE(name)
/**
   @brief The level of memory consumption, where '0' means as low as possible
   @remarks If set uses the maximises the speed by using most of the avilable memory on your system. By default set to 1
**/
#define MEMORY_CONSUMPTION_LEVEL 1 // 
#define LOG_FOLDER_NAME report_orthAgogue //report_orthaGogue // Defines the string name
#ifndef LOG_FOLDER_NAME
LOG_FOLDER_NAME report_orthaGogue
#endif
#define FUNC_LOG_FOLDER_NAME FUNCTION_NAME(LOG_FOLDER_NAME) // The variable to be called.

//! The number of Bytes to read from the disk on eeach call to fread(..)
#define BLOCK_FILE_READ_SIZE 5242880
#ifndef BLOCK_FILE_READ_SIZE
#define BLOCK_FILE_READ_SIZE 10485760
//#error A BLOCK_FILE_READ_SIZE value is required!
#endif

/**
   The size in Bytes to send from one pipe to another. Value would affect time- and memory usage.
   @warning If this value is set too low, a non-implemented merging of overlapping blocks will be required!
**/
// #define DISK_BUFFER_SIZE 20971520 
// #ifndef DISK_BUFFER_SIZE
// #error A DISK_BUFFER_SIZE value is required!
// #endif
/**
   @brief Builds the man page for the softwares usage.
   @remarks  Dot not set if buidling this man-page is out of intrest. Else the result file is named turboOrtho.1
**/
/* #undef BUILD_MAN_PAGE */
//! Builds the log files documenting the running time.
#define BUILD_LOG_FILES // If set builds the log files
//! Set if asserts of the functions is to be compiled in.
/* #undef assert_code */
//! If set disables the extra asserts to be run during execution: practical if software behaves in an unexpected manner.
//#define NDEBUG
//#define NDEBUG 
//! Includes extra parameters to the software for debugging purposes (named DEUBG_*).
#define INCLUDE_CMD_DEBUG_PARAMS 

//! Some extra outprints if this is activated.
#define DEBUG 0

/**
	@brief If set write a log file holdgint eh list of memory allocations.
	@remarks: If set generates a file describing the allocation blocks in our parallel building of the blastp data into internal structures used during the filtering process.
**/
#define LOG_WRITE_MEMORY_ALLOCATIONS_DURING_BLASTP_PARSING
/**
	@brief If set write a log file holdgint eh list of memory allocations.
	@remarks: If set generates a file describing the allocation blocks in our parallel building of the blastp data into internal structures used during the filtering process.
**/
#define LOG_WRITE_MEMORY_ALLOCATIONS_DURING_BLASTP_PARSING_FOR_PARSE_BLOCKS

/**
   @brief Produces (writes) a log file holding the list of memory allocations for this object.
   @remarks Useful for verifying the parsing- and uptimizing the internal data structures.
**/
#define LOG_WRITE_MEMORY_ALLOCATIONS_DURING_BLASTP_STRINGS_STORED_IN_MEMORY

/**
   @brief Produces (writes) a log file holding the list of memory allocations for this object.
   @remarks Useful for verifying the parsing- and uptimizing the internal data structures.
**/
#define LOG_WRITE_LIST_FILE_PARSE_CONTENT

/**
   @brief Includes the abc files into the log parsing files
**/
#define LOG_WRITE_LIST_FILE_PARSE_CONTENT_ABC
#ifdef NDEBUG // If in Release mode, discards this option anyhow.
#undef LOG_WRITE_LIST_FILE_PARSE_CONTENT_ABC
#endif


/**
	@brief If set write files describing the allocation table.
	@remarks: If set generates a file describing the allocation blocks in our parallel building of orthology- inparalog- and co-ortholog sequences.
**/
#define LOG_WRITE_SCHEDULING_LIST


/**
	@brief If set write files describing the taxa in the collection.
	@remarks: If set generates a file describing the allocation blocks in our parallel building of orthology- inparalog- and co-ortholog sequences.
**/
#define LOG_WRITE_TAXA_COLLECTION

//! If set prints to stderr those proteins not having an reciprocal pair, and thereby disregarded.
#undef PRINT_PROTEINS_DISREGARDED_DUE_TO_NON_RECIPROCABILITY


/**
   @brief The maximal number of chars a protein label may be identified by.
   @remarks  Estimation of memory consumptiion depends on this variable.
**/
#define SIZE_PROTEIN 23
// Below procedure for future scope.
#ifndef SIZE_PROTEIN 
#define SIZE_PROTEIN 20
#endif

/**
   @brief The maximal number of chars a taxon label may be identified by
   @remarks  Estimation of memory consumptiion depends on this variable.
**/
#define SIZE_TAXO 10 
// Below procedure for future scope.
#ifndef SIZE_TAXO
#define SIZE_TAXO 10
#endif

#ifdef NDEBUG
#undef assert_code
#else
#define assert_code
#endif

// \}
