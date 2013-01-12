#ifndef pipe_t_h
#define pipe_t_h
/**! Shows the types of operations to use
 */
typedef enum pipe {
  ORTH, // The creation of the possible Ortholog pairs
  INPA, // The creation of the list of inparalogs, based upon the data found in the 'ORTH' step
  INPA_ORTH, // The buling of a set of inpralog-ORTHOLOg_pair set
  INPA_INPA, // The buliding of inpralogs->inparalog set
  DUMP // The dumping of data to a file in ASCII mode (aug. 17: of MCL-format-type)
} pipe_t;


#endif
