/*!
  @file
*/
#ifndef o_rel_h
#define o_rel_h
/**
   @brief An 'outer-pair' description of a (pair--{pair}) relation.
   @struct o_rel
   @author Ole Kristian Ekseth (oekseth)
 */
struct o_rel { 
  //! Holds id of outer relation.
  int ind_out; 
  //! The similarity score.
  float distance; 
  //! The constructor.
  o_rel() : ind_out(0), distance(0.0) {};
  //! The constructor.
  o_rel(int _out, float _dist) : ind_out(_out), distance(_dist) {};
};
/**
   @brief An 'outer-pair' description of a (pair--{pair}) relation.
   @author Ole Kristian Ekseth (oekseth)
 */
typedef struct o_rel o_rel_t;
#endif
