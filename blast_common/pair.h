#ifndef pair_h
#define pair_h
/**
   @file
   @struct pair_
   @brief Defines a the pair type.
   @ingroup common
   @author Ole Kristian Ekseth (oekseth)
   @date 19.08.2011 by oekseth (initial)
**/
struct pair_ {
  //! The index of the inner protein.
  int ind_in;
  //! The index of the outer protein.
  int ind_out;
  //! Constructor.
  pair_(int _in, int _out) : ind_in(_in), ind_out(_out) {};
  //! Constructor.
  pair_() :  ind_in(0), ind_out(0) {};
};
/**
   @brief Defines a the pair type.
   @ingroup common
   @author Ole Kristian Ekseth (oekseth)
   @date 19.08.2011 by oekseth (initial)
**/
typedef pair_ pair_t;
#endif
