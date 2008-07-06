/*
 *  Function.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/27/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef FUNCTION_H
#define FUNCTION_H
 
#include <ClassId.h> 
  
namespace LO { 

//! Base class of a function. 
template <typename D, typename I>
class Function : public ClassId {

public:
  
  //! Class identification.
  DefineClassIdentity("LO::Function");
  
  //! Domain object type.
  typedef D domain;
  
  //! Image object type.
  typedef I image;
  
  //! Void constructor 
  Function()
  {
    _counter = 0;
  }
  
  //! Evaluate function.
  /*! 
      Operator () overload for function evaluation.
                  
      \param[in] point where to evaluate the function.
      \return value of the function for the given point. 
  */  
  image const & operator()(domain const &x)
  {
    _counter++;
    return _value = evaluate(x);
  }

  //! Reset function evaluation counter.
  void resetCounter()
  {
    _counter = 0;
  }

  //! Number of function evaluations.
  long counter()
  {
    return _counter;
  }
  
protected:  
  
  virtual image evaluate(domain const &) = 0;  
		    
private:

  image _value;
  long _counter;
  
};

}

#endif