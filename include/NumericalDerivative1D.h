/*
 *  NumericalDerivative.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/27/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef NUMERICAL_DERIVATIVE_1D_H
#define NUMERICAL_DERIVATIVE_1D_H

#include <NumericalDerivative.h>
 
namespace LO {

//! Interface class for numerical derivatives.
template <typename D, typename I>
class  NumericalDerivative1D : public Function<D,I> {

public:

  //! Class identification. 
  DefineClassIdentity("LO::NumericalDerivative1D");

  //! Constructor by function.
  /*!
      This object can only be created by given at least function.
      \param function to calculate the derivative. 
  */
  NumericalDerivative1D(
    Function<D,I> & f
  ) : Function<D,I>(), _function(f) 
  {
    _dfunction = new RiddersDerivative<D,I>;
    _userAlgorithmFlag = false;
    _error = 0;
  }

  //! Constructor by function and derivative algorithm.
  /*!
      This object can only be created by given at least function. 
      In particular this constructor allows to change the algorithm
      use for calculating the derivatives.
      \param function to calculate the derivative.
      \param algorithm use to calculate numerical derivatives. 
  */
  NumericalDerivative1D (
    Function<D,I> & f,
	  DerivativeAlgorithm<D,I> & df
  ) : Function<D,I>(), _function(f) 
  {
    _userAlgorithmFlag = true;
    _dfunction = & df;
    _error = 0;
  }

  //! Destructor.
  ~NumericalDerivative1D()
  {
    if (!_userAlgorithmFlag) delete _dfunction;
  }

  //! Return the error associated to the derivative calculation.
  I error()
  {
    return _error;
  }

protected:

  virtual I evaluate(D const & x)
  {	
    typename DerivativeAlgorithm<D,I>::pair_t result;
    result = _dfunction->evaluate(_function, x);
    _error = result.second;
    return result.first;
  }

private:

  D _error;
  bool _userAlgorithmFlag;
  Function<D,I> & _function;
  DerivativeAlgorithm<D,I> * _dfunction;
  
}; 

typedef NumericalDerivative1D<double,double> NumericalDerivativeDouble1D;
typedef NumericalDerivative1D<float,float> NumericalDerivativeFloat1D;

}

#endif
