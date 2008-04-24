/*
 *  NumericalDerivateND.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/27/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef NUMERICAL_DERIVATIVE_ND_H
#define NUMERICAL_DERIVATIVE_ND_H

#include <FunctionInLine.h>
#include <NumericalDerivative.h>
 
namespace LO {

//! Interface class for numerical gradients.
template <typename D, typename I>
class  NumericalDerivativeND : public Function<std::vector<D>,std::vector<I> > {

public:

  //! Class identification. 
  DefineClassIdentity("LO::NumericalDerivativeND");

  //! Constructor by function.
  /*!
      This object can only be created by given at least a function.
      \param function to calculate the derivative. 
  */
  NumericalDerivativeND(
    Function<std::vector<D>,I> & f
  ) : Function<std::vector<D>,std::vector<I> >(), _function(f) 
  {
    _dfunction = new RiddersDerivative<D,I>;    
    _userAlgorithmFlag = false;
  }

  //! Constructor by function and derivative algorithm.
  /*!
      This object can only be created by given at least function. 
      In particular this constructor allows to change the algorithm
      use for calculating the derivatives.
      \param function to calculate the derivative.
      \param algorithm use to calculate numerical derivatives. 
  */
  NumericalDerivativeND (
    Function<std::vector<D>,I> & f,
    DerivativeAlgorithm<D,I> & df
  ) : Function<std::vector<D>,std::vector<I> >(), _function(f) 
  {
    _userAlgorithmFlag = true;
    _dfunction = & df;
  }

  //! Destructor.
  ~NumericalDerivativeND()
  {
    if (!_userAlgorithmFlag) delete _dfunction;
  }

protected:

  virtual std::vector<I> evaluate(std::vector<D> const & x);

private:

  bool _userAlgorithmFlag;
  DerivativeAlgorithm<D,I> * _dfunction;
  Function<std::vector<D>,I> & _function;
  
}; 

typedef NumericalDerivativeND<float,float> NumericalDerivativeFloatND;
typedef NumericalDerivativeND<double,double> NumericalDerivativeDoubleND;

template<typename D, typename I>
std::vector<I> NumericalDerivativeND<D,I>::evaluate(std::vector<D> const & point)
{
  std::size_t dimension = point.size();
  
  std::vector<D> direction(dimension,0);
  std::vector<I> gradient(dimension);

  typename DerivativeAlgorithm<D,I>::pair_t result;  
  FunctionInLine<D,I> functionInLine(_function);
  
  for (std::size_t i=0; i<dimension; i++)
  {
    if (!i)
      direction[i] = 1;
    else
    {
      direction[i-1] = 0;
      direction[i] = 1;
    }
    
    functionInLine.line(direction, point);
    result = _dfunction->evaluate(functionInLine, 0);
    gradient[i] = result.first;
  }
  return gradient;  
}

}

#endif
