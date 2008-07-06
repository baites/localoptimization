/*
 *  NumericalDerivative.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/27/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef NUMERICAL_DERIVATE_H
#define NUMERICAL_DERIVATE_H

#include <math.h>
#include <utility>
#include <vector>

#include <ClassId.H>
#include <Function.h> 
 
namespace LO {

//! Base class for implementing numerical derivatives.
template <typename D, typename I>
class  DerivativeAlgorithm : public ClassId {

public:

  //! Pair type returned by evaluate function.
  typedef std::pair<I,I> pair_t;
  
  //! Class identification.
  DefineClassIdentity("LO::DerivativeAlgorithm");  
  
  //! Void constructor.
  DerivativeAlgorithm(){}
  
  //! Evaluate a derivative.
  /*! 
      This is a pure virtual function for evaluating the derivatives. Therefore
      any numerical implementation of a dervative is forced to define this function.
      
      \param[in] function to calculate the derivative.
      \param[in] point in where the derivative is calculated.
      \return pair_t where first has the actual value of the derivative and second
      has a estimation of the error of the calculation. 
  */
  virtual pair_t evaluate(Function<D,I> &, D const &) = 0;
};


//! Implementation of the three and five point numerical differentation.
template <typename D, typename I>
class NPointDerivative : public DerivativeAlgorithm<D,I> {
public:

  //! Step size.
  /*!
      By default this value is 0.1.
  */
  D step;

  //! Number of point in the calculation.
  /*!
      By default this value is 0.1.
  */
  int npoints;

public:

  //! Class identification. 
  DefineClassIdentity("LO::NPointDerivative");  
  
  //! Void constructor.
  NPointDerivative() : DerivativeAlgorithm<D,I>()
  {
    step = 0.1;
    npoints = 3;
  };

  //! Evaluate a derivative.
  /*! 
      This is a function that implements the derivatives.
        
      \param[in] function to calculate the derivative.
      \param[in] point in where the derivative is calculated.
      \return pair_t where first has the actual value of the derivative and second
      has a estimation of the error of the calculation. 
  */  
  virtual typename DerivativeAlgorithm<D,I>::pair_t 
  evaluate(Function<D,I> &, D const &);

};


//! Implementation of the Ridders differentation.
template <typename D, typename I>
class  RiddersDerivative : public DerivativeAlgorithm<D,I> {
public:

  //! Number of interval "step size" use in the calculation.
  /*!
      By default this value is 10.
  */
  int ntabs;
  
  //! Initial step size.
  /*!
      By default this value is 1.0.
  */
  D step;
  
  //! Convergency safety parameter.
  /*!
      By default this value is 2.0.
  */
  I safety;
  
  //! Scaling factor of the step size.
  /*!
      By default this value is 1.4.
  */
  I factor;
  
public:

  //! Class identification. 
  DefineClassIdentity("LO::RiddersDerivative");  
  
  //! Void constructor.
  RiddersDerivative() : DerivativeAlgorithm<D,I>()
  {
    step = 1.0;
    ntabs = 10;
  	factor = 1.4;
    safety = 2.0;
    _bigNumber = 1e30;
    _factor2 = factor * factor;
  };

  //! Evaluate a derivative.
  /*! 
      This is a function that implements the derivatives.
        
      \param[in] function to calculate the derivative.
      \param[in] point in where the derivative is calculated.
      \return pair_t where first has the actual value of the derivative and second
      has a estimation of the error of the calculation. 
  */  
  virtual typename DerivativeAlgorithm<D,I>::pair_t 
  evaluate(Function<D,I> &, D const &);
  
protected:

  I _factor2;
  I _bigNumber;

  typedef std::vector<D> vector_t;
  typedef typename std::vector<vector_t> matrix_t;

};


template<typename D, typename I>
typename DerivativeAlgorithm<D,I>::pair_t 
NPointDerivative<D,I>::evaluate(
  Function<D,I> & f,
  D const & x
)
{
  D h = step;
  I error = 0;
  I value = 0;
  if (npoints == 3)
    value = ( f(x + h) - f(x - h) ) / (2 * h);
  else if (npoints == 5)
    value = ( f(x - 2 * h) - 8 * f(x - h) + 8 * f(x + h) - f(x + 2 * h) ) / (12 * h);

  return typename DerivativeAlgorithm<D,I>::pair_t(value, error);
}


template<typename D, typename I>
typename DerivativeAlgorithm<D,I>::pair_t 
RiddersDerivative<D,I>::evaluate(
  Function<D,I> & f,
  D const & x
)
{
  D h;
  I value, error, nerror, g;

  vector_t v(ntabs); 
  matrix_t m(ntabs, v);
  
  h = step;
  
  m[0][0] = ( f(x+h) - f(x-h) ) / (2.0*h);

  error = _bigNumber;
  
  for (int i=1; i<ntabs; i++) 
  {
    h /= factor; g = _factor2;
    
    m[0][i] = ( f(x+h) - f(x-h) ) / (2.0*h);
    
    for (int j=1; j<i; j++) 
    {
      m[j][i] = ( m[j-1][i] * g - m[j-1][i-1] ) / (g - 1.0);
      g = _factor2 * g;
      
      nerror = std::max(
        fabs( m[j][i] - m[j-1][i] ),
        fabs( m[j][i] - m[j-1][i-1] )
      );
      
      if (nerror <= error)
      {
        error = nerror;
        value = m[j][i];
      }
	  }
	  if( fabs( m[i][i] - m[i-1][i-1] ) >= safety * error ) break;
  }
  
  return typename DerivativeAlgorithm<D,I>::pair_t(value, error);
}


}

#endif