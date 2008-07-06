/*
 *  AnalyticFunction.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/27/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef ANALYTICFUNCTION_H
#define ANALYTICFUNCTION_H

#include <vector>
#include <Function.h> 
 
namespace LO {

//! Define function with close formula.
template <typename D, typename I>
class AnalyticFunction : public Function<D,I> {

public:

  //! Class identification.   
  DefineClassIdentity("LO::AnalyticFunction");
  
  //! Constructor by pointer to function.
  /*! 
      This cinstructor allows to the define a close function. The formula
      is defined by a common function that is passed as a pointer.
      
      \param[in] pointer to a C function.
  */
  AnalyticFunction(I (*f) (D const &)) : Function<D,I>(), _function(f) {};

protected:

  virtual I evaluate(D const & x)
  {
    return (*_function)(x);
  }

private:

  I (*_function) (D const &);

}; 

typedef AnalyticFunction<double,double> AnalyticFunctionDouble1D;
typedef AnalyticFunction<float,float> AnalyticFunctionFloat1D;
typedef AnalyticFunction<long,float> AnalyticFunctionLong1D;
typedef AnalyticFunction<int,int> AnalyticFunctionInt1D;

typedef AnalyticFunction<std::vector<double>,double> AnalyticFunctionDoubleND;
typedef AnalyticFunction<std::vector<float>,float> AnalyticFunctionFloatND;
typedef AnalyticFunction<std::vector<long>,float> AnalyticFunctionLongND;
typedef AnalyticFunction<std::vector<int>,int> AnalyticFunctionIntND;

typedef AnalyticFunction<std::vector<double>,std::vector<double> > AnalyticGradientDouble;
typedef AnalyticFunction<std::vector<float>,std::vector<float> > AnalyticGradientFloat;
typedef AnalyticFunction<std::vector<long>,std::vector<float> > AnalyticGradientLong;
typedef AnalyticFunction<std::vector<int>,std::vector<int> > AnalyticGradientInt;

}

#endif