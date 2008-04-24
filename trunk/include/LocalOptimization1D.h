/*
 *  LocalOptimization1D.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/31/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef LOCAL_OPTIMIZATION_1D_H
#define LOCAL_OPTIMIZATION_1D_H

#include <BrentOptimization.h>
#include <DBrentOptimization.h>
#include <Function.h>
#include <LocalOptimization.h>
#include <SmallestBracket.h>

namespace LO {

//! Interface class for 1D local optimization.
template <typename D, typename I>
class  LocalOptimization1D : public Function<D,typename OptimizationAlgorithm1D<D,I>::pair_t> {

public:

  //! Class identification.  
  DefineClassIdentity("LO::LocalOptimization1D");

  //! Constructor by function.
  /*! 
      This constructor creates a local optimization object associated to the given function. 
      
      \param[in] function to be optimized.
  */
  LocalOptimization1D (
    Function<D,I> & f
  ) : Function<D,typename OptimizationAlgorithm1D<D,I>::pair_t>(), _function(f) 
  {
    _optimization = new BrentOptimization<D,I>;
    _userAlgorithmFlag = false;
  }

  //! Constructor by function and local optimization algorithm.
  /*! 
      This constructor creates a local optimization object associated to the given function. 
      
      \param[in] function to be optimized.
      \param[in] algorithm use in the optimization.
  */
  LocalOptimization1D (
    Function<D,I> & f,
	  OptimizationAlgorithm1D<D,I> & o
  ) : Function<D,typename OptimizationAlgorithm1D<D,I>::pair_t>(), _function(f) 
  {
    _userAlgorithmFlag = true;
    _optimization = & o;
  }

  //! Destructor.
  ~LocalOptimization1D()
  {
    if (!_userAlgorithmFlag) delete _optimization;
  }

protected:

  virtual typename OptimizationAlgorithm1D<D,I>::pair_t evaluate(D const & x)
  {	
    return _optimization->evaluate(_function, x);
  }

private:

  bool _userAlgorithmFlag;
  Function<D,I> & _function;
  OptimizationAlgorithm1D<D,I> * _optimization;
  
}; 


typedef LocalOptimization1D<float,float> LocalOptimizationFloat1D;
typedef LocalOptimization1D<double,double> LocalOptimizationDouble1D;


//! Interface class for 1D local optimization with derivatives.
template <typename D, typename I>
class  DLocalOptimization1D : public Function<D,typename OptimizationAlgorithm1D<D,I>::pair_t> {

public:

  //! Class identification.
  DefineClassIdentity("LO::DLocalOptimization1D");

  //! Constructor by function and derivative.
  /*! 
      This constructor creates a local optimization object associated to the given function. 
      
      \param[in] function to be optimized.
      \param[in] derivative of the function to be optimized.
  */
  DLocalOptimization1D (
    Function<D,I> & f,
    Function<D,I> & df    
  ) : Function<D,typename DOptimizationAlgorithm1D<D,I>::pair_t>(), _function(f), _dfunction(df) 
  {
    _optimization = new DBrentOptimization<D,I>;
    _userAlgorithmFlag = false;
  }

  //! Constructor by function and local optimization algorithm.
  /*! 
      This constructor creates a local optimization object associated to the given function. 
      
      \param[in] function to be optimized.
      \param[in] derivative of the function to be optimized.
      \param[in] algorithm use in the optimization.
  */
  DLocalOptimization1D (
    Function<D,I> & f,
    Function<D,I> & df,        
	  DOptimizationAlgorithm1D<D,I> & o
  ) : Function<D,typename DOptimizationAlgorithm1D<D,I>::pair_t>(), _function(f), _dfunction(df) 
  {
    _userAlgorithmFlag = true;
    _optimization = & o;
  }

  //! Destructor.
  ~DLocalOptimization1D()
  {
    if (!_userAlgorithmFlag) delete _optimization;
  }

protected:

  virtual typename DOptimizationAlgorithm1D<D,I>::pair_t evaluate(D const & x)
  {	
    return _optimization->evaluate(_function, _dfunction, x);
  }

private:

  bool _userAlgorithmFlag;
  Function<D,I> & _function;
  Function<D,I> & _dfunction;
  DOptimizationAlgorithm1D<D,I> * _optimization;
  
}; 

typedef DLocalOptimization1D<float,float> DLocalOptimizationFloat1D;
typedef DLocalOptimization1D<double,double> DLocalOptimizationDouble1D;

}

#endif
