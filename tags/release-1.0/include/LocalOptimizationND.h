/*
 *  LocalOptimizationND.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 6/2/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  LocalOptimization1D.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/31/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef LOCAL_OPTIMIZATION_ND_H
#define LOCAL_OPTIMIZATION_ND_H

#include <Function.h>
#include <CGOptimization.h>
#include <PowellOptimization.h>

namespace LO {

//! Interface class for ND local optimization.
template <typename D, typename I>
class  LocalOptimizationND : public Function<std::vector<D>,typename OptimizationAlgorithmND<D,I>::pair_t> {

public:

  //! Class identification.
  DefineClassIdentity("LO::LocalOptimizationND");

  //! Constructor by function.
  /*! 
      This constructor creates a local optimization object associated to the given function. 
      
      \param[in] function to be optimized.
  */
  LocalOptimizationND (
    Function<std::vector<D>,I> & f
  ) : Function<std::vector<D>,typename OptimizationAlgorithmND<D,I>::pair_t>(), _function(f)
  {
    _optimization = new PowellOptimization<D,I>;
    _userAlgorithmFlag = false;
  }

  //! Constructor by function and local optimization algorithm.
  /*! 
      This constructor creates a local optimization object associated to the given function. 
      
      \param[in] function to be optimized.
      \param[in] algorithm use in the optimization.
  */
  LocalOptimizationND (
    Function<std::vector<D>,I> & f,
	  OptimizationAlgorithmND<D,I> & o
  ) : Function<std::vector<D>,typename OptimizationAlgorithmND<D,I>::pair_t>(), _function(f) 
  {
    _userAlgorithmFlag = true;
    _optimization = & o;
  }

  //! Destructor.
  ~LocalOptimizationND()
  {
    if (!_userAlgorithmFlag) delete _optimization;
  }

protected:

  virtual typename OptimizationAlgorithmND<D,I>::pair_t evaluate(std::vector<D> const & x)
  {	
    return _optimization->evaluate(_function, x);
  }

private:

  bool _userAlgorithmFlag;
  Function<std::vector<D>,I> & _function;
  OptimizationAlgorithmND<D,I> * _optimization;
  
}; 

typedef LocalOptimizationND<float,float> LocalOptimizationFloatND;
typedef LocalOptimizationND<double,double> LocalOptimizationDoubleND;

//! Interface class for ND local optimization with derivatives.
template <typename D, typename I>
class  DLocalOptimizationND : public Function<std::vector<D>,typename DOptimizationAlgorithmND<D,I>::pair_t> {

public:

  //! Class identification.
  DefineClassIdentity("LO::DLocalOptimizationND");

  //! Constructor by function and derivative.
  /*! 
      This constructor creates a local optimization object associated to the given function. 
      
      \param[in] function to be optimized.
      \param[in] derivative of the function to be optimized.
  */
  DLocalOptimizationND (
    Function<std::vector<D>,I> & f,
    Function<std::vector<D>,std::vector<I> > & df    
  ) : Function<std::vector<D>,typename DOptimizationAlgorithmND<D,I>::pair_t>(), _function(f), _dfunction(df) 
  {
    _optimization = new CGOptimization<D,I>;
    _userAlgorithmFlag = false;
  }

  //! Constructor by function and local optimization algorithm.
  /*! 
      This constructor creates a local optimization object associated to the given function. 
      
      \param[in] function to be optimized.
      \param[in] derivative of the function to be optimized.
      \param[in] algorithm use in the optimization.
  */
  DLocalOptimizationND (
    Function<std::vector<D>,I> & f,
    Function<std::vector<D>,std::vector<I> > & df,    
	  DOptimizationAlgorithmND<D,I> & o
  ) : Function<std::vector<D>,typename DOptimizationAlgorithmND<D,I>::pair_t>(), _function(f), _dfunction(df) 
  {
    _userAlgorithmFlag = true;
    _optimization = & o;
  }

  //! Destructor.
  ~DLocalOptimizationND()
  {
    if (!_userAlgorithmFlag) delete _optimization;
  }

protected:

  virtual typename DOptimizationAlgorithmND<D,I>::pair_t evaluate(std::vector<D> const & x)
  {	
    return _optimization->evaluate(_function, _dfunction, x);
  }

private:

  bool _userAlgorithmFlag;
  Function<std::vector<D>,I> & _function;
  Function<std::vector<D>,std::vector<I> > & _dfunction;
  DOptimizationAlgorithmND<D,I> * _optimization;
  
}; 

typedef DLocalOptimizationND<float,float> DLocalOptimizationFloatND;
typedef DLocalOptimizationND<double,double> DLocalOptimizationDoubleND;

}

#endif
