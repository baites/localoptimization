/*
 *  LocalOptimization.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/31/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef LOCAL_OPTIMIZATION_H
#define LOCAL_OPTIMIZATION_H

#include <utility>
#include <ClassId.h>

namespace LO {


//! Base class for implemanting 1D local optimization algorithms without derivatives.
template<typename D, typename I>
class OptimizationAlgorithm1D : ClassId {

public:

  typedef std::pair<D,I> pair_t;

  DefineClassIdentity("LO::OptimizationAlgorithm1D");

  //! Void constructor.
  OptimizationAlgorithm1D(){};

  //! Evaluate the optimization.
  virtual pair_t evaluate(Function<D,I> &, D const & x) = 0;

};

typedef OptimizationAlgorithm1D<float,float> OptimizationAlgorithmFloat1D;
typedef OptimizationAlgorithm1D<double,double> OptimizationAlgorithmDouble1D;


//! Base class for implemanting 1D local optimization algorithms with derivatives.
template<typename D, typename I>
class DOptimizationAlgorithm1D : ClassId {

public:

  typedef std::pair<D,I> pair_t;

  DefineClassIdentity("LO::DOptimizationAlgorithm1D");

  //! Void constructor.
  DOptimizationAlgorithm1D(){};

  //! Evaluate the optimization.
  virtual pair_t evaluate(Function<D,I> &, Function<D,I> &, D const & x) = 0;
  
};


typedef DOptimizationAlgorithm1D<float,float> DOptimizationAlgorithmFloat1D;
typedef DOptimizationAlgorithm1D<double,double> DOptimizationAlgorithmDouble1D;


//! Base class for implemanting ND local optimization algorithms with derivatives.
template<typename D, typename I>
class OptimizationAlgorithmND : ClassId {

public:

  typedef std::pair<std::vector<D>,I> pair_t;

  DefineClassIdentity("LO::OptimizationAlgorithmND");

  //! Void constructor.
  OptimizationAlgorithmND(){};

  //! Evaluate the optimization.
  virtual pair_t evaluate(Function<std::vector<D>,I> &, std::vector<D> const & x) = 0;
  
};


typedef OptimizationAlgorithmND<float,float> OptimizationAlgorithmFloatND;
typedef OptimizationAlgorithmND<double,double> OptimizationAlgorithmDoubleND;


//! Base class for implemanting ND local optimization algorithms with derivatives.
template<typename D, typename I>
class DOptimizationAlgorithmND : ClassId {

public:

  typedef std::pair<std::vector<D>,I> pair_t;

  DefineClassIdentity("LO::DOptimizationAlgorithmND");

  //! Void constructor.
  DOptimizationAlgorithmND(){};

  //! Evaluate the optimization.
  virtual pair_t evaluate(Function<std::vector<D>,I> &, Function<std::vector<D>,std::vector<I> > &, std::vector<D> const & x) = 0;
  
};

typedef DOptimizationAlgorithmND<float,float> DOptimizationAlgorithmFloatND;
typedef DOptimizationAlgorithmND<double,double> DOptimizationAlgorithmDoubleND;

}

#endif
