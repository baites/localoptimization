/*
 *  ConjugateGradient.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 6/2/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CG_OPTIMIZATION_H
#define CG_OPTIMIZATION_H

#include <FunctionInLine.h>
#include <LocalOptimization.h>
#include <DBrentOptimization.h>

namespace LO {

//! Implementation of CG optimization algorithm.
template<typename D, typename I>
class CGOptimization : public DOptimizationAlgorithmND<D,I> {

public:

  //! Maximum number of iterations.
  /*!
      By default this value is 1000.
  */    
  long maxIterations;

  //! Tolerance to function changes.
  /*!
      By default this value is 0.001.
  */    
  D fTolerance;
  
public:

  DefineClassIdentity("LO::CGOptimization");

  //! Void constructor.
  CGOptimization() : DOptimizationAlgorithmND<D,I>() 
  {
    _dLineSearch = new DBrentOptimization<D,I>();
    _lineSearchType = HardCoded;    
    maxIterations = 1000;
    fTolerance = 0.001;
  }

  //! Constructor by OptimizationAlgorithm1D.
  /*!
      Create a CG implementation with a line search algorithm defined by 
      1D optimization without derivative.
      
      \param[in] algorithm for 1D optimization without derivative
  */
  CGOptimization(OptimizationAlgorithm1D<D,I> & lineSearch) : DOptimizationAlgorithmND<D,I>()
  {
    _lineSearchType = LineSearch;
    _lineSearch = & lineSearch;
    maxIterations = 1000;
    fTolerance = 0.001;
  }

  //! Constructor by DOptimizationAlgorithm1D.
  /*!
      Create a CG implementation with a line search algorithm defined by 
      1D optimization with derivative.
      
      \param[in] algorithm for 1D optimization with derivative
  */
  CGOptimization(DOptimizationAlgorithm1D<D,I> & lineSearch) : DOptimizationAlgorithmND<D,I>()
  {
    _lineSearchType = DLineSearch;
    _dLineSearch = & lineSearch;
    maxIterations = 1000;
    fTolerance = 0.001;
  }
  
  //! Destructor.
  ~CGOptimization()
  {
    if (_lineSearchType == HardCoded) 
      delete _dLineSearch;
  }
  
  //! Evaluate the optimization.
  virtual typename DOptimizationAlgorithmND<D,I>::pair_t evaluate(Function<std::vector<D>,I> &, Function<std::vector<D>,std::vector<I> > &, std::vector<D> const &);

private:
 
  enum {
    HardCoded,
    LineSearch,
    DLineSearch
  };
  
  unsigned _lineSearchType;
  OptimizationAlgorithm1D<D,I> * _lineSearch;
  DOptimizationAlgorithm1D<D,I> * _dLineSearch;
    
};


template<typename D, typename I>
typename DOptimizationAlgorithmND<D,I>::pair_t
CGOptimization<D,I>::evaluate(
  Function<std::vector<D>,I> & f, 
  Function<std::vector<D>,std::vector<I> > & df, 
  std::vector<D> const & x
)
{  
  D EPS = 1.0e-10;
	I gg,gam,fp,dgg;

	std::vector<D> g, h, p, xi;
  
  std::size_t dimension = x.size();

  p = x;
	fp = f(p);
  xi = df(p);

	h.resize(dimension);
  g.resize(dimension);
  
	for (std::size_t i = 0; i < dimension; i++)
  {
		g[i] = - xi[i];
		xi[i] = h[i] = g[i];
	}
  
  FunctionInLine<D,I> flin(f);
  DFunctionInLine<D,I> dflin(df);

  typename DOptimizationAlgorithm1D<D,I>::pair_t fo;

	for (int its=1; its <= maxIterations; its++)
  {    
    flin.line(xi, p);
    dflin.line(xi, p);
    
    if(_lineSearchType == LineSearch)
      fo = _lineSearch->evaluate(flin, 0);
    else
      fo = _dLineSearch->evaluate(flin, dflin, 0);

  	for (std::size_t i = 0; i < dimension; i++)
		  p[i] += fo.first * xi[i];
      
		if (2.0 * fabs(fo.second - fp) <= fTolerance * (fabs(fo.second) + fabs(fp) + EPS)) 
			return typename DOptimizationAlgorithmND<D,I>::pair_t(p,fo.second);

		fp = f(p);
    xi = df(p);

		dgg = gg = 0.0;
    for (std::size_t i = 0; i < dimension; i++) 
    {
			gg += g[i] * g[i];
			dgg += (xi[i] + g[i]) * xi[i];
			//dgg += xi(j)*xi(j);
		}

		if (gg == 0.0) return typename DOptimizationAlgorithmND<D,I>::pair_t(p, fp);

		gam = dgg / gg;
    for (std::size_t i=0; i < dimension; i++)
    {
			g[i] = -xi[i];
			xi[i] = h[i] = g[i] + gam * h[i];
		}
	}

  ClassError(LOCATION, id(), "evaluate", "Too many iterations.");
  return typename DOptimizationAlgorithmND<D,I>::pair_t(p,fo.second);
}


}

#endif