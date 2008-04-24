/*
 *  PowellOptimization.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 6/6/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef POWELL_OPTIMIZATION_H
#define POWELL_OPTIMIZATION_H

#include <FunctionInLine.h>
#include <LocalOptimization.h>
#include <BrentOptimization.h>

namespace LO {

//! Base class for implemanting 1D local optimization algorithms without derivatives.
template<typename D, typename I>
class PowellOptimization : public OptimizationAlgorithmND<D,I> {

public:
  
  long maxIterations;
  D fTolerance;
  D radius;
  D cgold;

public:

  DefineClassIdentity("LO::PowellOptimization");

  //! Void constructor.
  PowellOptimization() : OptimizationAlgorithmND<D,I>() 
  {
    _lineSearch = new BrentOptimization<D,I>();
    _lineSearchType = HardCoded;    
    maxIterations = 1000;
    fTolerance = 0.01;
  }

  //! Constructor by DOptimizationAlgorithm1D.
  PowellOptimization(OptimizationAlgorithm1D<D,I> & lineSearch) : OptimizationAlgorithmND<D,I>()
  {
    _lineSearchType = LineSearch;
    _lineSearch = & lineSearch;
    maxIterations = 1000;
    fTolerance = 0.01;
  }
  
  ~PowellOptimization()
  {
    if (_lineSearchType == HardCoded) 
      delete _lineSearch;
  }
  
  //! Evaluate the optimization.
  virtual typename OptimizationAlgorithmND<D,I>::pair_t evaluate(Function<std::vector<D>,I> &, std::vector<D> const &);

private:
 
  enum {
    HardCoded,
    LineSearch
  };
  
  unsigned _lineSearchType;
  OptimizationAlgorithm1D<D,I> * _lineSearch;
    
};


template<typename D, typename I>
typename OptimizationAlgorithmND<D,I>::pair_t
PowellOptimization<D,I>::evaluate(
  Function<std::vector<D>,I> & f, 
  std::vector<D> const & x
)
{  
	std::size_t ibig;
  I fp, fptt, fret, del,t;
  
  std::size_t dimension = x.size();
  
  std::vector<D> ptt(dimension);
  std::vector<D> xit(dimension);
  std::vector<D> pt(dimension);
  std::vector<D> p(dimension,0);

  std::vector<std::vector<D> > xi(dimension,p);
  
  for(std::size_t i = 0; i < dimension; i++)
    xi[i][i] = 1;

	fret = f(x);
	pt = x;
  p = x;
 
  FunctionInLine<D,I> flin(f);

  typename OptimizationAlgorithm1D<D,I>::pair_t fo;
    
  for (int iter = 1;;iter++)
  {
		fp = fret;
		ibig = 0;
		del = 0.0;
    
		for (std::size_t i=0; i < dimension; i++) 
    {
      for(std::size_t j = 0; j < dimension; j++)
			  xit[j] = xi[j][i];

			fptt = fret;
      
      flin.line(xit, p);
      fo = _lineSearch->evaluate(flin, 0);
      
      for (std::size_t j = 0; j < dimension; j++)
		    p[j] += fo.first * xit[j];

      fret = fo.second;
	
      if (fabs(fptt - fret) > del)
      {
				del = fabs(fptt - fret);
				ibig = i;
			}
		}
    
		if (2.0 * fabs(fp - fret) <= fTolerance * (fabs(fp) + fabs(fret)))
      return typename OptimizationAlgorithmND<D,I>::pair_t(p, fret);
      
    if (iter == maxIterations)
    {
      ClassError(LOCATION, id(), "evaluate", "Too many iterations.");
      return typename OptimizationAlgorithmND<D,I>::pair_t(p, fret);
    }

    for (std::size_t j = 0; j < dimension; j++)
		{
      ptt[j] = 2.0 * p[j] - pt[j];
		  xit[j] = p[j] - pt[j];
		  pt[j] = p[j];
    }
    
		fptt = f(ptt);
           
		if (fptt < fp) 
    {
			t = 2.0 * (fp - 2.0 * fret + fptt) * sqrt(fp - fret - del) - del * sqrt(fp - fptt);
			if (t < 0.0) 
      {
        flin.line(xit, p);
        fo = _lineSearch->evaluate(flin, 0);
        
        for (std::size_t j = 0; j < dimension; j++)
          p[j] += fo.first * xit[j];

        fret = fo.second;

				for (std::size_t j=0; j < dimension; j++) 
        {
					xi[j][ibig] = xi[j][dimension - 1];
					xi[j][dimension - 1] = xit[j];
				}
			}
		}
	}
}


}

#endif