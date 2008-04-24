/*
 *  BrentOptimization.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 6/1/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DBRENT_OPTIMIZATION_H
#define DBRENT_OPTIMIZATION_H

#include <math.h>

#include <ErrorHandler.h>
#include <Function.h>
#include <LocalOptimization.h>
#include <SmallestBracket.h>

namespace LO {

//! Base class for implemanting 1D local optimization algorithms without derivatives.
template<typename D, typename I>
class DBrentOptimization : public DOptimizationAlgorithm1D<D,I> {

public:

  //! Maximum number of iterations.
  /*!
      By default this value is 300.
  */  
  long	maxIterations;

  //! Tolerance to coordinate changes.
  /*!
      By default this value is 0.001.
  */    
  D xTolerance;

  //! Radius of the interval to sample enviroment of a point.
  /*!
      By default this value is 0.0001.
  */
  D radius;

public:

  //! Class identification.
  DefineClassIdentity("LO::DBrentOptimization");

  //! Void constructor.
  DBrentOptimization() : DOptimizationAlgorithm1D<D,I>() 
  {
    _bracket = new SmallestBracket<D,I>();
    _userBracketFlag = false;    
    maxIterations = 300;
    xTolerance = 0.001;
    radius = 0.0001;
  }

  //! Constructor by SmallestBracket
  /*! 
      This constructor allows to replace bracket algorithm.
        
      \param[in] bracket object that trap the minimum.
  */
  DBrentOptimization(SmallestBracket<D,I> & bracket) : DOptimizationAlgorithm1D<D,I>()
  {
    _userBracketFlag = true;
    _bracket = & bracket;
    maxIterations = 300;
    xTolerance = 0.001;
    radius = 0.0001;
  }

  //! Destructor.  
  ~DBrentOptimization()
  {
    if(!_userBracketFlag) delete _bracket;
  }
  
  //! Evaluate the optimization.
  /*!
      This function implements the optimization.
            
      \param function to optimize.
      \param derivative of the function to optimize.
      \param value of the initial point of the optimization. 
  */
  virtual typename DOptimizationAlgorithm1D<D,I>::pair_t evaluate(Function<D,I> &, Function<D,I> &, D const &);

private:
    
  bool _userBracketFlag;
  SmallestBracket<D,I> * _bracket;
  
  typename DOptimizationAlgorithm1D<D,I>::pair_t brent (Function<D,I> &, Function<D,I> &, D &, D &, D &);
  
};


typedef DBrentOptimization<float,float> DBrentFloatOptimization;
typedef DBrentOptimization<double,double> DBrentDoubleOptimization;


template<typename D, typename I>
typename DOptimizationAlgorithm1D<D,I>::pair_t DBrentOptimization<D,I>::evaluate(Function<D,I> & f, Function<D,I> & df, D const & x)
{  
  D a = x + radius;
  D b = x - radius;  
  D c;

  _bracket->evaluate(f, a, b, c); 
  return brent(f, df, a, b, c);
}


template<typename D, typename I>
typename DOptimizationAlgorithm1D<D,I>::pair_t DBrentOptimization<D,I>::brent(Function<D,I> & f, Function<D,I> & df, D & ax, D & bx, D & cx)
{  
	int ok1, ok2;
  
  D a, b, x, w, v, xm, tol1, tol2, d1, d2, e=0.0;
  D u, u1, u2, olde;
  
	I d,du,dv,dw,dx;
	I fu,fv,fw,fx;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);

	x = w = v = bx;
	fw = fv = fx = f(x);
	dw = dv = dx = df(x);

	for (long iter = 1; iter <= maxIterations; iter++) 
  {
		xm = 0.5 * (a + b);
		tol1 = xTolerance * fabs(x) + ZEPS;
		tol2 = 2.0 * tol1;
    
		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) 
      return typename DOptimizationAlgorithm1D<D,I>::pair_t(x, fx);

		if (fabs(e) > tol1) 
    {
			d1 = 2.0 * (b - a);
			d2 = d1;

			if (dw != dx) d1 = (w - x) * dx / (dx - dw);
			if (dv != dx) d2 = (v - x) * dx / (dx - dv);

			u1 = x + d1;
			u2 = x + d2;
			ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
			ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
			olde = e;
			e = d;

			if (ok1 || ok2)
      {
				if (ok1 && ok2)
					d = (fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d = d1;
				else
					d = d2;
          
				if (fabs(d) <= fabs(0.5 * olde)) 
        {
					u = x + d;
					if (u - a < tol2 || b - u < tol2)
						d = SIGN(tol1, xm - x);
				} 
        else
        {
					d = 0.5 * (e = (dx >= 0.0 ? a-x : b-x));
				}
			} 
      else 
      {
				d = 0.5 * (e = (dx >= 0.0 ? a-x : b-x));
			}
		} 
    else
    {
			d = 0.5 * (e = (dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) 
    {
			u = x + d;
			fu = f(u);
		}
    else
    {
			u = x + SIGN(tol1, d);
			fu = f(u);
			if (fu > fx) 
        return typename DOptimizationAlgorithm1D<D,I>::pair_t(x, fx);
		}

		du = df(u);

		if (fu <= fx) 
    {
			if (u >= x) a=x; else b=x;
			MOV3(v, fv, dv, w, fw, dw)
			MOV3(w, fw, dw, x, fx, dx)
			MOV3(x, fx, dx, u, fu, du)
		} 
    else 
    {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) 
      {
				MOV3(v, fv, dv, w, fw, dw)
				MOV3(w, fw, dw, u, fu, du)
			} 
      else if (fu < fv || v == x || v == w)
 				MOV3(v, fv, dv, u, fu, du)
		}
	}

  ClassError(LOCATION, id(), "evaluate", "Too many iterations.");
  return typename OptimizationAlgorithm1D<D,I>::pair_t(x, fx);
}


}

#endif
