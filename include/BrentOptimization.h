/*
 *  BrentOptimization.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 6/1/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef BRENT_OPTIMIZATION_H
#define BRENT_OPTIMIZATION_H

#include <math.h>

#include <ErrorHandler.h>
#include <Function.h>
#include <LocalOptimization.h>
#include <SmallestBracket.h>

namespace LO {

//! Implementation of the Brent algorithm.
template<typename D, typename I>
class BrentOptimization : public OptimizationAlgorithm1D<D,I> {

public:

  //! Maximum number of iterations.
  /*!
      By default this value is 300.
  */  
  long maxIterations;
  
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
  
  //! Scale factor.
  /*!
      By default this value is 0.1.
  */  
  D cgold;

public:

  //! Class identification.
  DefineClassIdentity("LO::BrentOptimization");

  //! Void constructor.
  BrentOptimization() : OptimizationAlgorithm1D<D,I>() 
  {
    _bracket = new SmallestBracket<D,I>();
    _userBracketFlag = false;    
    maxIterations = 300;
    xTolerance = 0.001;
    radius = 0.0001;
    cgold = 0.1;
  }

  //! Constructor by SmallestBracket
  /*! 
      This constructor allows to replace bracket algorithm.
        
      \param[in] bracket object that trap the minimum.
  */
  BrentOptimization(SmallestBracket<D,I> & bracket) : OptimizationAlgorithm1D<D,I>()
  {
    _userBracketFlag = true;
    _bracket = & bracket;
    maxIterations = 300;
    xTolerance = 0.001;
    radius = 0.0001;
    cgold = 0.1;
  }

  //! Destructor.  
  ~BrentOptimization()
  {
    if(!_userBracketFlag) delete _bracket;
  }
  
  //! Evaluate the optimization.
  /*!
      This function implements the optimization.
            
      \param function to optimize.
      \param value of the initial point of the optimization. 
  */
  virtual typename OptimizationAlgorithm1D<D,I>::pair_t evaluate(Function<D,I> &, D const &);

private:
    
  bool _userBracketFlag;
  SmallestBracket<D,I> * _bracket;
  
  typename OptimizationAlgorithm1D<D,I>::pair_t brent (Function<D,I> &, D &, D &, D &);
  
};

typedef BrentOptimization<float,float> BrentFloatOptimization;
typedef BrentOptimization<double,double> BrentDoubleOptimization;

template<typename D, typename I>
typename OptimizationAlgorithm1D<D,I>::pair_t BrentOptimization<D,I>::evaluate(Function<D,I> & f, D const & x)
{  
  D a = x + radius;
  D b = x - radius;  
  D c;

  _bracket->evaluate(f, a, b, c);
  return brent(f, a, b, c);
}


template<typename D, typename I>
typename OptimizationAlgorithm1D<D,I>::pair_t BrentOptimization<D,I>::brent(Function<D,I> & f, D & ax, D & bx, D & cx)
{  
  D e = 0.0;
  D a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);

  x = w = v = bx;
  fw = fv = fx = f(x);

  for (long iter=1; iter <= maxIterations; iter++) 
  {
    xm = 0.5 * (a + b);
    tol2 = 2.0 * (tol1 = xTolerance * fabs(x) + ZEPS);

    if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) 
      return typename OptimizationAlgorithm1D<D,I>::pair_t(x, fx);

    if (fabs(e) > tol1) 
    {
      r = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v) * q - (x-w) * r;
      q = 2.0 * (q-r);
      
      if (q > 0.0) p = -p;

      q = fabs(q);
      etemp = e;
      e = d;
      
      if (fabs(p) >= fabs(0.5* q * etemp) || p <= q * (a - x) || p >= q * (b - x))
        d = cgold * (e = (x >= xm ? a-x : b-x));
      else 
      {
        d = p / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d=SIGN(tol1,xm-x);
      }
    } 
    else 
    {
      d = cgold * ( e = (x >= xm ? a - x : b - x));
    }
    
    u = (fabs(d) >= tol1 ? x+d : x + SIGN(tol1, d));
    fu=f(u);
    
    if (fu <= fx) 
    {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    } 
    else 
    {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) 
      {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } 
      else if (fu <= fv || v == x || v == w)
      {
        v = u;
        fv = fu;
      }
    }
  }

  ClassError(LOCATION, id(), "evaluate", "Too many iterations.");
  return typename OptimizationAlgorithm1D<D,I>::pair_t(x, fx);
}


}

#endif
