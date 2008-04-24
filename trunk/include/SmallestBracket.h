/*
 *  SmallestBracket.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 6/1/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SMALLEST_BRACKET_H
#define SMALLEST_BRACKET_H

#include <iostream>

#include <math.h>
#include <Function.h>

namespace LO {

#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
#define MOV3(a,b,c,d,e,f) (a)=(d);(b)=(e);(c)=(f);

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

//! Class that calculate the smallest bracket with minimum value inside. 
template<typename D, typename I>
struct SmallestBracket
{
  //! Scale factor.
  /*!
      By default this value is 1.618043.
  */  
  D scale;
  
  //! Limit in the number of iterations.
  /*!
      By default this value is 100.
  */    
  D limit;
 
  //! Void constructor
  SmallestBracket()
  {
    scale = 1.618034;
    limit = 100.0;
  }

  //! Find the smallest bracket
  /*
      This function evaluates the smallest bracket with the minimum of the function inside.
      
      \param[in] function in the optimization. 
      \param[in] (a,b,c) initial braket were a < b < c.
  */
  void evaluate(Function<D,I> & f, D & a, D & b, D & c);
};

template<typename D, typename I>
void SmallestBracket<D,I>::evaluate(Function<D,I> & f, D & a, D & b, D & c)
{
  D ulim,u,r,q,fu,dum;
  D TINY = 1.0e-20;

  I fa = f(a);
  I fb = f(b);
  
  if (fb > fa)
  {
    SHFT(dum, a, b, dum);
    SHFT(dum, fb, fa, dum);
  }
  
  c = b + scale * (b - a);
  I fc = f(c);

  while (fb > fc)
  {
    r = (b - a) * (fb - fc);
    q = (b - c) * (fb - fa);
    u = b - ((b - c) * q - (b - a) * r)/(2.0 * SIGN(FMAX(fabs(q-r),TINY), q-r));
    
    ulim = b + limit * (c - b);
    
    if ((b - u) * (u - c) > 0.0) 
    {
      fu = f(u);
      if (fu < fc)
      {
        a = b;
	      b = u;
        fa = fb;
        fb = fu;
        return;
      } 
      else if (fu > fb) 
      {
        c = u;
        fc = fu;
        return;
      }
      u = c + scale * (c - b);
      fu = f(u);
    } 
    else if ((c - u) * (u - ulim) > 0.0) 
    {
      fu = f(u);
      if (fu < fc)
      {
        SHFT(b, c, u, c + scale * (c - b));
        SHFT(fb, fc, fu, f(u));
      }
    } 
    else if ((u - ulim) * (ulim - c) >= 0.0)
    {
      u = ulim;
      fu = f(u);
    } 
    else 
    {
      u = c + scale * (c - b);
      fu = f(u);
    }
    SHFT(a, b, c, u);
    SHFT(fa, fb, fc, fu);
  }
}

}

#endif