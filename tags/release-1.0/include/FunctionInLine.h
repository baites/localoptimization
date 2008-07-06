/*
 *  FunctionInLine.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/28/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

 
#ifndef FUNCTION_IN_LINE_H
#define FUNCTION_IN_LINE_H

#include <vector>

#include <ErrorHandler.h>
#include <Function.h>
 
namespace LO {

//! Function values in a given line.
template <typename D, typename I>
class FunctionInLine : public Function<D,I> {
public:

  //! Class identification.
  DefineClassIdentity("LO::FunctionInLine");
  
  //! Constructor by a function.
  /*!
      This constructor only acceps ND function.
      
      \param[in] function to be evaluated in a line.
  */
  FunctionInLine(Function<std::vector<D>,I> & f) : Function<D,I>(), _function(f)
  {
    _lineFlag = false;
  };

  //! Define initial line point and direction.
  /*!
      Initial point and direction must be vectors of the 
      same type and size than the function.
      
      \param[in] direction of the line.
      \param[in] point od the line.
  */  
  void line(
    std::vector<D> const & direction,
    std::vector<D> const & point
  )
  {
    if (direction.size() != point.size())
    {
      ClassError(LOCATION, id(), "evaluate", "Dimesion mismatch between point and direction vectors.");
      return;
    }
    _direction = direction;
    _point = point;
    _lineFlag = true;
  }

protected:

  virtual I evaluate(D const & x);
  
private:

  bool _lineFlag;
  std::vector<D> _point, _direction;
  Function<std::vector<D>,I> & _function;

};

//! Derivative values in a given line.
template <typename D, typename I>
class DFunctionInLine : public Function<D,I> {
public:

  //! Class identification.
  DefineClassIdentity("LO::DFunctionInLine");
  
  //! Constructor by a derivative (gradient).
  /*!
      This constructor only acceps ND function.
      
      \param[in] function to be evaluated in a line.
  */
  DFunctionInLine(Function<std::vector<D>,std::vector<I> > & f) : Function<D,I>(), _dfunction(f)
  {
    _lineFlag = false;
  };

  //! Define initial line point and direction.
  /*!
      Initial point and direction must be vectors of the 
      same type and size than the function.
      
      \param[in] direction of the line.
      \param[in] point od the line.
  */  
  void line(
    std::vector<D> const & direction,
    std::vector<D> const & point
  )
  {
    if (direction.size() != point.size())
    {
      ClassError(LOCATION, id(), "evaluate", "Dimesion mismatch between point and direction vectors.");
      return;
    }
    _direction = direction;
    _point = point;
    _lineFlag = true;
  }

protected:

  virtual I evaluate(D const & x);
  
private:

  bool _lineFlag;
  std::vector<D> _point, _direction;
  Function<std::vector<D>,std::vector<I> > & _dfunction;

};


template<typename D, typename I>
I FunctionInLine<D,I>::evaluate(D const & x)
{
  if (!_lineFlag) ClassError(LOCATION, id(), "evaluate", "No line has been set.");

  std::vector<D> vector(_direction.size());

  for(std::size_t i=0; i<_direction.size(); i++) 
    vector[i] = x * _direction[i] + _point[i];
  
  return _function(vector);
}


template<typename D, typename I>
I DFunctionInLine<D,I>::evaluate(D const & x)
{
  if (!_lineFlag) ClassError(LOCATION, id(), "evaluate", "No line has been set.");

  std::vector<D> vector(_direction.size());

  for(std::size_t i=0; i<_direction.size(); i++) 
    vector[i] = x * _direction[i] + _point[i];
  
  std::vector<I> gradient = _dfunction(vector);
  
  I result = 0;
  for(std::size_t i=0; i<_direction.size(); i++) 
    result += gradient[i] * _direction[i];
  
  return result;
}


}

#endif
