/*
 *  NumericalDerivativeND.cpp
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/28/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>

#include <AnalyticFunction.h>
#include <NumericalDerivativeND.h>

float module_formula(std::vector<float> const & x)
{
  float value, result = 0;
  for (std::size_t i=0; i<x.size(); i++)
  {
	  value = x[i];
	  result += value * value;
  }
  return result;
}

int main()
{
  LO::AnalyticFunctionFloatND module(module_formula);

  std::vector<float> point(2,2.0);

  std::cout << "module(2.0,2.0) is equal to " << module(point) << std::endl;

  LO::NumericalDerivativeFloatND dmodule(module);
  
  std::vector<float> result = dmodule(point);
  
  std::cout << "Gradient dmodule(2.0,2.0) calculated by numerical recipe is equal to ";
  std::cout << "("<< result[0] << "," << result[1] << ")" << std::endl;

  LO::NPointDerivative<float,float> derivative;
   
  LO::NumericalDerivativeFloatND dpmodule(module, derivative);
  
  result = dpmodule(point);
  
  std::cout << "Gradient dmodule(2.0,2.0) calculated by three point differentiation is equal to ";
  std::cout << "("<< result[0] << "," << result[1] << ")" << std::endl;
  
  return 0;
}
