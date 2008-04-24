/*
 *  FunctionInLine.cpp
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/28/07.
 */

#include <vector>
#include <iostream>

#include <AnalyticFunction.h>
#include <FunctionInLine.h>


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


std::vector<float> dmodule_formula(std::vector<float> const & x)
{
  std::vector<float> result(x.size());

  for (std::size_t i=0; i<x.size(); i++)
  {
	  result[i] = 2 * x[i];
  }
  return result;
}


int main()
{
  LO::AnalyticFunctionFloatND module(module_formula);

  std::cout << "module(2.0,2.0) is equal to " << module(std::vector<float>(2,2.0)) << std::endl;

  LO::FunctionInLine<float,float> lmodule(module);
  
  std::vector<float> direction(2,1.0);
  std::vector<float> point(2, 0.0);
  
  lmodule.line(direction,point);

  std::cout << "lmodule(2.0) in (1.0,1.0) direction is equal to " << lmodule(2.0) << std::endl;

  LO::AnalyticGradientFloat dmodule(dmodule_formula);

  LO::DFunctionInLine<float,float> ldmodule(dmodule);

  ldmodule.line(direction,point);

  std::cout << "ldmodule(2.0) in (1.0,1.0) direction is equal to " << ldmodule(2.0) << std::endl;
  
  direction[1] = 0.0;
  point[1] = 2.0;
  
  lmodule.line(direction,point);

  std::cout << "lmodule(2.0) in (1.0,0.0) direction and point (0.0,2.0) is equal to " << lmodule(2.0) << std::endl;

  return 0;
}
