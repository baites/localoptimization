/*
 *  AnalyticFunction.cpp
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/27/07.
 *
 */

#include <vector>
#include <iostream>
#include <AnalyticFunction.h>

float module_formula(float const & x)
{
  return x * x;
}

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
  LO::AnalyticFunctionFloat1D pow2(module_formula);

  std::cout << pow2(2.0) << std::endl;

  LO::AnalyticFunctionFloatND module(module_formula);

  std::cout << module(std::vector<float>(2.0,2)) << std::endl;

  return 0;
}

