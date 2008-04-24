/*
 *  LocalOptimizationND.cpp
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 6/2/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <AnalyticFunction.h>
#include <LocalOptimizationND.h>
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

std::vector<float> dmodule_formula(std::vector<float> const & x)
{
  std::vector<float> result(x.size());

  for (std::size_t i=0; i<x.size(); i++)
  {
	  result[i] = 2 * x[i];
  }
  return result;
}

float f_formula(std::vector<float> const & x)
{
  float result = 0;
  for (std::size_t i=0; i<x.size(); i++)
  {
	  result += cos(x[i]);
  }
  return result;
}

std::vector<float> df_formula(std::vector<float> const & x)
{
  std::vector<float> result(x.size());

  for (std::size_t i=0; i<x.size(); i++)
  {
	  result[i] = - sin(x[i]);
  }
  return result;
}


int main()
{
  LO::AnalyticFunctionFloatND f(f_formula);
  LO::AnalyticGradientFloat df(df_formula);

  LO::LocalOptimizationFloatND opt(f);
  LO::DOptimizationAlgorithmFloatND::pair_t result1 = opt(std::vector<float>(2,0.1));

  std::cout << "Optimization OPT(0.1,0.1) is equal to (x,y) = " << "(" << result1.first[0] << "," << result1.first[1] << ")" << std::endl;
  std::cout << "  number function evaluations : " << f.counter() << std::endl;
  std::cout << "  number derivative evaluations : " << df.counter() << std::endl;
  
  f.resetCounter();
  df.resetCounter();
  
  LO::DLocalOptimizationFloatND dopt(f,df);
  result1 = dopt(std::vector<float>(2,0.1));

  std::cout << "Optimization OPT(0.1,0.1) is equal to (x,y) = " << "(" << result1.first[0] << "," << result1.first[1] << ")" << std::endl;
  std::cout << "  number function evaluations : " << f.counter() << std::endl;
  std::cout << "  number derivative evaluations : " << df.counter() << std::endl;

  return 0;
}
