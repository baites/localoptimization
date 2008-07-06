/*
 *  LocalOptimization1D.C
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 6/1/07.
 */

#include <math.h>
#include <iostream>

#include <AnalyticFunction.h>
#include <LocalOptimization1D.h>
#include <NumericalDerivative1D.h>

float pow2_formula(float const & x)
{
  return (x - 2) * (x - 2);
}

float dpow2_formula(float const & x)
{
  return 2 * (x - 2);
}

float cos_formula(float const & x)
{
  return cos(x);
}

float dcos_formula(float const & x)
{
  return - sin(x);
}

int main()
{
  LO::AnalyticFunctionFloat1D pow2(pow2_formula);
  
  LO::LocalOptimizationFloat1D fopt(pow2);
  
  LO::OptimizationAlgorithmFloat1D::pair_t result1 = fopt(1.0);
    
  std::cout << "Optimization FOPT(1.0) is equal to (x,fx) = " << "(" << result1.first << "," << result1.second << ")" << std::endl;
  std::cout << "  number function evaluations : " << pow2.counter() << std::endl;
  
  pow2.resetCounter();
  
  LO::AnalyticFunctionFloat1D dpow2(dpow2_formula);
  
  LO::DLocalOptimizationFloat1D dopt(pow2,dpow2);

  LO::DOptimizationAlgorithmFloat1D::pair_t result2 = dopt(1.0);
    
  std::cout << "Optimization DOPT(1.0) is equal to (x,fx) = " << "(" << result2.first << "," << result2.second << ")" << std::endl;
  std::cout << "  number function evaluations : " << pow2.counter() << std::endl;
  std::cout << "  number derivative evaluations : " << dpow2.counter() << std::endl;

  pow2.resetCounter();

  //LO::NumericalDerivativeFloat1D npow2(pow2);

  LO::NPointDerivative<float,float> derivative;
  LO::NumericalDerivativeFloat1D npow2(pow2,derivative);
 
  LO::DLocalOptimizationFloat1D nopt(pow2,npow2);

  result2 = nopt(1.0);
    
  std::cout << "Optimization DOPT(1.0) is equal to (x,fx) = " << "(" << result2.first << "," << result2.second << ")" << std::endl;
  std::cout << "  number function evaluations : " << pow2.counter() << std::endl;
  std::cout << "  number derivative evaluations : " << npow2.counter() << std::endl;
   
  return 0;
}
