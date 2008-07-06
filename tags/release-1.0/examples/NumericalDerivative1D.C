/*
 *  NumericalDerivative1D.cpp
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/27/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <iomanip>
#include <iostream>

#include <AnalyticFunction.h>
#include <NumericalDerivative1D.h>

float module_formula(float const & x)
{
  return x * x;
}

int main()
{
  LO::AnalyticFunctionFloat1D pow2(module_formula);

  std::cout << "POW2(4.0) is equal to " << pow2(4.0) << std::endl;

  LO::NumericalDerivativeFloat1D dpow2(pow2);
  
  std::cout << "Derivative DPOW2(4.0) calculated by numerical recipe is equal to " << dpow2(4.0) << std::endl;
  
  LO::NPointDerivative<float,float> derivative;
   
  LO::NumericalDerivativeFloat1D dppow2(pow2, derivative);
  
  std::cout << "Derivative DPOW2(4.0) calculated by three point algorithm is equal to " << dppow2(4.0) << std::endl;    

  derivative.npoints = 5;
  
  std::cout << "Derivative DPOW2(4.0) calculated by three point algorithm is equal to " << dppow2(4.0) << std::endl;    
  
  return 0;
}
