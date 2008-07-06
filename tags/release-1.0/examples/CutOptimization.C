/*
 *  CutOptimization.cpp
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 2/5/08.
 */

#include <math.h>
#include <iostream>

#include <AnalyticFunction.h>
#include <LocalOptimization1D.h>
#include <LocalOptimizationND.h>
#include <NumericalDerivative1D.h>
#include <NumericalDerivativeND.h>

#define SIG(a) ((a) >= 0.0 ? 1 : -1)

#include "CutOptimization.h"

double R_tip(double const & x)
{
   return V_tip(x)/B_tip(x);
}

double tight = 0.90;

int main()
{
    
  double c, b, v;
  
  std::vector<double> x(3);

  LO::BrentDoubleOptimization brent;

  brent.maxIterations = 3000;
  brent.xTolerance = 1e-8;
  brent.radius = 1e-5;
  
  b = pow(tight,1./4);
  c = inv_B_tip(b);
  v = V_tip(c);
  std::cout << "Cut efficiencies for tip." << std::endl;
  std::cout << "Inversion(Tight,x,B,V,B/V) : " << "(" << tight << "," << c << "," << b << "," << v << "," << (b/v) << ")"  << std::endl;
 
  LO::AnalyticFunctionDouble1D ftip(R_tip);
  LO::NumericalDerivativeDouble1D dftip(ftip);
  LO::LocalOptimizationDouble1D lftip(ftip,brent);
  LO::OptimizationAlgorithmFloat1D::pair_t rtip = lftip(c);
    
  b = B_tip(rtip.first);
  v = V_tip(rtip.first);

  std::cout << "Optimization(x,B,V,B/V) : " << "(" << rtip.first << "," << b << "," << v << "," << (b/v) << ")"  << std::endl;
  return 0;
}
