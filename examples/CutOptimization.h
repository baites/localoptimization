/*
 *  CutOptimization.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 2/5/08.
 */

double B_sdl(double const & x)
{
  double a = 1.178;
  double b = 1.415;
  double c = 0.4135;
  return pow(x/pow(a+pow(x,b),1.0/b),c);
}

double inv_B_sdl(double const & x)
{
  double a = 1.178;
  double b = 1.415;
  double c = 0.4135;
  return pow(a*pow(x,b/c)/(1-pow(x,b/c)),1/b);
}

double B_dta(double const & x)
{
  double a = 0.003208;
  double b = 1.566;
  double c = 0.6481;
  return pow(x/pow(a+pow(x,b),1.0/b),c);
}

double inv_B_dta(double const & x)
{
  double a = 0.003208;
  double b = 1.566;
  double c = 0.6481;
  return pow(a*pow(x,b/c)/(1-pow(x,b/c)),1/b);
}

double B_lip(double const & x)
{
  double a = 0.03151;
  double b = 1.291;
  double c = 0.426;
  return pow(x/pow(a+pow(x,b),1.0/b),c);
}

double inv_B_lip(double const & x)
{
  double a = 0.03151;
  double b = 1.291;
  double c = 0.426;
  return pow(a*pow(x,b/c)/(1-pow(x,b/c)),1/b);
}

double B_tip(double const & x)
{
  double a = 0.007388;
  double b = 1.546;
  double c = 0.6148;
  return pow(x/pow(a+pow(x,b),1.0/b),c);
}

double inv_B_tip(double const & x)
{
  double a = 0.007388;
  double b = 1.546;
  double c = 0.6148;
  return pow(a*pow(x,b/c)/(1-pow(x,b/c)),1/b);
}

/* double B_pt(double const & x)
{
  if ( x < 0.8 ) return 1;

  double a = 130.8;
  double b = 2.402;
  double c = 0.6867;

  return 1-pow((x-0.8)/pow(a+pow(x-0.8,b),1.0/b),c);
}

double inv_B_pt(double const & x)
{
  if( x == 1 ) return 0.8;

  double a = 130.8;
  double b = 2.402;
  double c = 0.6867;

  return pow(a*pow(1-x,b/c)/(1-pow(1-x,b/c)),1/b)+0.8;
}

double B_chi2(double const & x)
{
  double a = 0.4585;
  double b = 2.393;
  double c = 7.279;
  return pow(x/pow(a+pow(x,b),1.0/b),c);
}

double inv_B_chi2(double const & x)
{
  double a = 0.4585;
  double b = 2.393;
  double c = 7.279;
  return pow(a*pow(x,b/c)/(1-pow(x,b/c)),1/b);
}
*/

double V_sdl(double const & x)
{
  double a = 23.42;
  double b = 1.238;
  double c = 0.18;
  return pow(x/pow(a+pow(x,b),1.0/b),c);
}

double V_dta(double const & x)
{
  double a = 0.1578;
  double b = 1.331;
  double c = 0.3713;
  return pow(x/pow(a+pow(x,b),1.0/b),c);
}

double V_lip(double const & x)
{
  double a = 2.299;
  double b = 1.139;
  double c = 0.157;
  return pow(x/pow(a+pow(x,b),1.0/b),c);
}

double V_tip(double const & x)
{
  double a = 0.05801;
  double b = 1.076;
  double c = 1.231;
  return pow(x/pow(a+pow(x,b),1.0/b),c);
}

/*double V_pt(double const & x)
{
  if ( x < 0.8 ) return 1;

  double a = 3.62;
  double b = 1.849;
  double c = 0.7138;

  return 1-pow((x-0.8)/pow(a+pow(x-0.8,b),1.0/b),c);
}

double V_chi2(double const & x)
{
  double a = 0.09158;
  double b = 1.22;
  double c = 17.46;
  double d = 0.8894;
    
  return d * pow(x/pow(a+pow(x,b),1.0/b),c);
}
*/
