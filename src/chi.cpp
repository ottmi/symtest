////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// NOTE PERTAINING TO THE FOLLOWING FOUR FUNCTIONS                            //
//                                                                            //
// Source: gaussian_distribution_tail.c                                       //
//                                                                            //
// Source: chi-square_distribution_tail.c                                     //
//                                                                            //
// Author: Dick Horn (mathretprogr@gmail.com)                                 //
//                                                                            //
// Note: Used with permission from the author (Wednesday, 9 July 2014 4:30)   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>

// This function returns the probability that a random variable with
// a standard Normal (Gaussian) distribution has a value greater than "x".

long double xGaussian_Distribution_Tail( long double x )
{
   long double sqrt2 = 0.7071067811865475244008443621048490L;
   return  0.5L * erfcl(sqrt2 * x );
}


// The number of degrees of freedom, nu, is an even integer, nu = 2*n.
// The argument x is chi^2 / 2.

static long double Sum_Poisson_Terms(long double x, int n) {
   int k;
   long double term;
   long double sum;

   term = 1.0L;
   sum = 1.0L;
   for (k = 1; k < n; k++) {
      term *= (x / k);
      sum += term;
   }
   return expl(-x)*sum;
}

static long double Sum_Over_Odd_Terms(long double x, int dof)
{
   int k;
   int n;
   long double term;
   long double sum;
   long double sqrtx;
   long double twooverpi;

   twooverpi = 0.6366197723675813430755350534900574L;

   if (dof == 1) return 2.0L * xGaussian_Distribution_Tail( sqrtl(x) );
   n = (dof - 1) / 2;
   sqrtx = sqrtl(x);
   term = sqrtx;
   sum = sqrtx;
   for (k = 2; k <=n; k++) {
      term *= ( x / (k + k - 1) );
      sum += term;
   };

   return 2.0L * xGaussian_Distribution_Tail(sqrtx)
                                      + sqrtl(twooverpi) * expl(-x/2.0L)*sum;
}

// This function returns the probability that a random variable with a standard
// Normal (Gaussian) distribution has a value greater than "x"

long double xChi_Square_Distribution_Tail(long double x, int dof) {

   if (dof <= 0) return 0.0L;

   if (dof % 2 == 0)
		return Sum_Poisson_Terms(x/2.0L, dof/2);
   else
		return Sum_Over_Odd_Terms(x,dof);
}
