
#include "RandomNumbers.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// Code from http://www.firstpr.com.au/dsp/rand31/
//
// adapted by Jean-Francois Cote (jean-francois.cote@fsa.ulaval.ca)
//

/* +++Date last modified: 02-Nov-1995 */

/*
**  longrand() -- generate 2**31-2 random numbers
**
**  public domain by Ray Gardner
**
**  based on "Random Number Generators: Good Ones Are Hard to Find",
**  S.K. Park and K.W. Miller, Communications of the ACM 31:10 (Oct 1988),
**  and "Two Fast Implementations of the 'Minimal Standard' Random
**  Number Generator", David G. Carta, Comm. ACM 33, 1 (Jan 1990), p. 87-88
**
**  linear congruential generator f(z) = 16807 z mod (2 ** 31 - 1)
**
**  uses L. Schrage's method to avoid overflow problems
*/

#define consta 16807
#define PI 3.14159265

int RandomNumbers::randInt(int min, int max)
{
  if(min == max) return min;
	int rnd = next_rand();
	return rnd % (max - min) + min;
}
int RandomNumbers::randIntDifferentThan(int min, int max, int diff)
{
  if(min == max) return min;
  int cpt = 0;
  int num = randInt(min, max);
  while(num == diff && cpt++ < 20)
    num = randInt(min, max);
  return num;
}
double RandomNumbers::rand01()
{
  return next_rand() / 2147483647.0;
}

long unsigned int RandomNumbers::next_rand()
{
  long unsigned int hi, lo;
  lo = consta * (seed31 & 0xFFFF);
  hi = consta * (seed31 >> 16);
  lo += (hi & 0x7FFF) << 16;
  lo += hi >> 15;
  if (lo > 0x7FFFFFFF) lo -= 0x7FFFFFFF;
  return ( seed31 = (long)lo );
}
//===========================================================================
//=  Function to generate normally distributed random variable using the    =
//=  Box-Muller method                                                      =
//=    - Input: mean and standard deviation                                 =
//=    - Output: Returns with normally distributed random variable          =
//===========================================================================
//Taken from: https://cse.usf.edu/~kchriste/tools/gennorm.c
double RandomNumbers::randNormal(double mean, double std_dev)
{
  double   u, sqr, theta;           // Variables for Box-Muller method
  double   x;                     // Normal(0, 1) rv
  double   norm_rv;               // The adjusted normal rv

  // Generate u
  u = 0.0;
  while (u == 0.0)
    u = rand01();

  // Compute r
  sqr = sqrt(-2.0 * log(u));

  // Generate theta
  theta = 0.0;
  while (theta == 0.0)
    theta = 2.0 * PI * rand01();

  // Generate x value
  x = sqr * cos(theta);

  // Adjust x value for specified mean and variance
  norm_rv = (x * std_dev) + mean;

  // Return the normally distributed RV value
  return norm_rv;
}

double RandomNumbers::randBeta(double alpha, double beta)
{
	double x;
	if(alpha <= 1.0 && beta <= 1.0)
		x = randBetaJohnk(alpha,beta);
	else
		x = randBetaMarsagliaTsang(alpha,beta);
	return x;
}	

//Johnk, 1964.
//pp. 416 in: Devroye, L. (1986). Non-Uniform Random Variate Generation. doi:10.1007/978-1-4613-8643-8 
double RandomNumbers::randBetaJohnk(double alpha, double beta) 
{
	double u1, u2, b;
	while(1)
	{
		u1 = rand01();
		u2 = rand01();
		if(std::pow(u1, (1.0 / alpha)) + std::pow(u2, (1.0 / beta)) <= 1.0)
		{
			b = std::pow(u1, (1.0 / alpha)) + std::pow(u2, (1.0 / beta)); break;
		}
	}
    return b;
}

double RandomNumbers::randBetaMarsagliaTsang(double alpha, double beta) 
{

	double gammaAlpha, gammaBeta; //The gamma variables from alpha and beta doubles
	gammaAlpha = randGamma(alpha);
	gammaBeta = randGamma(beta);
	
	//Theorem 4.1 in Devroye, L. (1986) pp. 430: relationship betweer gamma(a), gamma(b) and beta(a,b)
	return gammaAlpha / (gammaAlpha + gammaBeta);

}

// Marsaglia and Tsang's algorithm [2000]
double RandomNumbers::randGamma(double alpha)
{
	double U, x; //The random varibles: uniform and normal.
	double d, c, v; //The values of the formulas
	double gamma; //The gamma variables from alpha and beta doubles
	bool increasedAlpha = false;
	
	//Note in Marsaglia and Tsang's algorithm [2000]: pp.371
	if(alpha < 1.0)
	{
		alpha += 1.0; increasedAlpha = true;
	}
	d = alpha - (1.0/3.0); c=1/(sqrt(9*d));	
	while(1)
	{
		x = randNormal(0,1);  // Generating a standard normal variate (x)
		U = rand01(); // Generating a uniform variate (U)
		v = (1 + c*x) * (1 + c*x) * (1 + c*x);
		if(v <= 0) continue;
		// Check the two conditions for acceptance	
		if( U < 1-0.0331*x*x*x*x )
		{
			gamma = d*v; break;
		}
		if( log(U) < 0.5 * x * x + d * (1 - v + log(v) ))
		{
			gamma = d*v; break;
		}			
	}
	if(increasedAlpha)
	{
		gamma = gamma * ( std::pow(U,1.0/(alpha-1.0)) );
	}
	//printf("G_alpha:%.2lf v:%.2lf d:%.2lf c:%.2lf\n\n",gamma,v,d,c);	
	return gamma;
}
