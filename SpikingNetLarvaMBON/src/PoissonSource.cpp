/*
############## POISSON SOURCE #########################################
# Implements a spike event source timing trigger, where on each simulation step 
# the source returns if a spike occured stochastically according to a poisson distribution
# parametarised by lamda (ë) the rate of spike occurance in Hz (no of spikes per second). 
# The probability of each event is exponentially distributed in time.
# Kostantinos Lagogiannis 07/2007
########################################################################
*/

#include "StdAfx.h"
#include "PoissonSource.h"

PoissonSource::PoissonSource(int lamda,float timeStep,double noiseStdev)
{
	
	mlamda	= (lamda< 0)?0:(float)lamda;
	h		= timeStep;
	sigma	= noiseStdev;

	//Seed random number generator
	//srand(unsigned(time(&t)));
	 rng_r = gsl_rng_alloc (gsl_rng_mt19937);
	 if (!rng_r) throw "GSL RNG Init Failed. Out Of Memory?";

	 unsigned int seed = unsigned(time(&t)) + rand()*100;
//	 unsigned int seed = rand()*100;


	 gsl_rng_set(rng_r,seed);

	 //Seed random number generator
	 srand(seed); 

}

bool PoissonSource::drawSpikeEvent()
{
	float noise =0; // PoissonSource::randGauss(0,4.0f*sigma,sigma,2.0f*sigma);

	//poisson spike draw
	//double r = (rand()/(double)RAND_MAX);
	double r = gsl_rng_uniform(rng_r);
	if ((mlamda*h + noise)< r) return false;
	else return true;
}

float PoissonSource::getRate()
{
	return mlamda;
}

//gaussian distribution
double PoissonSource::randGauss( double min, double max, double sigma, double centre)
{
/*double random = (min + (max-min) * (double)rand()/RAND_MAX); //create random domain between [min,max]

double tmp = (random-centre)/sigma; 
double gauss = exp(-tmp*tmp/2); //gaussian formula
*/
//Use Chi Square distribution with 2 degrees of freedom
double r1 = (double)rand()/RAND_MAX;
double r2 = (double)rand()/RAND_MAX;

//Note Log =ln
double gauss = sqrt(-2.0f*log(r1))*cos(2*M_PI*r2);

return gauss*sigma;
}

//Returns the time of the next spike under an exponential distribution
double PoissonSource::gsl_ran_exponential ()
{
  double mu = (double)(1.0/(mlamda));
  double u = gsl_rng_uniform_pos (rng_r);

  return -mu * log (u);
}



PoissonSource::~PoissonSource(void)
{
}
