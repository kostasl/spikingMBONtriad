#pragma once

class PoissonSource
{
public:
	PoissonSource(int lamda,float timeStep,double noiseStdev);
	bool drawSpikeEvent(void);
	double randGauss( double min, double max, double sigma, double centre);
	double gsl_ran_exponential ();
	int getLamda();
	float getRate();
	~PoissonSource(void);

private:
	float mlamda; //Rate of event Arrival
	float h; //Simulation Timestep
	double sigma; //Gaussian Noise StdDev in Sec
	gsl_rng * rng_r; //Used by GSL Rand Num Generator

	time_t t; //Used for random num generation

};
