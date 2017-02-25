/*
 * CRSNeuron.cpp
 * This is a Izhikevich Regular Spiking  neuron -
 *  Created on: 19 Feb 2013
 *      Author: kostasl
 */

#include "CRSNeuron.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

extern gsl_rng* r;

using namespace SpikeN;

CRSNeuron::CRSNeuron(unsigned int ID):CNeuron(ID)  {
	// TODO Auto-generated constructor stub
	double rnd =  gsl_rng_uniform(r);
	mfa = 0.02;

	mfb = 0.2;

	mfc = -65.0+15.0*rnd*rnd;

	rnd =  gsl_rng_uniform(r);
	mfd = 8.0 - 6.0*rnd*rnd;

	//initial Values
	mfv = -65.0;
	mfu = mfb*mfv;

}

CRSNeuron::~CRSNeuron() {
	// TODO Auto-generated destructor stub
}

