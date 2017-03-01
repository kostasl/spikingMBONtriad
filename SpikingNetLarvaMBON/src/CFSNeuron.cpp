///*
/// * \class CFSNeuron.cpp
/// * \brief This is a Izhikevich Fast Spiking  neuron -
/// * \date 19 Feb 2013
/// * \author: kostasl
/// */

///\todo Convert to implement INeuron
#include "CFSNeuron.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

extern gsl_rng* r;

using namespace SpikeN;

CFSNeuron::CFSNeuron(unsigned int ID):CNeuron(ID) {
	// TODO Auto-generated constructor stub


    double rnd = gsl_rng_uniform(r);
	mfa = 0.02+0.08*rnd;

	rnd =  gsl_rng_uniform(r);
	mfb = 0.25-0.05*rnd;

	mfc = -65.0;
	mfd = 2.0f;


	mfv = -65.0;
	mfu = mfb*mfv;


}

CFSNeuron::~CFSNeuron() {
	// TODO Auto-generated destructor stub
}

