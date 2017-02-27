
#include "StdAfx.h"
#include "synapse.h"
#include "synapseEnsemble.h"
#include "synapticTransmission.h"
#include "INeuron.h"
#include "PoissonNeuron.h"

/* An implementation of a Poisson type neuron. The presynaptic rate activity is set and a method
is provided to change it setFireRate during runtime. Averaging would introduce complexities during testing.
The change in rate should be done externally and set to 
// newrate =Pre Rate*Synaptic strength

As every pre neuron is assumed to have the same rate, the update occurs in ActionPotentialEvent
and is mLastFireRate*totalSynapticStrength
Kostantinos Lagogiannis 2007
*/
PoissonNeuron::PoissonNeuron(float timestep,short ID,int StartFireRate,bool FixedRate)
{
	h					= timestep;	
	mID					= ID;
	iLastSynapseIndex	= 0;
	msNumberOfSpikesInPeriod = 0;
	mfPeriodOfSpikeCount	= 0.0f;
	msLastFireRate			= StartFireRate;
	mbFixedRate				= FixedRate;

	//Seed random number generator
	 rng_r = gsl_rng_alloc (gsl_rng_mt19937);
	 if (!rng_r) throw "GSL RNG Init Failed. Out Of Memory?";

//	 unsigned int seed = unsigned(time(&t)) + rand()*100;
	 unsigned int seed = rand()*100;

	 gsl_rng_set(rng_r,seed);

	 //Seed random number generator
	 srand(seed); 

	//srand(unsigned(time(&t)));
}

//adds Afferent to array and registers target neuron with it
//Same As IFNeuron
void PoissonNeuron::RegisterAfferent(synapseEnsemble* pSyn)
{
	//Check that the array is not Maxed.
	if (iLastSynapseIndex == (MAX_AFFERENTS-1))
	{
		std::cout << "Max afferent number reached, Increase MAX_AFFERENTS ";
		abort();
		throw "Max afferent number reached, Increase MAX_AFFERENTS ";
		//Can have code to increase array size with malloc()
	}

	//Register Back to the Synapse Ensemble
	pSyn->RegisterNeuron(this);
	//Add to private array
	mSynapses[iLastSynapseIndex] = pSyn;
	iLastSynapseIndex++;

}

void PoissonNeuron::SpikeArrived(synapticTransmission* s)
{

	//Clear memory from new object
	delete s; 
}
//Called on every Simulation step
bool PoissonNeuron::StepDrawSpike()
{

	//poisson spike draw
	double r = (rand()/(double)RAND_MAX);
	if ((msLastFireRate*h) < r) return false;
	else
	{
		ActionPotentialEvent();
		return true;
	}
}

//Called when neuron Produces output
void PoissonNeuron::ActionPotentialEvent() 
{
	double nRate= 0;
	for(int i=0;i<iLastSynapseIndex;i++)
	{
		if (!mSynapses[i]) break; //Null so Next
		//Notify all synapses (Time steps Fwd?)
		mSynapses[i]->SpikeArrived(synapse::SPIKE_POST);
		
		nRate+= mSynapses[i]->getSourceFireRate()*mSynapses[i]->getAvgStrength();
	}

	//Assume all pre neurons have the same strength
	if (!mbFixedRate) setFireRate(nRate);
	//Simple Linear Relation As in STDP Paper
	//setFireRate(20);
#ifdef VERBOSE
		std::cout << std::endl << "-* Neuron POST Spike new Rate: " << nRate;
#endif

}

//Fire Rate Reset to new Value.
//The change in rate should be done externally and set to 
// newrate =Pre Rate*Synaptic strength
void PoissonNeuron::setFireRate(int newFireRate)
{
	msLastFireRate = newFireRate;
}

//Returns the number of Spikes In the previous Elapsed second 
unsigned short PoissonNeuron::getFireRate()
{
	return msLastFireRate;
}


int PoissonNeuron::getID()
{
	return mID;
}

PoissonNeuron::~PoissonNeuron(void)
{
}
