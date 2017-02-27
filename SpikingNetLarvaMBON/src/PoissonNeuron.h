#pragma once
#include "ineuron.h"
//Fwd Declaration
class synapseEnsemble; //Fwd Declaration
class synapticTransmission;
//class INeuron;

class PoissonNeuron:
		public INeuron

{
public:
	PoissonNeuron(float timestep,short ID=0,int StartFireRate=0,bool FixedRate=false);
	void RegisterAfferent(synapseEnsemble* s); //adds Afferent to array and registers target neuron with it
	void SpikeArrived(synapticTransmission* s); //Called by SynapseEnsemble
	bool StepDrawSpike();
	void ActionPotentialEvent(); //Called when neuron reaches threshold
	void setFireRate(int newFireRate);
	unsigned short getFireRate(); //Returns the number of Spikes In the previous Elapsed second 
	int getID(void);

public:
	~PoissonNeuron(void);

private:
	short mID; //A Number to distinguish this Neuron
	float h; //The Simulation Timestep
	float mfPeriodOfSpikeCount; //Used to increment up to A second to measure post rate
	unsigned short msNumberOfSpikesInPeriod; //Count of Post Spike Occurances
	unsigned short msLastFireRate; //The number of spikes during the last second
	synapseEnsemble* mSynapses[MAX_AFFERENTS]; //pointer to array of afferent synapses to this neuron
	unsigned int iLastSynapseIndex;
	bool mbFixedRate;
	time_t t; //Used for random num generation
	gsl_rng * rng_r; //Used by GSL Rand Num Generator

};
