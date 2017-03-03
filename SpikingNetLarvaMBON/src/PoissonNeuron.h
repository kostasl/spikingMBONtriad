#pragma once
#include "INeuron.h"
#include "isynapseensemble.h"

//Fwd Declaration
//class synapseEnsemble; //Fwd Declaration
class synapticTransmission;
//class INeuron;

class PoissonNeuron:
		public INeuron

{
public:
	PoissonNeuron(float timestep,short ID=0,int StartFireRate=0,bool FixedRate=false);
    virtual void RegisterAfferent(ISynapseEnsemble* s); //adds Afferent to array and registers target neuron with it
    virtual void RegisterEfferent(ISynapseEnsemble* s); //adds Afferent to array and registers target neuron with it
    virtual void SpikeArrived(synapticTransmission* s); //Called by SynapseEnsemble
    virtual void ActionPotentialEvent(); //Called when neuron reaches threshold
    virtual bool ActionPotentialOccured();
    virtual float getFireRate(); //Returns the number of Spikes In the previous Elapsed second
    virtual void StepSimTime();
    virtual int getID(void);

    void setFireRate(float newFireRate);

public:
	~PoissonNeuron(void);

private:
	short mID; //A Number to distinguish this Neuron
	float h; //The Simulation Timestep
    float mlamda; /// \variable mlamda Rate of event Arrival
    double sigma; //Gaussian Noise StdDev in Sec

	float mfPeriodOfSpikeCount; //Used to increment up to A second to measure post rate
	unsigned short msNumberOfSpikesInPeriod; //Count of Post Spike Occurances
    unsigned short msLastFireRate; //The actual since the last call to getFireRate/ number of spikes during the last second
    ISynapseEnsemble* mAfferents[MAX_AFFERENTS]; //pointer to array of afferent synapses to this neuron
	unsigned int iLastSynapseIndex;

    ISynapseEnsemble* mEfferents[MAX_AFFERENTS]; //pointer to array of Efferent (source) synapses to this neuron
    unsigned int iLastSourceSynapseIndex;

	bool mbFixedRate;
    bool bPostspiketoLog; //True When POst event occurs, Reset when StepNoSpike Happens
	time_t t; //Used for random num generation
	gsl_rng * rng_r; //Used by GSL Rand Num Generator

};
