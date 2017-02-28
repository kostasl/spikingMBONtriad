
#pragma once
#include "INeuron.h"

#include "isynapseensemble.h"

//Fwd Declaration
//class synapseEnsemble; //Fwd Declaration
class synapticTransmission;
//class INeuron;

class IFNeuron: public INeuron
{
public:
	IFNeuron(float timestep,short ID=0);
    virtual double StepRK_UpdateVm(void);
    virtual void StepSimTime();
    virtual void RegisterAfferent(ISynapseEnsemble* s); //adds Afferent to array and registers target neuron with it
    virtual void SpikeArrived(synapticTransmission* s); //Called by SynapseEnsemble
    virtual void ActionPotentialEvent(); //Called when neuron reaches threshold
    virtual bool PostSpikeOccured();
    virtual void Reset();
    virtual unsigned short getFireRate(); //Returns the number of Spikes In the previous Elapsed second
    virtual int getID(void);
	~IFNeuron(void);
	
private:
	short ClearSpikeList();//Searches and removes expired spikes
	void PushSpikeStack(); //Copy each Spike 1 place down towards the end of the array.
	double SumExInjections(); //Calc Sum(Sj*I(t)).Contribution of each spike from synapse
	double SumExInjections(float Dt); //Calc at a time+Dt
	double SumInhInjections();
	double SumInhInjections(float Dt); //Calc at a time+Dt
	short mID; //A Number to distinguish this Neuron
	double Eex,Ein;
	double Vrest,Vm,Vreset,Vthres;
	double tafm;
	float tafs; //conductance time constant

	float h; //The Simulation Timestep
	float mfPeriodOfSpikeCount; //Used to increment up to A second to measure post rate
	unsigned short msNumberOfSpikesInPeriod; //Count of Post Spike Occurances
	unsigned short msLastFireRate; //The number of spikes during the last second
    ISynapseEnsemble* mAfferents[MAX_AFFERENTS]; //pointer to array of afferent synapses to this neuron
	unsigned int iLastSynapseIndex;
	synapticTransmission* mSpikes[MAX_SPIKES]; //Pointer to Spike Members array
	int iLastSpikeIndex;
	bool bPrespiketoLog; //True When A pre Spike Occurs, reset when the spike is logged
	bool bPostspiketoLog; //True When POst event occurs, Reset when StepNoSpike Happens
	//SONG Model - 
	double gex_song,gin_song; // Conductance state variables as described by song gex and g_inhibitory
	double gmax;

#ifdef DEBUG_LOG
	
	typedef std::ofstream ofile;
	ofile out;

#endif
};
