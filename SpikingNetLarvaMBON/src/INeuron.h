#pragma once

class synapseEnsemble;
class synapticTransmission;

class INeuron
{
public:
	INeuron(void);
	virtual void RegisterAfferent(synapseEnsemble* s); //adds Afferent to array and registers target neuron with it
	virtual void SpikeArrived(synapticTransmission* s); //Called by SynapseEnsemble
	virtual void ActionPotentialEvent(); //Called when neuron reaches threshold
	virtual unsigned short getFireRate(); //Returns the number of Spikes In the previous Elapsed second 
	void StepSimTime();
	virtual int getID(void);

public:
	~INeuron(void);
};
