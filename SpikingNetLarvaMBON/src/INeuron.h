///* Interface Class which all neuron types need implement
/// \author Kostas Lagogiannis
#pragma once

//class synapseEnsemble<>;
#include "isynapseensemble.h"
#include "synapseEnsemble.h"

//Forward declarations needed
class synapticTransmission;
class ISynapseEnsemble;

class INeuron
{
public:
	INeuron(void);
    virtual void RegisterAfferent(ISynapseEnsemble* s)=0; //adds Afferent to array and registers target neuron with it
    virtual void SpikeArrived(synapticTransmission* s)=0; //Called by SynapseEnsemble
    virtual void ActionPotentialEvent()=0; //Called when neuron reaches threshold
    virtual unsigned short getFireRate()=0; //Returns the number of Spikes In the previous Elapsed second
    virtual void StepSimTime()=0;
    virtual int getID(void)=0;

	~INeuron(void);
};
