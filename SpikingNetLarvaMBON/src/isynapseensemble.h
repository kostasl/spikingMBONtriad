///* \interface ISynapseEnsemble
/// \brief Allows Neurons to be synapse type independent
/// \note Needs empty virtual fn implementation to resolve linker's vtable error
/// \author Kostas Lagogiannis
/// \date 27/2/2007
///
///

#ifndef ISYNAPSEENSEMBLE_H
#define ISYNAPSEENSEMBLE_H

#include "INeuron.h"

//Forward Declaration
class INeuron;

class ISynapseEnsemble
{
public:
    ISynapseEnsemble();

    virtual float SpikeArrived(ISynapse::SPIKE_SITE type); //Returns the avg strenght of the ensemble
    virtual void StepNoSpike()=0; //Called when no spike event on the simulation step
    virtual void StepNoSpike(double t)=0;
    virtual ISynapse* getSynapseArray()=0;
    virtual int getSynapsesCount()=0;
    virtual double getAvgStrength()=0;
    virtual unsigned short getSourceFireRate()=0;
    virtual void Reset(void)=0;

    ///  \brief When a Neuron registers this SynapseEnsemble it will register also pass a pointer to its self so the SynapseEnsemble can notify the neuron of a spike arrival
    virtual void RegisterNeuron(INeuron* pNeuron)=0;
    virtual short getsourceID()=0;
    virtual short gettargetID()=0;

};

#endif // ISYNAPSEENSEMBLE_H
