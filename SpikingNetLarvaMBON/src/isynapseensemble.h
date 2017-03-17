///* \interface ISynapseEnsemble
/// \brief Allows Neurons to be synapse type independent
/// \note Needs empty virtual fn implementation to resolve linker's vtable error
/// \author Kostas Lagogiannis
/// \date 27/2/2007
///
///

#ifndef ISYNAPSEENSEMBLE_H
#define ISYNAPSEENSEMBLE_H

#include "isynapse.h"

//Forward Declaration
class INeuron;

class ISynapseEnsemble
{
public:
    ISynapseEnsemble();
//Pure Virtual Functions
    /// \returns the avg strenght of the all synapses in the ensemble
    virtual float SpikeArrived(ISynapse::SPIKE_SITE type)=0;
    /// \brief Called when no spike event on the simulation step
    virtual void StepNoSpike()=0;
    ///\brief Advances time / used in the RK algorithm
    virtual void StepNoSpike(double t)=0;
    virtual void StepSimTime()=0;
    virtual ISynapse* getSynapseArray()=0;
    /// \returns Number of synapses in Ensemble
    virtual int getSynapsesCount()=0;
    /// \returns the avg strenght of the all synapses in the ensemble
    virtual double getAvgStrength()=0;
    //virtual unsigned short getSourceFireRate()=0;
    virtual void Reset(void)=0;

    virtual ~ISynapseEnsemble();

    ///  \brief When a Neuron registers this SynapseEnsemble it will register also pass a pointer to its self so the SynapseEnsemble can notify the neuron of a spike arrival
    INeuron* getSourceNeuron();
    INeuron* getTargetNeuron();

    void RegisterAfferentNeuron(INeuron* pNeuron); //Target
    void RegisterEfferentNeuron(INeuron* pNeuron); //Source

protected:

    short msourceID; //An Id for the source of this synapse - Target ID is taken from RegisterTarget
    short mtargetID;
    INeuron* mpTargetNeuron;
    INeuron* mpSourceNeuron;


};

#endif // ISYNAPSEENSEMBLE_H
