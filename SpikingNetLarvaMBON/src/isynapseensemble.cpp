/// \interface ISynapseEnsemble - For collection of some synapse-type bundle
/// \brief Base & Interface class used as a reference point for specific synapse ensemble template instantiations
/// \author kostas lagogiannis
/// \date 28/2/17
///
#include "isynapse.h"
#include "INeuron.h"
#include "isynapseensemble.h"

ISynapseEnsemble::ISynapseEnsemble()
{
 mpTargetNeuron =   0;
 mpSourceNeuron =   0;

 msourceID      =   0;
 mtargetID      =   0;

}

/// \brief Pointer to Source Neuron of this synaptic connection
/// \returns pointer to Afferent Neuron of this synapse
 INeuron* ISynapseEnsemble::getSourceNeuron(){

    return mpSourceNeuron;
}

/// \brief Pointer to Targer Neuron of this synaptic connection
/// \returns pointer to Efferent Neuron of this synapse
 INeuron* ISynapseEnsemble::getTargetNeuron()
{

    return mpTargetNeuron;
}


/// \brief Register a pointer to the Target neuron receiving the impulses - When a Neuron registers this SynapseEnsemble it will register also pass a pointer to its self so the SynapseEnsemble can notify the neuron of a spike arrival
/// \note Called By TargetNeuron
void  ISynapseEnsemble::RegisterAfferentNeuron(INeuron* pTargetN)
{
    mpTargetNeuron = pTargetN; //Synapse now connected to post Neuron
    mtargetID = pTargetN->getID();
}

/// \brief Register Pointer To Source Neuron Providing the impulses - When a Neuron registers this SynapseEnsemble it will register also pass a pointer to its self so the SynapseEnsemble can notify the neuron of a spike arrival
/// \note Called By
void  ISynapseEnsemble::RegisterEfferentNeuron(INeuron* pSourceN)
{
    mpSourceNeuron = pSourceN; //Synapse now connected to post Neuron
    msourceID = pSourceN->getID();
}

ISynapseEnsemble::~ISynapseEnsemble()
{
 mpTargetNeuron =   0;
 mpSourceNeuron =   0;

 msourceID      =   0;
 mtargetID      =   0;

}
