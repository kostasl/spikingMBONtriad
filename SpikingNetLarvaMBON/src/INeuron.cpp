#include "stdafx.h"
#include "isynapse.h"
#include "isynapseensemble.h"
#include "INeuron.h"


////Abstract Class
INeuron::INeuron(void)
{
}

//int INeuron::getID()
//{
//	return 0;
//}

//void INeuron::RegisterAfferent(ISynapseEnsemble* s)
//{
//}
//void INeuron::SpikeArrived(synapticTransmission* s)
//{
//	//Abstract
//}

//void INeuron::ActionPotentialEvent()
//{
//}

////Should be made virtual and implemented by each Neuron Class to step simulation time
////Not implemented
//void INeuron::StepSimTime()
//{
//}

//unsigned short INeuron::getFireRate()
//{
//return 0;
//}



INeuron::~INeuron(void)
{
}
