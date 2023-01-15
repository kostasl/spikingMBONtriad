/*
##### SynapticTransmission #####
# Handles the dynamics of synaptic transmission due to a pre synaptic spike event. 
# It represents the current injection at the post-synaptic terminal at time t which is used by the
# I.F neuron to calculate the potential.
# The injection dynamics are governed by alpha(t)=1/taf_s * exp(-t/taf_s)
# gs(t)=epsilon(t)= [tm/(tm-ts)]*[exp(-t/tm)-exp(-t/ts)]
# I(t)=gs(t)(Es-Vm(t)) [Es=Reversal Potential. Excitatory synapses (sodium) Es~0mV]
# In this approach input weights are described as in terms of peak synaptic conductances in pS
# Weights scale the conductance gs(t) which has a peak gmax, 
# sj*gs(t) is returned by getCurrentStep()
############-----Note---------########
# Originaly implemented to represent the lifetime of the conductance channel and added to a stack in the IFNeuron
if the USE_SONG_CONDUCTANCE is defined in the IFNEURON.h then the stack is not used and simpler conductance channel 
update and decay calculation is used as found on the Song et al. paper of 2000. In that case the synaptic transmission 
is object is only used to obtain the synaptic strength at the time of transmission.
*/

#include "stdafx.h"
#include "math.h"
#include "synapticTransmission.h"
#include "isynapseensemble.h"


synapticTransmission::synapticTransmission(double SynapseStrength)
{
	mbTransmitted = false;
    mfSimTimestep = TIMESTEP;
	mdtime= 0; //Start at time 0 and increment by timestep
	mbTransmitted = false;
	tafm = 0.020f;
	tafs = 0.005f;
	//tau = tafm/(tafm-tafs);
	mfSynapticStrength = SynapseStrength;
}

/// \brief Use this constructor passing the information on which Synapse is propagating the spike
synapticTransmission::synapticTransmission(ISynapseEnsemble* pSourceSynapse)
{
    mbTransmitted = false;
    mfSimTimestep = TIMESTEP;
    mdtime= 0; //Start at time 0 and increment by timestep
    mbTransmitted = false;
    tafm = 0.020f;
    tafs = 0.005f;
    //tau = tafm/(tafm-tafs);
    mfSynapticStrength = pSourceSynapse->getAvgStrength();
    mpSourceSynapse     = pSourceSynapse;
}

//Returns the value of the current injection*Synaptic Strength
//Note: Synaptic strength is assumed to be constant during transmission until the injection is finished.
double synapticTransmission::getCurrentStep()
{
	//Step the time fwd
	mdtime+=mfSimTimestep;

	///If 3 time const have past then the injection has decayed by 90%, So the spike has expired
	mbTransmitted = (mdtime/tafm) > 3;

	return getValueAt(0);
}
//Returns the Value of gex evaluated at a future time Dt
double synapticTransmission::getValueAt(float Dt)
{
	//The Epsilon (wrong function)
	//double dgs = tau*(exp(-(mdtime+Dt)/tafm)-exp(-(mdtime+Dt)/tafs));
	///Changed to Scale to G_MAX
	double dgs  = G_MAX*exp(-(mdtime+Dt)/(double)tafs);
	//Return the spike dynamics 
	return mfSynapticStrength*dgs;

}

double synapticTransmission::getSynapseStrength()
{
	return mfSynapticStrength;
}
 //Returns true if the spike has almost returned to 0
bool synapticTransmission::transmitExpired(void)
{
		return mbTransmitted;
}

synapticTransmission::~synapticTransmission(void)
{

}
