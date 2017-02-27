/* ###################  SYNAPSE ENSEMBLE ##############################################
# The collection of single synapses. The Ensemble propagates spike events to synapses, 
# and holds the time since the last spike event. On every simulation time tick, either 
#  SpikeArrived or StepNoSpike() must be called.
# Kostas 23/6/07
######################################################################################
*/
#include "StdAfx.h"
#include "synapse.h"
#include "IFNeuron.h"
#include "PoissonNeuron.h"
#include "INeuron.h"
#include "synapticTransmission.h"
#include "synapseEnsemble.h"

//Constructor Not used
synapseEnsemble::synapseEnsemble(float simTimeStep=0.0001f,short sourceID=0,short ID)
{
	mfSimTimestep = simTimeStep;
	//Reset Time since last spike
	mdTimeSinceLastSpike = 0;
	mSourceFireRate = 0;
	msourceID=sourceID;
	mtargetID = 0;
	mpTargetNeuron		=0; //Init Pointer to 0

}


//Main Constructor
//If StartStrength is -ve then this converts the Synapse to an Inhibitory one
//By Default Plasticity is On
synapseEnsemble::synapseEnsemble(float simTimeStep,int SynapsesNumber,float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float StartStrength,short sourceID,short ID,bool bNoPlasticity,unsigned short SourceFireRate)
{
		mbNoPlasticity = bNoPlasticity; //True : Do not Modify mfAvgStrength
		msourceID=sourceID;
		mSourceFireRate  = SourceFireRate;
		mfStartStrength = StartStrength; //Holds The Constructor given value, Used when Reset

		#ifdef USE_SONG_LEARNING
		//Song Paper States Start Strength Maximum is g_max ~0.02
			if (StartStrength < 0) StartStrength = -(float)G_INH;//Inhibitory (will be multip by -1 later in IFN)
			else StartStrength = (float)G_MAX;
			//mfAvgStrength = StartStrength; 
			//SONG Learning
			gex_a		= G_MAX; //Start From Max Conductance
			Aspot		= 0.0050; //0.009; //Finite Conductance step when new spike arrives
			Asdep		= Aspot*1.06; //0.00945;//The Change amount used for learning A_-/A+=1.05
			songMdepot	= 0.0; //Current State of M(t). Learning M(t)-Decrease Strength when PreSynaptic Activity
			songPot		= 0.0; //Current State of P(t) Increase strength after Neuron fire Event
			tafdep		=0.020;//0.020f; F(Dt) Timeconst for Learned depression to decay tau_+ 
			tafpot		=0.020;//0.020f; Timeconst for Learned potentiation to decay tau_-
		#endif


		//Fast way to init memory
		arrSynapses = new synapse[SynapsesNumber];
        //Call Constructor
        for (int i=0;i<SynapsesNumber;i++)
            arrSynapses[i] = *(new synapse(A1,A2,tafPOT,tafDEP,nPOT,nDEP,(StartStrength),mbNoPlasticity));

		//Save Simulation timestep
		mfSimTimestep = simTimeStep;
		miSynapsesCount = SynapsesNumber;

		//Reset Time since last spike
		mdTimeSinceLastSpike = 0;
		mpTargetNeuron		=0; //Init Pointer to 0
		mfAvgStrength = StartStrength*SynapsesNumber;

}

 //Called when no spike event on the simulation step
void synapseEnsemble::StepNoSpike()
{
//Decay the Synaptic Modification Functions
#ifdef USE_SONG_LEARNING
	UpdateSongModFunctions();
#endif
	///Increment Time since last spike
	mdTimeSinceLastSpike+=mfSimTimestep;
}

///Used for Spike to Spike Interaction, Parameter is time since last spike
void synapseEnsemble::StepNoSpike(double t)
{
//Decay the Synaptic Modification Functions
#ifdef USE_SONG_LEARNING
	UpdateSongModFunctions();
#endif
	///Add Time since last spike
	mdTimeSinceLastSpike+=t;
}

void synapseEnsemble::UpdateSongModFunctions()
{
	songMdepot	= songMdepot*(double)exp(-mfSimTimestep/tafdep);
	songPot		= songPot*(double)exp(-mfSimTimestep/tafpot);
}

void synapseEnsemble::Reset()
{
	if (!arrSynapses) return;

    for(int i=0;i<miSynapsesCount;i++)
        arrSynapses[i].Reset();

	//mfAvgStrength = G_MAX; //Reset to Max COnductance
	mfAvgStrength = mfStartStrength*miSynapsesCount; //0; //Start from 0 for the Poisson Train test.
	mdTimeSinceLastSpike = 0;
}

///Propage the spike to all synapse members, reset the time since last spike
//Return the Change in Ds in strength
float synapseEnsemble::SpikeArrived(synapse::SPIKE_SITE type)
{
        double totalWeight=0;
                
        
		
		///Only obtain the Average When using Switch Rule
#ifndef USE_SONG_LEARNING
        //NOTE: Probably need to ignore PRE spikes when in State POT
        //And POST spikes when in DEP
		mfAvgStrength = 0;
        for(int i=0;i<miSynapsesCount;i++)
		{
			//Spike Arrived Returns the Change in Weight. Need to Sum to Obtain Aggregate
            totalWeight += arrSynapses[i].SpikeArrived(mdTimeSinceLastSpike,type);
			mfAvgStrength+=arrSynapses[i].getStrength();
		}
		//If PLasticity is Off The strength is not Saved in Synapses So AvgStrength=0;
		//Save & Return Average Strength
		//if (!mbNoPlasticity) ///If running in Plasticity
		//mfAvgStrength = totalWeight; // No Averaging  ///miSynapsesCount;
#endif
		//Reset Time since last spike
		mdTimeSinceLastSpike = 0;

		//Propagate Spike to Neuron if Target exists
		if (mpTargetNeuron && type == synapse::SPIKE_SITE::SPIKE_PRE)
		{
			synapticTransmission* sp = new synapticTransmission((float)mfSimTimestep,mfAvgStrength);
			if (mpTargetNeuron != 0) mpTargetNeuron->SpikeArrived(sp); //Spike Is passed to Neuron
		}

		#ifdef USE_SONG_LEARNING
				

		if (type == synapse::SPIKE_SITE::SPIKE_PRE)
		{ //Presynaptic Event Update P
			songPot += Aspot;
			gex_a = gex_a + songMdepot*G_MAX; //Decrease Synaptic Strength

			if (gex_a <0) gex_a = 0;

		}
		else//Postsynaptic Event
		{
			songMdepot -= Asdep; //Increase potential Decrease (M(t)in synaptic Strength
			gex_a = gex_a + songPot*G_MAX; //Increase Synaptic Strength

			if (gex_a > G_MAX) 
				gex_a = G_MAX; //IMpose limits 0>gex<gmax
		}

		if (!mbNoPlasticity) 
			mfAvgStrength = gex_a;

		//Do time Step Decay Modification of Learning State
		UpdateSongModFunctions(); //Testing if order significant

		#endif

		return (float)totalWeight;  //Return Change in s //OLD:Changed from AvgStrength (Always return New strength)
}
double synapseEnsemble::getAvgStrength()
{
	return mfAvgStrength;
}

int synapseEnsemble::getSynapsesCount()
{
	return miSynapsesCount;
}

//When a Neuron registers this SynapseEnsemble it will register also pass
// a pointer to its self so the SynapseEnsemble can notify the neuron of a spike arrival
//Called By TargetNeuron
void synapseEnsemble::RegisterNeuron(INeuron* pTargetN)
{
	mpTargetNeuron = pTargetN; //Synapse now connected to post Neuron	
	mtargetID = pTargetN->getID();
}

short synapseEnsemble::getsourceID()
{
	return msourceID;
}
short synapseEnsemble::gettargetID()
{
	return mtargetID;
}
//Logs Synapse State to open File stream
void synapseEnsemble::logtofile(std::ofstream &fs)
{
	fs << mID << " " << msourceID << " " <<  mtargetID << " " << mfAvgStrength << std::endl;
}

unsigned short synapseEnsemble::getSourceFireRate()
{
	return mSourceFireRate;
}

synapseEnsemble::~synapseEnsemble(void)
{
	delete [] arrSynapses;
}
