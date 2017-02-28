#pragma once
#include "stdafx.h"
#include "math.h"
#include "isynapse.h"
#include "synapticTransmission.h"
#include "INeuron.h"

class IFNeuron; //Fwd Declaration
class PoissonNeuron;
class INeuron;

template <class T, int N>
class synapseEnsemble
{

public:
    synapseEnsemble(float simTimeStep=0.0001f,short sourceID=0,short ID=0);
    //synapseEnsemble(float simTimeStep,int SynapsesNumber,float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float StartStrength,short sourceID=0,short ID=0,bool bNoPlasticity=true,unsigned short SourceFireRate=0);
    synapseEnsemble(float simTimeStep,ISynapse& osynapse); //This uses the synapse Instance To create copies which are held in the ensemble
	//Function called by afferent and by target neuron when a post synaptic spike occurs
    virtual float SpikeArrived(ISynapse::SPIKE_SITE type); //Returns the avg strenght of the ensemble
    virtual void StepNoSpike(); //Called when no spike event on the simulation step
    virtual void StepNoSpike(double t);
    virtual ISynapse* getSynapseArray();
    virtual int getSynapsesCount();
    virtual double getAvgStrength();
    virtual unsigned short getSourceFireRate();
    virtual void Reset(void);
	//When a Neuron registers this SynapseEnsemble it will register also pass
	// a pointer to its self so the SynapseEnsemble can notify the neuron of a spike arrival
    virtual void RegisterNeuron(INeuron* pNeuron);
    virtual short getsourceID();
    virtual short gettargetID();
	void logtofile(std::ofstream &fs);//Passed an open filestream it will log the synapse strength
	~synapseEnsemble(void);


private:
	short mID; //An ID externaly created to identify the synapse in a population synapsing on the same Neuron
	short msourceID; //An Id for the source of this synapse - Target ID is taken from RegisterTarget
	short mtargetID; 
	double mfSimTimestep;//The simulation time step in seconds
	double mdTimeSinceLastSpike; //In Seconds
	double mfAvgStrength; //Was Float made into Double
	unsigned short mSourceFireRate; //Holds the Afferent Fire Rate, used by PoissonNeuron
	int	  miSynapsesCount;
	float mfStartStrength; //Holds The Constructor given value, Used when Reset

	bool mbNoPlasticity;
	INeuron* mpTargetNeuron;
    T arrSynapses[N]; //Initialize memory Array

	///Song Model
	void UpdateSongModFunctions(void);
	double gex_a,g_max; //Finite Conductance step when new spike arrives
	double Aspot,Asdep; //The Change amount used for learning
	double tafdep,tafpot;
	double songMdepot; //Current State of M(t). Learning M(t)-Decrease Strength when PreSynaptic Activity
	double songPot; //Current State of P(t) Increase strength after Neuron fire Event
};


/// \brief Main Constructor. If StartStrength is -ve then this converts the Synapse to an Inhibitory one
/// \note By Default Plasticity is On
template<class T,int N>
synapseEnsemble<T,N>::synapseEnsemble(float simTimeStep, ISynapse& osynapse)
{
       // mbNoPlasticity      = bNoPlasticity; //True : Do not Modify mfAvgStrength
        //msourceID           = sourceID;
        //mSourceFireRate     = SourceFireRate;
        mfStartStrength      = osynapse.getStrength(); //Holds The Constructor given value, Used when Reset

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


        //Already Initialized T x N- Fast way to init memory
        //arrSynapses = new synapse[SynapsesNumber];
        //Call Constructor
        for (int i=0;i<N;i++)
            arrSynapses[i] = *(new T(static_cast<T>(osynapse)));

        //Save Simulation timestep
        mfSimTimestep           = simTimeStep;
        miSynapsesCount         = N;

        //Reset Time since last spike
        mdTimeSinceLastSpike    = 0;
        mpTargetNeuron          = 0; //Init Pointer to 0
        mfAvgStrength           = mfStartStrength*N;

}



//Called when no spike event on the simulation step
template<class T,int N>
void synapseEnsemble<T,N>::StepNoSpike()
{
//Decay the Synaptic Modification Functions
#ifdef USE_SONG_LEARNING
   UpdateSongModFunctions();
#endif
   ///Increment Time since last spike
   mdTimeSinceLastSpike+=mfSimTimestep;
}

///Used for Spike to Spike Interaction, Parameter is time since last spike
template<class T,int N>
void synapseEnsemble<T,N>::StepNoSpike(double t)
{
//Decay the Synaptic Modification Functions
#ifdef USE_SONG_LEARNING
   UpdateSongModFunctions();
#endif
   ///Add Time since last spike
   mdTimeSinceLastSpike+=t;
}

/// \brief For Song Plasticity Specific
/// \todo remove it and place in SongSynapse Object
template<class T,int N>
void synapseEnsemble<T,N>::UpdateSongModFunctions()
{
   songMdepot	= songMdepot*(double)exp(-mfSimTimestep/tafdep);
   songPot		= songPot*(double)exp(-mfSimTimestep/tafpot);
}


template<class T,int N>
void synapseEnsemble<T,N>::Reset()
{
   if (!arrSynapses) return;

   for(int i=0;i<miSynapsesCount;i++)
       arrSynapses[i].Reset();

   //mfAvgStrength = G_MAX; //Reset to Max COnductance
   mfAvgStrength = mfStartStrength*miSynapsesCount; //0; //Start from 0 for the Poisson Train test.
   mdTimeSinceLastSpike = 0;
}

