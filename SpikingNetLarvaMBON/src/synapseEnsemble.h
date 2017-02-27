#pragma once

class IFNeuron; //Fwd Declaration
class PoissonNeuron;
class INeuron;

class synapseEnsemble
{

public:
	synapseEnsemble(float simTimeStep,short sourceID,short ID=0);
	synapseEnsemble(float simTimeStep,int SynapsesNumber,float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float StartStrength,short sourceID=0,short ID=0,bool bNoPlasticity=false,unsigned short SourceFireRate=0);
	//Function called by afferent and by target neuron when a post synaptic spike occurs
	float SpikeArrived(synapse::SPIKE_SITE type); //Returns the avg strenght of the ensemble
	void StepNoSpike(); //Called when no spike event on the simulation step
	void StepNoSpike(double t);
	synapse* getSynapseArray();
	int getSynapsesCount();
	double getAvgStrength();
	unsigned short getSourceFireRate();
	void Reset(void);
	//When a Neuron registers this SynapseEnsemble it will register also pass
	// a pointer to its self so the SynapseEnsemble can notify the neuron of a spike arrival
	void RegisterNeuron(INeuron* pNeuron);

	short getsourceID();
	short gettargetID();
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
	synapse* arrSynapses;

	///Song Model
	void UpdateSongModFunctions(void);
	double gex_a,g_max; //Finite Conductance step when new spike arrives
	double Aspot,Asdep; //The Change amount used for learning
	double tafdep,tafpot;
	double songMdepot; //Current State of M(t). Learning M(t)-Decrease Strength when PreSynaptic Activity
	double songPot; //Current State of P(t) Increase strength after Neuron fire Event
};
