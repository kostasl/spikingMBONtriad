#pragma once

class synapticTransmission
{
public:
	synapticTransmission(float timestep,double SynapseStrength);
	double getCurrentStep(void); //Move Fwd in time by h and return I level at this time step
	double getValueAt(float Dt);//Evaluate gex at a later time Dt 
	double getSynapseStrength(); //Returns the strength of the synapse at the time of spike initiation
	bool transmitExpired(void); //Returns true if the spike has almost returned to 0
	~synapticTransmission(void);

private:
	float mfSimTimestep; //The simulation time step
	double mdtime; //Life time of current spike transmission
	float tafm; //Membrane time constant
	float tafs; //conductance time constant
	//double tau;
	bool mbTransmitted; //Flagged when the transmission has expired 6*tafs
	double mfSynapticStrength; //Holds the strength of the synapticEnsemble connection at the time of SpikeEvent
};
