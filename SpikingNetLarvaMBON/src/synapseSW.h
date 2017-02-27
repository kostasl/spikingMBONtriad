#pragma once

class synapse
{
public:
	enum SW_STATE {POT, OFF, DEP};
	enum SPIKE_SITE {SPIKE_PRE, SPIKE_POST};

	synapse(void);
	synapse(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity);
	SW_STATE getState();
	float getStrength(); //Returns the Potentiation value
	float SpikeArrived(double t,SPIKE_SITE type); //throws Exception
	bool StochasticReturnedtoOFF(double t);
	double factorial(int k);///Calc factoria
	void Reset(); //Reset Strength and State
	~synapse(void);
    
	
	time_t t; //Used for random num generation

private:
	gsl_rng * rng_r; //Used by GSL Rand Num Generator

	float APOT,ADEP; //Step Change of Strength
    float tauPOT,tauDEP; //Time window for +/- in seconds
    int channelsPOT,channelsDEP;
	float fStrength,mfSreset;
	bool mbNoPlasticity; //When true the change in Strength is not saved
    SW_STATE mstate;
	


};
