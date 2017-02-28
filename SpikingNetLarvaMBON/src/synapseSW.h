#pragma once

#include "isynapse.h"

class synapseSW: public ISynapse
{
public:
	enum SW_STATE {POT, OFF, DEP};

    synapseSW(void);
    /// \brief Copy constructor
    synapseSW(const synapseSW& obj);
    synapseSW(const ISynapse& obj);
    synapseSW(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity);
    virtual SW_STATE getState();
    virtual float getStrength(); //Returns the Potentiation value
    virtual float SpikeArrived(double t,SPIKE_SITE type); //throws Exception
    float SwitchRulePlasticity(double t,SPIKE_SITE type); //implemets Switch Plasticity Model
	bool StochasticReturnedtoOFF(double t);
	double factorial(int k);///Calc factoria
    virtual void Reset(); //Reset Strength and State
    ~synapseSW(void);
    
	
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
