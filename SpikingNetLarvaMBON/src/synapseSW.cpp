/*###############	SYNAPSE	######
## Implements the individual Synaptic switch as described in the Appleby & Elliot 2005 STPD paper 
## The switch maintains state by holding the POT LEVEL (Strength) which is 0 > POT 
#################################
*/
#include "StdAfx.h"
#include "synapse.h"


//Default constructor
synapse::synapse(void)
{
	//Set default model parameters
	mstate		= OFF;//Start in the OFF state
    APOT        = 1;
    ADEP        = 0.95f;
    tauPOT      = 0.020f;
    tauDEP      = 0.013f;
    channelsPOT = 3;
    channelsDEP = 3;
	mfSreset	= 0.1f;
	fStrength	= 0.1f;

	///Setup the GSL Random number Generator
	 // select random number generator 
	 rng_r = gsl_rng_alloc (gsl_rng_mt19937);
	 if (!rng_r) throw "GSL RNG Init Failed. Out Of Memory?";

//	 unsigned int seed = unsigned(time(&t)) + rand()*100;
	 unsigned int seed = rand()*100;

	 gsl_rng_set(rng_r,seed);

	 //Seed random number generator
	 srand(seed); 
		
}

synapse::synapse(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity)
{
	///NOTE CHANGE TO OFF
	mstate		= OFF;//Start in the OFF state
    APOT        = A1;
    ADEP        = A2;
    tauPOT      = tafPOT;
    tauDEP      = tafDEP;
    channelsPOT = nPOT;
    channelsDEP = nDEP;

	mfSreset	= Sreset;
	fStrength	= Sreset; //The Start and reset Strength

	mbNoPlasticity = bNoPlasticity; //No Change in weights Occurs

 	///Setup the GSL Random number Generator
	 // select random number generator 
	 rng_r = gsl_rng_alloc (gsl_rng_mt19937);
	 if (!rng_r) throw "GSL RNG Init Failed. Out Of Memory?";

	 unsigned int seed = unsigned(time(&t)) + rand()*100;
	 //unsigned int seed = rand()*100;

	 gsl_rng_set(rng_r,seed);


	//Seed random number generator
	srand(seed); 
}


synapse::SW_STATE synapse::getState()
    {
        
        return mstate;
    }

float synapse::getStrength()
{
	return fStrength;
}

void synapse::Reset()
{
	///NOTE CHANGE TO OFF
	mstate		= OFF;//Start in the OFF state
	fStrength = mfSreset;
}

//T Time since last spike
//Returns the Change Ds due to spike Arrival.
float synapse::SpikeArrived(double t,SPIKE_SITE type)
{
        float Ds=0;//=fStrength;// Potentiation
        ///React according to current sw state
        switch (mstate)
        {
        case OFF:
            if (type==SPIKE_POST)
                //FROM OFF POst Spike arrived, switch to DEP
                mstate=DEP;
            else
                //FROM OFF a pre spike arrived switch on POT
                mstate=POT;
                ///PRE Synaptic Spike
                
            break;
            
         case POT:
                 if (type==SPIKE_POST)
                 {
                     //FROM POT a POSTSyn spike arrived 
                     //Has the switch closed?
                    if (StochasticReturnedtoOFF(t)) 
                    {
                        //Since Switch At off and POST received go to DEP
                        mstate=DEP;
                        ///Switch off, out of time
                    }
                    else
                    {//ON
                        //Switch operational In POT Mode, so potentiate
                        ///Normal POT
                        Ds+=APOT;//Potentiate A+
                        mstate=OFF; //Return to OFF
                    }   
                 }///PRE Synaptic Spike
                 else
                 {
                     ///Paper says ignore
                     //Ds=0;
                     mstate=POT;
                 }
             break;
             
         case DEP:
                 if (type==SPIKE_PRE)
                 {
                     //FROM POT a POSTSyn spike arrived 
                     //Has the switch closed?
                    if (StochasticReturnedtoOFF(t)) 
                    {
                        ///Switch was off and PRE spike arrived, Switch to POT
                        mstate=POT;
                        ///Switch off, out of time
                    }
                    else
                    {//ON -- DO Dep
                        //Switch operational In DEP Mode, so de-potentiate
                        ///Normal DEP
                        Ds-=ADEP;//Potentiate A-
                        mstate=OFF; //Return to OFF
                    }   
                 }///DEP AND POST Synaptic Spike
                 else
                 {
                     ///Paper says ignore
                     //Ds=0;
                     mstate=DEP;
                 }

             break;
                
          default:
                //throw new Exception("Invalid State");
            
			  throw "Synapse at Uknown State!";
              Ds=-1000;
                
        }//Close Switch
		//Ds=(Ds<0)?0:Ds; 
		if (!mbNoPlasticity)
		{
			//Plasticity is on
			//Low Bound is 0 
			fStrength += Ds; ///Save new Strength

			// Temp out if (fStrength<0) fStrength =0;
		}
		//return the change in strength
		return Ds;

}

    //Parameter t: Time since last spike
bool synapse::StochasticReturnedtoOFF(double t)
    {
        double taf;
		int	  n;
        double P;
        double en=0;


		if (mstate==OFF) return true;
        
        if (mstate==POT)
        {
            taf= tauPOT;
            n = channelsPOT;
        }
        else
        {
            taf= tauDEP;
            n = channelsDEP;            
        }
    

        //Calc en+/- Note: Paper states j=0..n-1
        for (int j=0;j<n;j++)
			en+=pow((t/taf),j)/factorial(j);
        //Calc Probability of switch returned to OFF
        P=1-exp((double)(-t/taf))*en;
        
        //switch is OFF if Random Num is drawn within OFF probability range
		double r = gsl_rng_uniform(rng_r);// OLD: rand()/(double)RAND_MAX;
        if (r<=P) 
			return true;
        else 
			return false;

}

double synapse::factorial(int k)
{
	if (k==0) return 1;
	for (int j=(k-1);j>0;j--) k*=j;
    return k;
}


synapse::~synapse(void)
{
	//Free resources Used by GSL RNG
	gsl_rng_free(rng_r);
}
