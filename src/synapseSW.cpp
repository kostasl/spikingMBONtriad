/*###############	SYNAPSE	######
## Implements the individual Synaptic switch as described in the Appleby & Elliot 2005 STPD paper 
## The switch maintains state by holding the POT LEVEL (Strength) which is 0 > POT 
#################################
*/
#include "stdafx.h"
#include "isynapse.h"
#include "synapseSW.h"

#include "math.h"


//Default constructor
synapseSW::synapseSW(void)
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


synapseSW::synapseSW(const synapseSW& obj): ISynapse(obj)
{
    ///Setup the GSL Random number Generator
     // select random number generator
     rng_r = gsl_rng_alloc (gsl_rng_mt19937);
     if (!rng_r) throw "GSL RNG Init Failed. Out Of Memory?";

//	 unsigned int seed = unsigned(time(&t)) + rand()*100;
     unsigned int seed = rand()*100;

     gsl_rng_set(rng_r,seed);

  mstate        = obj.mstate;
  APOT          = obj.APOT;
  ADEP          = obj.ADEP;
  tauPOT        = obj.tauPOT;
  tauDEP        = obj.tauDEP;
  channelsPOT   = obj.channelsPOT;
  channelsDEP   = obj.channelsDEP;
  mfSreset      = obj.mfSreset;
  fStrength     = obj.fStrength;
  mbNoPlasticity = obj.mbNoPlasticity;

}


/// \brief Without this  derived constructor we get : undefined reference to `typeinfo for ISynapse
synapseSW::synapseSW(const ISynapse& obj):synapseSW((synapseSW&)obj)
{

}


synapseSW::synapseSW(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity)
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


synapseSW::SW_STATE synapseSW::getState()
    {
        
        return mstate;
    }

float synapseSW::getStrength()
{
	return fStrength;
}

void synapseSW::Reset()
{
	///NOTE CHANGE TO OFF
	mstate		= OFF;//Start in the OFF state
	fStrength = mfSreset;
}

float synapseSW::SwitchRulePlasticity(double t,SPIKE_SITE type)
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

    return Ds;
}

/// \brief Handles Spike events of either post or presynaptic origin - Passes these to the switch rule
/// \param time of event
/// \param Enum value of where the Spike occured, presynaptic or Postsynaptic
/// \return the Change Ds due to spike Arrival.
float synapseSW::SpikeArrived(double t,SPIKE_SITE type)
{
        float Ds=0;//=fStrength;// Potentiation -> This is now in SwitchPlasticity Rule

		//Ds=(Ds<0)?0:Ds; 
        if (mbNoPlasticity)
		{
            //Plasticity is off
			//Low Bound is 0 
            Ds = 0; //Save No Change
		}
        else //Run Switch-Rule Plasticity
        {
            Ds += SwitchRulePlasticity(t,type);
            fStrength += Ds;
        }

		//return the change in strength
		return Ds;

}

    //Parameter t: Time since last spike
bool synapseSW::StochasticReturnedtoOFF(double t)
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

double synapseSW::factorial(int k)
{
	if (k==0) return 1;
	for (int j=(k-1);j>0;j--) k*=j;
    return k;
}



synapseSW::~synapseSW(void)
{
	//Free resources Used by GSL RNG
    //try{
    gsl_rng_free(rng_r);
    //}catch (int e)
    //{

    //}
}
