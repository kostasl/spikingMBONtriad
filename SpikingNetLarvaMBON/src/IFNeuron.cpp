/// \class IFNeuron  Integrate and Fire Neuron ###
///\brief Implements a simple 1st order Euler Numerical approximation solution to the Diff Eqn
/// dV/dt=1/taf*(El-V(t)+R*Sumj[sj*Sum_tfr(epsilon(t-tf))]))
/// Details: Everytime a spike event occurs a synapticTransmission object is created and added to a list
/// in the target IFNeuron. The object is removed from the list when it's life time expires.
/// IF Neuron shows a spike potential of 20mv, for 1 timestep
///
/// \note EXTENSION: The Song et al. 2000 learning and conductance u\pdate method is also accommodated. The class can use
/// either the independent decay of conductances using a stack of synaptic response (SynapticTransmission) objects
/// or by Defining USE_SONG_CONDUCTANCE in the header file, all conductances are summed and decay is treated
/// as a single exponential.
/// \todo: When Using the Spike Stack need to extend to Have a seperate Stack for Inhibitory Injections

#include "stdafx.h"
#include <gsl/gsl_math.h>
#include "synapseSW.h"
#include "synapseEnsemble.h"
#include "synapticTransmission.h"
#include "INeuron.h"
#include "IFNeuron.h"

IFNeuron::IFNeuron(float timestep,short ID)
{
	mID					= ID;
    iLastEfferentSynapseIndex	= 0;
    iLastSpikeIndex		= 0.0;
    Eex					= 0.0; //mV //Excitatory Reversal Potential
    Ein					= -0.070;//mV //inhibitory Reversal Potential
	Vrest				= -0.070;//mV
    Vspike              = 0.020;//Brief t=1 Spike Potential (Delta like)
#ifdef USE_SONG_LEARNING
	Vreset				= -0.060; //mV -0.060 Proposed by Song, should use -90mV
#else
	Vreset				= -0.090; //mV -0.060 Proposed by Song, should use -90mV
#endif
	Vm					= Vrest;//Membrane Potential -Start From
	Vthres				= -0.054; //54mV in Song Paper
	tafm				= 0.020; //ms Membrane Time Const
	tafs				= 0.005f; //ms of Conductance channel (tau_ex and tau_in on Paper)
	h					= timestep;	

	//SONG Model - Conductance Variables
    gex_song			= 0.0;//Start Value
    gin_song			= 0.0;
	gmax				= G_MAX;		
	//gex_a				= 0.02; 
	
	//Statistics

    uiNumberOfSpikesInPeriod = 0;
    uiPeriodOfSpikeCount	= 0.0f;
    fMeanFireRate			= 0.0f;

#ifdef DEBUG_LOG
	char *File = new char[50];
	char *buff = new char[4];
	strcpy(File,"outputIFN_");
	strcat(File,_itoa(mID,buff,10));//Convert Id to String
	strcat(File,".csv");
	out.open(File,std::ios::out);
	out << "#ID,Pre,Post,Vm,gex,stack" << std::endl;

#endif

}

/// \brief Registers an incoming pathway (afferent) synapse that provides impulses to this neuron
/// notifies the synapse ensemble that this neuron is its target
void IFNeuron::RegisterAfferent(ISynapseEnsemble* pSyn)
{
	//Check that the array is not Maxed.
    if (iLastAfferentSynapseIndex == (MAX_AFFERENTS-1))
	{
		std::cout << "Max afferent number reached, Increase MAX_AFFERENTS ";
		abort();
		throw "Max afferent number reached, Increase MAX_AFFERENTS ";
		//Can have code to increase array size with malloc()
	}

    if (!pSyn)
    {
        std::cerr << "missing synapse ensemble pointer!" << std::endl;
        abort();
    }
	//Register Back to the Synapse Ensemble
    pSyn->RegisterAfferentNeuron(this);
	//Add to private array
    mAfferents[iLastAfferentSynapseIndex] = pSyn;
    iLastAfferentSynapseIndex++;
}

/// \name RegisterEfferent
/// \brief Save a pointer to an outgoing synapse  It connects to an efferent (target) neuron via the ensemble Synapse to this neuron
/// \param A synapse Ensemble/set connecing the Afferent neuron to this neuron
void IFNeuron::RegisterEfferent(ISynapseEnsemble* pSyn)
{
    //Check that the array is not Maxed.
    if (iLastEfferentSynapseIndex == (MAX_AFFERENTS-1))
    {
        std::cout << "Max efferent number reached, Increase MAX_AFFERENTS ";
        abort();
        throw "Max efferent number reached, Increase MAX_AFFERENTS ";
        /// \todo Can have code to increase array size with malloc() or convert to std vector

    }

    //Register Back to the Synapse Ensemble
    pSyn->RegisterEfferentNeuron(this); //The Receiving Neuron
    //Add to private array
    mEfferents[iLastEfferentSynapseIndex] = pSyn;
    iLastEfferentSynapseIndex++;

}


 //Called by SynapseEnsemble
//Adds the Spike to the list of current transmitting spikes

void IFNeuron::SpikeArrived(synapticTransmission* Spike)
{
	bPrespiketoLog = true;

	//SOng combined with Elliot.
	#ifdef USE_SONG_CONDUCTANCE
		#ifndef USE_SONG_LEARNING 
			//Using Switch Rule
        if (Spike->getSynapseStrength() > 0.0)
            gex_song += Spike->getSynapseStrength()*G_MAX; //Add conductance step
        else
            gin_song += Spike->getSynapseStrength()*G_INH; //Add conductance step

		#else
			//Using Song Learning
			double fInjection;
			
			fInjection = Spike->getSynapseStrength(); 
			//Check if Inhibitory Or Ex
			if (fInjection > 0)
				gex_song += fInjection; //Will use the M(t),P(t) update rule
			else
                gin_song -= fInjection; //Substract cause Injection is -ve
				//SONG Learning will be found in the SynapseEnsemble


		#endif
			//Clear Memory from syntransmission object
			delete Spike;
		return; //Dont execute further - increase speed
	#endif

	PushSpikeStack();
	mSpikes[0] = Spike;	
	
	#ifdef VERBOSE
		std::cout << std::endl << "-* Neuron Spike Arrived SpikeStack Length:" << iLastSpikeIndex << std::endl ;
	#endif
}

///Copies each spike to the next item in the array effectively pushing the stack
void IFNeuron::PushSpikeStack()
{
	if (iLastSpikeIndex == (MAX_SPIKES-1)) 
	{
		std::cout << "Error: Max Spike Stack Size reached, Can't Push!";
//		throw "Max Spike Stack Size reached, Can't Push!";
		abort();
	}
		//if (iLastSpikeIndex ==0) return; //Nothing to Do, Empty Stack
	
	//Start from end and copy each spike one up
	for (int i=(iLastSpikeIndex-1);i>=0;i--)
	{
		if (!mSpikes[i]) 
		{
			std::cout << "Error: PushSpikeStack Failed: Null Spike found";
			abort();
			throw "PushSpikeStack Failed: Null Spike found";
		}
			mSpikes[i+1]=mSpikes[i]; //Copy to next
	}
	///Now pointer at 0 and 1 point to the same object.[0] is about to be replaced by new Spike
	iLastSpikeIndex++; //Increment the last item index to reflect the change

}

//Check each spike, and delete the Expired ones. Expired ones should be found at the end of the array
//Should be called once after each time step
short IFNeuron::ClearSpikeList()
{
	short cntDeleted=0;
	//Works because the expired ones are always at one end of the array
	//this is taken care of by the SpikeArrived, which creates a stack of Spikes
	for (int i=(iLastSpikeIndex-1);i>=0;i--)
	{
		if (!mSpikes[i])
			break; //Stop if a null found

		if(mSpikes[i]->transmitExpired())
		{
			cntDeleted++;
			iLastSpikeIndex--;
			delete mSpikes[i];
		}
	}
return cntDeleted;
}
//Run Through All spikes in stack and get the step contribution to the V
double IFNeuron::SumExInjections()
{
#ifdef USE_SONG_CONDUCTANCE
	double sgex = gex_song; //Save state 
	//Update State for next simulation step
	gex_song = gex_song*exp(-h/tafs);
	return sgex; //Return Current State after updating for next time step
#endif

	double Vex=0;
	for(int i=0;i<iLastSpikeIndex;i++)
	{
		if (!mSpikes[i]) break; //Null so end of Stack
		Vex+=mSpikes[i]->getCurrentStep();
	}
	return Vex;
}

//Run Through All spikes in stack and get the step contribution to the V at a future time Dt
//Note: The time is not stepped fwd for the Spikes.
//Sums All independent of Inhibitory Or Exhitatory,
double IFNeuron::SumExInjections(float Dt)
{
#ifdef USE_SONG_CONDUCTANCE	
	 
	return gex_song*exp((double)-Dt/(double)tafs);
#endif

	double Vex=0;
	//Sums All independent of Inhibitory Or Exhitatory, No Prediction For Inhibitory On Stack Conductance
	for(int i=0;i<iLastSpikeIndex;i++)
	{
		if (!mSpikes[i]) break; //Null so end of Stack
		Vex+=mSpikes[i]->getValueAt(Dt);
	}
	return Vex;
}

//Inhibitory Injection on Current Time and update the state for the next simulation step
double IFNeuron::SumInhInjections()
{
#ifdef USE_SONG_CONDUCTANCE
	double sgin = gin_song; //Save state 
	//Update State for next simulation step
	gin_song = gin_song*exp((double)-h/(double)tafs);
	return sgin; //Return Current State after updating for next time step
#endif

	return 0;
}

//Inhibitoty Injection Decay in Future Time Dt
double IFNeuron::SumInhInjections(float Dt)
{
#ifdef USE_SONG_CONDUCTANCE	
	 
    return gin_song*exp((double)-Dt/(double)tafs);
#endif
    ///\todo Missing inh Spikes Summation
	return 0;
}

void IFNeuron::StepSimTime()
{
    StepRK_UpdateVm();
}

///\todo I noticed some Overflows in Vm
//Solve the Diff Equation , Step Time Fwd Returns the new Membrane V
double IFNeuron::StepRK_UpdateVm(void)
{
	double gex,gin,k1,k2,k3,k4;
    double dV=0.0;
	bPostspiketoLog = false;

	ClearSpikeList();
    if (Vm > Vthres) //If Spike Occured on previous step - Then set back to Reset potential
        Vm=Vreset;
	
	///DO the RK Update Algorithm
	//1 k1=hf(xn,yn)
	gex = SumExInjections(0); //Spikes Don't Step fwd in time
	gin = SumInhInjections(0);
    dV=(1/tafm)*(Vrest-Vm+gex*(Eex-Vm)+gin*(Ein-Vm));
	k1 = h*dV;

	//2 k2=hf(xn+0.5h,yn+0.5k1)
	gex = SumExInjections((float)(0.5*h));
	gin = SumInhInjections((float)(0.5*h));
    dV=(1/tafm)*(Vrest-(Vm+0.5*k1)+gex*(Eex-Vm+0.5*k1) + gin*(Ein-Vm+0.5*k1));
	k2 = h*(dV);

	//3 k3=hf(xn+0.5h,yn+0.5k2)
    dV=(1/tafm)*(Vrest-(Vm+0.5*k2)+gex*(Eex-Vm+0.5*k2) + gin*(Ein-Vm+0.5*k2));
	k3= h*dV;

	//4 k4=hf(xn+h,yn+k3)
	gex = SumExInjections(); //Step Spike time Fwd by h /Like Calling SumExInjections(h)
	gin = SumInhInjections();
    dV=(1/tafm)*(Vrest-(Vm+k3)+gex*(Eex-Vm+k3) + gin*(Ein-Vm+k3));
	k4 = h*dV;

	//increment Time -> No need, done by spike time increment
	//Update yn -> yn+1
	Vm+=(1/6.0f)*(k1+2*k2+2*k3+k4); //Update the membrane Potential

	///END OF RUNGE-KUTTA - Process repeats at every time step h
	//Check if threshold reached
	
	#ifdef VERBOSE
		std::cout << "Vm: " << Vm <<std::endl;
	#endif

	///Log OUtput
		#ifdef DEBUG_LOG
			out << mID << " "; //Print this neurons ID to distinguish in the log
		   if (bPrespiketoLog) out << "1"; else out << "0"; //PreSpike
		   bPrespiketoLog = false; //Reset until next pre Spike
			
		   if (Vm>Vthres) out << " 1"; else out << " 0"; //Post spike
		   out << " " << Vm << " " << gex << " " << iLastSpikeIndex << std::endl;
		#endif

        //COunt POst Rate Timer
         uiPeriodOfSpikeCount++;

        //Mean Period Elapsed? Reset Period

         //Exponential Decay Of Running Average //Bin # Spike
         if (uiPeriodOfSpikeCount > IFFIRERATE_PERIOD)
         {
            uiPeriodOfSpikeCount = 0;
            uiNumberOfSpikesInPeriod = 0;

        }

	///Threshold Reached - Fire Action Potential
	if (Vm>Vthres) 
	{
        uiNumberOfSpikesInPeriod++;
        //Vm=Vreset;
        Vm =Vspike;
		ActionPotentialEvent(); //Send Post Syn Spike event back to Synapses
	
        //return Vthres;
	}
	

	return Vm;
}

double IFNeuron::getMembraneVoltage()
{
    return Vm;
}

///Sends the post Synaptic event Back to the Registered Synapses
void IFNeuron::ActionPotentialEvent()
{

	bPostspiketoLog = true;

    //double nRate= 0; //Loop until the largest index is passed between the Afferent / Efferent Synapse Lists
    for(unsigned int i=0;i< GSL_MAX(iLastEfferentSynapseIndex,iLastAfferentSynapseIndex);i++)
    {
        if (!mAfferents[i] &&  !mEfferents[i]) break; //Null so Next - So Skip Loop

        //Notify all synapses (Time steps Fwd?)
        if (mAfferents[i]) //Notify Input pathways
            mAfferents[i]->SpikeArrived(ISynapse::SPIKE_POST); //incoming synapses receive A spike (BAP) at the post synaptic site

        if (mEfferents[i]) //Notify Output Pathway
            mEfferents[i]->SpikeArrived(ISynapse::SPIKE_PRE); //OutGoing Synapses Receive A pre Spike

        ///Deprecated :nRate+= mSynapses[i]->getSourceFireRate()*mSynapses[i]->getAvgStrength();
    }


#ifdef VERBOSE
	std::cout << std::endl << "-* Neuron POST Spike ";
#endif

}

int IFNeuron::getID()
{
	return mID;
}

float IFNeuron::getFireRate()
{
    //msLastFireRate = msNumberOfSpikesInPeriod/mfPeriodOfSpikeCount;

    //msNumberOfSpikesInPeriod=0; //Reset Counter
    //mfPeriodOfSpikeCount = 0.0;
    double alpha = 1/IFFIRERATE_PERIOD;
    fMeanFireRate = (alpha*uiNumberOfSpikesInPeriod) + (1.0 - alpha) * fMeanFireRate;


    return fMeanFireRate;
}

//Returns true if neuron has fired during previous Step
bool IFNeuron::ActionPotentialOccured()
{
	return bPostspiketoLog;
}
IFNeuron::~IFNeuron(void)
{
#ifdef DEBUG_LOG
	out.close(); //Close File Handle
#endif

	//Clear Memory - Delete All Spikes
	for(int i=0;i<iLastSpikeIndex;i++)
	{
		if (!mSpikes[i]) continue; //Null so Next
		delete mSpikes[i];
	}

//Delete SynapseEnsembles //Already Deleted
    for(unsigned int i=0;i<iLastEfferentSynapseIndex;i++)
	{
 //       if (!mAfferents[i]) continue; //Null so Next
//        delete mAfferents[i];
        mAfferents[i] = 0;
	}

}
	//Reset Neuron To start conditions
	void IFNeuron::Reset()
	{
		ClearSpikeList();
		Vm	= Vrest;
		gex_song = 0;
		gin_song = 0;
	}
